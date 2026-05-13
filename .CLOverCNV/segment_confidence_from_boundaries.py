#!/usr/bin/env python3
import argparse
import gzip
import math
from typing import Dict, Optional, Tuple

import numpy as np
import pandas as pd


def read_tsv_any(path: str, header="infer") -> pd.DataFrame:
    if path.endswith(".gz"):
        return pd.read_csv(path, sep="\t", header=header, compression="gzip")
    return pd.read_csv(path, sep="\t", header=header)


def read_segments_bed(path: str) -> pd.DataFrame:
    """
    segments.bed(.gz) like:
      chrom start end name score strand y ratio depth n_probes
    Extra columns ok.
    """
    df = read_tsv_any(path, header=None)
    if df.shape[1] < 4:
        raise RuntimeError("segments.bed(.gz) must have at least 4 columns: chrom start end name")

    df = df.rename(columns={0: "chrom", 1: "start", 2: "end", 3: "name"})
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)

    # Optional known columns (based on your file)
    if df.shape[1] > 4:
        df = df.rename(columns={4: "bed_score"})
    if df.shape[1] > 5:
        df = df.rename(columns={5: "strand"})
    if df.shape[1] > 6:
        df = df.rename(columns={6: "y"})
    if df.shape[1] > 7:
        df = df.rename(columns={7: "ratio"})
    if df.shape[1] > 8:
        df = df.rename(columns={8: "depth_mean"})
    if df.shape[1] > 9:
        df = df.rename(columns={9: "n_probes"})

    return df


def read_boundary_support(path: str) -> pd.DataFrame:
    """
    boundary_support.tsv(.gz) from your Step 3.
    Must have at least:
      chrom, pos, name, combined_z_chr, combined_support (support optional but useful)
    """
    df = read_tsv_any(path, header="infer")
    required = {"chrom", "pos", "name", "combined_z_chr"}
    missing = required - set(df.columns)
    if missing:
        raise RuntimeError(f"boundary_support TSV missing required columns: {sorted(missing)}")

    df["pos"] = df["pos"].astype(int)
    df["combined_z_chr"] = df["combined_z_chr"].astype(float)

    if "combined_support" in df.columns:
        df["combined_support"] = df["combined_support"].astype(float)

    return df


def robust_score_0_1000_from_z(z: float, z_max: float = 8.0) -> int:
    """
    IGV-like BED 'score' expects 0..1000.
    We clamp negative z to 0 because negative means 'below typical boundary evidence'.
    """
    if z is None or (isinstance(z, float) and (math.isnan(z) or math.isinf(z))):
        return 0
    zc = max(0.0, min(float(z), float(z_max)))
    return int(round(1000.0 * (zc / float(z_max))))


def robust_score_0_1000_from_raw(x: float, x_p99: float) -> int:
    """
    Map raw combined_support to 0..1000 using a chromosome/global percentile scale.
    """
    if x is None or (isinstance(x, float) and (math.isnan(x) or math.isinf(x))):
        return 0
    if x_p99 <= 0:
        return 0
    xc = max(0.0, min(float(x), float(x_p99)))
    return int(round(1000.0 * (xc / float(x_p99))))


def nearest_pos_value(pos_sorted: np.ndarray, val_by_pos: Dict[int, float], target: int, tol: int) -> Tuple[Optional[int], Optional[float]]:
    """
    Find nearest position within tol; return (pos, value) else (None, None).
    """
    if pos_sorted.size == 0:
        return (None, None)
    i = np.searchsorted(pos_sorted, target)
    candidates = []
    if i < pos_sorted.size:
        candidates.append(int(pos_sorted[i]))
    if i > 0:
        candidates.append(int(pos_sorted[i - 1]))

    best_pos = None
    best_dist = None
    for p in candidates:
        d = abs(p - target)
        if best_dist is None or d < best_dist:
            best_dist = d
            best_pos = p

    if best_pos is not None and best_dist is not None and best_dist <= tol:
        return (best_pos, float(val_by_pos[best_pos]))
    return (None, None)


def main():
    ap = argparse.ArgumentParser(
        description="Part 4: compute segment confidence from flanking boundary support (name-based pairing; fallback to nearest-pos matching)."
    )
    ap.add_argument("--segments-bed", required=True, help="segments.bed(.gz) from Step 2.")
    ap.add_argument("--boundary-support-tsv", required=True, help="boundary_support.tsv(.gz) from Step 3.")
    ap.add_argument("--out-tsv", required=True, help="Output TSV(.gz ok) with per-segment confidence.")
    ap.add_argument("--out-bed", default="", help="Optional BED(.gz ok) scored by confidence for IGV.")

    ap.add_argument("--method", choices=["min", "mean", "max"], default="mean",
                    help="Combine left/right boundary values into segment confidence. Default mean.")
    ap.add_argument("--missing-policy", choices=["one_side", "na"], default="one_side",
                    help="If one boundary missing (chrom end), use the other (one_side) or set NA.")
    ap.add_argument("--z-col", default="combined_z_chr", help="Boundary column to use for confidence. Default combined_z_chr.")
    ap.add_argument("--raw-col", default="combined_support", help="Boundary raw column for optional scoring. Default combined_support.")

    ap.add_argument("--fallback-pos-match", action="store_true",
                    help="If name-based boundary lookup fails, fall back to nearest position matching between segments.")
    ap.add_argument("--tol-bp", type=int, default=2000,
                    help="Tolerance for fallback position matching (bp). Default 2000 (safe for gapped segments).")

    ap.add_argument("--score-from", choices=["z", "raw"], default="z",
                    help="What to map into 0..1000 BED score. Default z (combined_z_chr; negatives -> 0).")
    ap.add_argument("--z-max", type=float, default=8.0,
                    help="Clamp z at this value for BED scoring. Default 8.")
    ap.add_argument("--raw-p99", type=float, default=0.0,
                    help="Override p99 scaling for raw score mapping. If 0, compute from boundary table globally.")

    args = ap.parse_args()

    seg = read_segments_bed(args.segments_bed)
    b = read_boundary_support(args.boundary_support_tsv)

    if args.z_col not in b.columns:
        raise RuntimeError(f"--z-col '{args.z_col}' not found in boundary support columns: {list(b.columns)}")

    # Name-based lookup: (chrom, name) -> boundary value
    # If duplicates exist, keep max (conservative for "support")
    b_name = (
        b.groupby(["chrom", "name"], as_index=False)[args.z_col]
         .max()
         .rename(columns={args.z_col: "bval_z"})
    )
    bval_by_chr_name: Dict[Tuple[str, str], float] = {
        (r["chrom"], r["name"]): float(r["bval_z"]) for _, r in b_name.iterrows()
    }

    # Also keep raw lookup if present
    bval_raw_by_chr_name: Dict[Tuple[str, str], float] = {}
    raw_available = args.raw_col in b.columns
    if raw_available:
        b_name_raw = (
            b.groupby(["chrom", "name"], as_index=False)[args.raw_col]
             .max()
             .rename(columns={args.raw_col: "bval_raw"})
        )
        bval_raw_by_chr_name = {
            (r["chrom"], r["name"]): float(r["bval_raw"]) for _, r in b_name_raw.iterrows()
        }

    # Fallback pos-based structures (per chrom): pos -> z, sorted pos
    pos_sorted_by_chr: Dict[str, np.ndarray] = {}
    z_by_chr_pos: Dict[str, Dict[int, float]] = {}
    raw_by_chr_pos: Dict[str, Dict[int, float]] = {}

    if args.fallback_pos_match:
        for chrom, sub in b.groupby("chrom", sort=False):
            tmpz = sub.groupby("pos", as_index=False)[args.z_col].max()
            zmap = dict(zip(tmpz["pos"].astype(int), tmpz[args.z_col].astype(float)))
            pos_sorted_by_chr[chrom] = np.array(sorted(zmap.keys()), dtype=int)
            z_by_chr_pos[chrom] = zmap

            if raw_available:
                tmpr = sub.groupby("pos", as_index=False)[args.raw_col].max()
                rawmap = dict(zip(tmpr["pos"].astype(int), tmpr[args.raw_col].astype(float)))
                raw_by_chr_pos[chrom] = rawmap

    # Raw scoring scale (global p99 unless overridden)
    raw_p99 = float(args.raw_p99)
    if args.score_from == "raw":
        if not raw_available:
            raise RuntimeError(f"--score-from raw requested but '{args.raw_col}' not found in boundary support TSV.")
        if raw_p99 <= 0:
            raw_p99 = float(np.nanpercentile(b[args.raw_col].to_numpy(dtype=float), 99))
            if not (raw_p99 > 0):
                raw_p99 = 1.0  # avoid div0

    out_rows = []

    # Process per chromosome in segment order so "previous segment" is defined
    for chrom, subseg in seg.groupby("chrom", sort=False):
        subseg = subseg.sort_values(["start", "end"]).reset_index(drop=True)

        for i in range(len(subseg)):
            r = subseg.loc[i].to_dict()
            seg_name = str(r["name"])

            # Right boundary after this segment: by name = seg_i (per your boundary_support TSV)
            right_z = bval_by_chr_name.get((chrom, seg_name), None)
            right_raw = bval_raw_by_chr_name.get((chrom, seg_name), None) if raw_available else None
            right_pos = None  # position not guaranteed in name-table; we can add later if you want

            # Left boundary before this segment: by name = previous segment
            if i == 0:
                left_z = None
                left_raw = None
                left_pos = None
            else:
                prev_name = str(subseg.loc[i - 1, "name"])
                left_z = bval_by_chr_name.get((chrom, prev_name), None)
                left_raw = bval_raw_by_chr_name.get((chrom, prev_name), None) if raw_available else None
                left_pos = None

            # Optional fallback: if name-based missing, approximate boundary positions between segments
            # boundary between seg_{i-1} and seg_i is midpoint of prev.end and this.start
            if args.fallback_pos_match:
                # fill missing left
                if i > 0 and left_z is None:
                    target = int((int(subseg.loc[i - 1, "end"]) + int(subseg.loc[i, "start"])) // 2)
                    pos_arr = pos_sorted_by_chr.get(chrom, np.array([], dtype=int))
                    zmap = z_by_chr_pos.get(chrom, {})
                    lp, lz = nearest_pos_value(pos_arr, zmap, target, args.tol_bp)
                    left_z = lz
                    left_pos = lp
                    if raw_available and chrom in raw_by_chr_pos and lp is not None:
                        left_raw = float(raw_by_chr_pos[chrom].get(lp, np.nan))

                # fill missing right
                if i < len(subseg) - 1 and right_z is None:
                    target = int((int(subseg.loc[i, "end"]) + int(subseg.loc[i + 1, "start"])) // 2)
                    pos_arr = pos_sorted_by_chr.get(chrom, np.array([], dtype=int))
                    zmap = z_by_chr_pos.get(chrom, {})
                    rp, rz = nearest_pos_value(pos_arr, zmap, target, args.tol_bp)
                    right_z = rz
                    right_pos = rp
                    if raw_available and chrom in raw_by_chr_pos and rp is not None:
                        right_raw = float(raw_by_chr_pos[chrom].get(rp, np.nan))

            # Combine into segment confidence (z-space)
            vals_z = [z for z in [left_z, right_z] if z is not None and not (isinstance(z, float) and math.isnan(z))]
            if len(vals_z) == 2:
                if args.method == "min":
                    seg_conf_z = float(min(vals_z))
                elif args.method == "max":
                    seg_conf_z = float(max(vals_z))
                else:
                    seg_conf_z = float(np.mean(vals_z))
                missing = ""
            elif len(vals_z) == 1:
                if args.missing_policy == "one_side":
                    seg_conf_z = float(vals_z[0])
                    missing = "right" if right_z is None else "left"
                else:
                    seg_conf_z = float("nan")
                    missing = "one_missing"
            else:
                seg_conf_z = float("nan")
                missing = "both_missing"

            # Also provide raw-space combined support rollup if available
            seg_conf_raw = float("nan")
            if raw_available:
                vals_raw = [x for x in [left_raw, right_raw] if x is not None and not (isinstance(x, float) and math.isnan(x))]
                if len(vals_raw) == 2:
                    seg_conf_raw = float(min(vals_raw)) if args.method == "min" else (
                        float(max(vals_raw)) if args.method == "max" else float(np.mean(vals_raw))
                    )
                elif len(vals_raw) == 1 and args.missing_policy == "one_side":
                    seg_conf_raw = float(vals_raw[0])

            # Score for BED
            if args.score_from == "z":
                bed_score = robust_score_0_1000_from_z(seg_conf_z, z_max=args.z_max) if not math.isnan(seg_conf_z) else 0
            else:
                bed_score = robust_score_0_1000_from_raw(seg_conf_raw, x_p99=raw_p99) if not (isinstance(seg_conf_raw, float) and math.isnan(seg_conf_raw)) else 0

            r.update(
                dict(
                    left_boundary_name=str(subseg.loc[i - 1, "name"]) if i > 0 else "",
                    right_boundary_name=seg_name,
                    left_boundary_pos=left_pos if left_pos is not None else "",
                    right_boundary_pos=right_pos if right_pos is not None else "",
                    left_boundary_z=left_z if left_z is not None else np.nan,
                    right_boundary_z=right_z if right_z is not None else np.nan,
                    seg_conf_z=seg_conf_z,
                    seg_conf_raw=seg_conf_raw,
                    seg_conf_method=args.method,
                    seg_conf_missing=missing,
                    seg_conf_score_0_1000=int(bed_score),
                    score_from=args.score_from,
                )
            )
            out_rows.append(r)

    out_df = pd.DataFrame(out_rows)

    # Write TSV
    if args.out_tsv.endswith(".gz"):
        with gzip.open(args.out_tsv, "wt") as f:
            out_df.to_csv(f, sep="\t", index=False)
    else:
        out_df.to_csv(args.out_tsv, sep="\t", index=False)

    # Optional BED for IGV
    if args.out_bed:
        bed = pd.DataFrame(
            {
                "chrom": out_df["chrom"],
                "start": out_df["start"].astype(int),
                "end": out_df["end"].astype(int),
                "name": out_df["name"].astype(str),
                "score": out_df["seg_conf_score_0_1000"].astype(int),
                "strand": out_df["strand"] if "strand" in out_df.columns else ".",
            }
        )
        if args.out_bed.endswith(".gz"):
            with gzip.open(args.out_bed, "wt") as f:
                bed.to_csv(f, sep="\t", header=False, index=False)
        else:
            bed.to_csv(args.out_bed, sep="\t", header=False, index=False)


if __name__ == "__main__":
    main()
