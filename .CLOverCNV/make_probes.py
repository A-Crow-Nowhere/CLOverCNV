#!/usr/bin/env python3
import argparse
import gzip
import math
import os
import sys
from typing import Dict, List, Tuple, Optional

import numpy as np
import pandas as pd
import pyBigWig


def open_text(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def parse_gff_attributes(attr: str) -> Dict[str, str]:
    # GFF3 attributes: key=value;key2=value2
    d = {}
    for part in attr.strip().split(";"):
        if not part:
            continue
        if "=" in part:
            k, v = part.split("=", 1)
            d[k] = v
        elif " " in part:
            # tolerate GTF-ish
            k, v = part.split(" ", 1)
            d[k] = v.strip('"')
    return d


def load_cds_intervals(gff_path: str, allowed_chroms: Optional[set] = None) -> pd.DataFrame:
    rows = []
    with open_text(gff_path) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            chrom, source, ftype, start, end, score, strand, phase, attrs = parts
            if ftype != "CDS":
                continue
            if allowed_chroms is not None and chrom not in allowed_chroms:
                continue
            try:
                s = int(start) - 1  # GFF is 1-based inclusive; convert to 0-based half-open
                e = int(end)
            except ValueError:
                continue
            if e <= s:
                continue
            ad = parse_gff_attributes(attrs)

            # Best-effort gene identifier:
            # - Parent often points to transcript; ID sometimes is the CDS feature ID
            # - gene / Name may exist
            gene_id = ad.get("gene", None) or ad.get("Name", None) or ad.get("Parent", None) or ad.get("ID", None) or "."
            rows.append((chrom, s, e, strand, gene_id))
    if not rows:
        raise RuntimeError(f"No CDS features found in GFF: {gff_path}")
    df = pd.DataFrame(rows, columns=["chrom", "start", "end", "strand", "gene_id"])
    df = df.sort_values(["chrom", "start", "end"]).reset_index(drop=True)
    return df


def load_bed_intervals(bed_path: str, allowed_chroms: Optional[set] = None) -> Dict[str, List[Tuple[int, int]]]:
    ex: Dict[str, List[Tuple[int, int]]] = {}
    with open_text(bed_path) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            chrom = parts[0]
            if allowed_chroms is not None and chrom not in allowed_chroms:
                continue
            try:
                s = int(parts[1]); e = int(parts[2])
            except ValueError:
                continue
            if e <= s:
                continue
            ex.setdefault(chrom, []).append((s, e))
    # merge overlaps per chrom
    for chrom, ivs in ex.items():
        ivs.sort()
        merged = []
        for s, e in ivs:
            if not merged or s > merged[-1][1]:
                merged.append([s, e])
            else:
                merged[-1][1] = max(merged[-1][1], e)
        ex[chrom] = [(a, b) for a, b in merged]
    return ex


def subtract_intervals(intervals: List[Tuple[int, int]], excludes: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """Subtract excludes from intervals (both sorted, half-open)."""
    if not excludes:
        return intervals
    out = []
    j = 0
    for s, e in intervals:
        cur_s, cur_e = s, e
        while j < len(excludes) and excludes[j][1] <= cur_s:
            j += 1
        k = j
        while k < len(excludes) and excludes[k][0] < cur_e:
            xs, xe = excludes[k]
            if xe <= cur_s:
                k += 1
                continue
            if xs > cur_s:
                out.append((cur_s, min(xs, cur_e)))
            cur_s = max(cur_s, xe)
            if cur_s >= cur_e:
                break
            k += 1
        if cur_s < cur_e:
            out.append((cur_s, cur_e))
    return out


def build_cds_tiles(cds_df: pd.DataFrame,
                    tile_bp: int,
                    min_tile_bp: int,
                    exclude_map: Optional[Dict[str, List[Tuple[int, int]]]] = None) -> pd.DataFrame:
    rows = []
    # group CDS by chrom; process each CDS interval independently (does not tile across introns)
    for (chrom, strand, gene_id), sub in cds_df.groupby(["chrom", "strand", "gene_id"], sort=False):
        # We tile each CDS interval (row) independently, preserving gene_id for later aggregation/reporting.
        for _, r in sub.iterrows():
            ivs = [(int(r["start"]), int(r["end"]))]
            if exclude_map is not None and chrom in exclude_map:
                ivs = subtract_intervals(ivs, exclude_map[chrom])
            for s, e in ivs:
                L = e - s
                if L < min_tile_bp:
                    continue
                pos = s
                while pos < e:
                    t_end = min(pos + tile_bp, e)
                    t_len = t_end - pos
                    if t_len >= min_tile_bp:
                        rows.append((chrom, pos, t_end, strand, gene_id))
                    pos = t_end
    if not rows:
        raise RuntimeError("No CDS tiles produced (check GFF, excludes, and min_tile_bp).")
    out = pd.DataFrame(rows, columns=["chrom", "start", "end", "strand", "gene_id"])
    out = out.sort_values(["chrom", "start", "end"]).reset_index(drop=True)
    out["len_bp"] = out["end"] - out["start"]
    out["probe_id"] = [f"cds_{i:07d}" for i in range(1, len(out) + 1)]
    out["feature"] = "CDS_tile"
    return out


def bw_weighted_mean(bw: pyBigWig.pyBigWig, chrom: str, start: int, end: int) -> Tuple[float, int]:
    """Overlap-weighted mean of bigWig intervals over [start,end). Returns (mean, n_spans)."""
    ivs = bw.intervals(chrom, start, end)
    if ivs is None or len(ivs) == 0:
        return (float("nan"), 0)
    total = 0.0
    denom = 0
    n = 0
    for s, e, v in ivs:
        if v is None:
            continue
        # clip to requested window
        os_ = max(start, int(s))
        oe_ = min(end, int(e))
        if oe_ <= os_:
            continue
        w = oe_ - os_
        total += float(v) * w
        denom += w
        n += 1
    if denom == 0:
        return (float("nan"), n)
    return (total / denom, n)


def write_bed_gz(df: pd.DataFrame, out_bed_gz: str):
    # IGV-friendly BED: chrom start end name score strand
    with gzip.open(out_bed_gz, "wt") as out:
        for _, r in df.iterrows():
            out.write(
                f"{r['chrom']}\t{int(r['start'])}\t{int(r['end'])}\t{r['probe_id']}\t0\t{r['strand']}\n"
            )


def main():
    ap = argparse.ArgumentParser(
        description="Step1: Build CDS tiles and compute overlap-weighted means from corrected primary bigWig."
    )
    ap.add_argument("--genome", required=True, help="GENOME_KEY or path to fasta[.gz]. (Not used in step1; kept for pipeline consistency.)")
    ap.add_argument("--gff", required=True, help="Path to GFF/GFF3 (.gz ok).")
    ap.add_argument("--primary-bw", required=True, help="Corrected primary coverage bigWig.")
    ap.add_argument("--exclude-bed", default="", help="Optional BED of regions to exclude (.gz ok).")
    ap.add_argument("--tile-bp", type=int, default=100, help="Tile size within CDS space.")
    ap.add_argument("--min-tile-bp", type=int, default=50, help="Minimum tile length to keep.")
    ap.add_argument("--out-prefix", required=True, help="Output prefix directory (files written under this).")
    ap.add_argument("--eps", type=float, default=1e-6, help="Pseudocount for log2 ratio.")
    args = ap.parse_args()

    out_dir = args.out_prefix
    os.makedirs(out_dir, exist_ok=True)

    # Open BW and get chrom whitelist (prevents GFF contigs not present in track)
    bw = pyBigWig.open(args.primary_bw)
    bw_chroms = set(bw.chroms().keys())

    print(f"[make_probes] primary_bw chroms: {len(bw_chroms)}", file=sys.stderr)
    print(f"[make_probes] loading CDS from: {args.gff}", file=sys.stderr)
    cds_df = load_cds_intervals(args.gff, allowed_chroms=bw_chroms)

    exclude_map = None
    if args.exclude_bed:
        print(f"[make_probes] loading exclude bed: {args.exclude_bed}", file=sys.stderr)
        exclude_map = load_bed_intervals(args.exclude_bed, allowed_chroms=bw_chroms)

    print(f"[make_probes] building CDS tiles: tile_bp={args.tile_bp} min_tile_bp={args.min_tile_bp}", file=sys.stderr)
    probes = build_cds_tiles(cds_df, args.tile_bp, args.min_tile_bp, exclude_map=exclude_map)
    print(f"[make_probes] probes: {len(probes)}", file=sys.stderr)

    # Compute primary mean per probe
    means = []
    nspans = []
    for _, r in probes.iterrows():
        m, n = bw_weighted_mean(bw, r["chrom"], int(r["start"]), int(r["end"]))
        means.append(m)
        nspans.append(n)
    probes["primary_mean"] = means
    probes["n_bw_spans"] = nspans

    # Baseline per chromosome (median of primary_mean)
    probes["chrom_median"] = np.nan
    for chrom, sub in probes.groupby("chrom", sort=False):
        med = np.nanmedian(sub["primary_mean"].values.astype(float))
        probes.loc[sub.index, "chrom_median"] = med

    # log2 ratio to chrom median
    eps = float(args.eps)
    probes["y"] = np.log2((probes["primary_mean"].astype(float) + eps) / (probes["chrom_median"].astype(float) + eps))

    # excluded_bp is tracked implicitly by subtraction; for v1 we report 0 (tile coords already reflect exclusion)
    probes["excluded_bp"] = 0
    probes["weight"] = 1.0  # v1 placeholder; v2 can estimate within-probe variability
    # Reorder columns
    cols = [
        "probe_id", "chrom", "start", "end", "strand", "gene_id", "feature", "len_bp",
        "primary_mean", "chrom_median", "y",
        "weight", "excluded_bp", "n_bw_spans"
    ]
    probes = probes[cols]

    # Write outputs
    bed_out = os.path.join(out_dir, "probes.cds_tiles.bed.gz")
    tsv_out = os.path.join(out_dir, "probe_table.tsv.gz")

    print(f"[make_probes] writing: {bed_out}", file=sys.stderr)
    write_bed_gz(probes, bed_out)

    print(f"[make_probes] writing: {tsv_out}", file=sys.stderr)
    probes.to_csv(tsv_out, sep="\t", index=False, compression="gzip")

    bw.close()
    print("[make_probes] done", file=sys.stderr)


if __name__ == "__main__":
    main()
