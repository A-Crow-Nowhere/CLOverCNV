#!/usr/bin/env python3
# ============================================================
# cnv_refuse.py
#
# "Force-fuser" + "Re-fuser" for CLOverCNV segments.calls.tsv(.gz)
#
# - force-fuse: fuse user-designated rows/groups regardless of tiers
# - re-fuse:    second-pass fusion using tunable eligibility knobs
#
# Key feature: direction-policy controls gain/loss mixing behavior
#   strict   : require same direction within a fused block
#   weighted : allow mixed; recompute direction from fused_ratio
#   ignore   : allow mixed; recompute direction from fused_ratio
#
# Outputs a *non-redundant* table by default (drops original rows
# that were fused, unless --keep-originals is set).
# ============================================================

import argparse
import gzip
import math
import os
import re
import tempfile
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd


# ----------------------------
# IO helpers
# ----------------------------

def _is_gzip_by_magic(path: str) -> bool:
    try:
        with open(path, "rb") as fh:
            return fh.read(2) == b"\x1f\x8b"
    except Exception:
        return False


def read_tsv_any(path: str, header="infer") -> pd.DataFrame:
    is_gz = _is_gzip_by_magic(path)
    if is_gz:
        return pd.read_csv(path, sep="\t", header=header, compression="gzip", low_memory=False)
    return pd.read_csv(path, sep="\t", header=header, low_memory=False)


def write_tsv_any(df: pd.DataFrame, path: str) -> None:
    if path.endswith(".gz"):
        with gzip.open(path, "wt") as f:
            df.to_csv(f, sep="\t", index=False, na_rep="NA")
    else:
        df.to_csv(path, sep="\t", index=False, na_rep="NA")


def _atomic_write_no_follow(write_fn, df: pd.DataFrame, out_path: str) -> None:
    out_dir = os.path.dirname(os.path.abspath(out_path)) or "."
    os.makedirs(out_dir, exist_ok=True)
    fd, tmp_path = tempfile.mkstemp(prefix=".tmp.", dir=out_dir)
    os.close(fd)
    try:
        write_fn(df, tmp_path)
        os.replace(tmp_path, out_path)
    finally:
        try:
            if os.path.exists(tmp_path):
                os.remove(tmp_path)
        except OSError:
            pass


# ----------------------------
# Core tier + direction logic (match finalize semantics)
# ----------------------------

def z_tier(z: float, weak_z: float, strong_z: float) -> str:
    if not math.isfinite(z):
        return "NA"
    if z >= strong_z:
        return "strong"
    if z >= weak_z:
        return "weak"
    return "low"


def ratio_tier(r: float, weak_ratio: float, strong_ratio: float) -> str:
    if r <= 0 or (not math.isfinite(r)):
        return "NA"
    if r >= strong_ratio or r <= (1.0 / strong_ratio):
        return "strong"
    if r >= weak_ratio or r <= (1.0 / weak_ratio):
        return "weak"
    return "low"


def call_direction(r: float) -> str:
    if not math.isfinite(r) or r <= 0:
        return "NA"
    if r > 1.0:
        return "gain"
    if r < 1.0:
        return "loss"
    return "neutral"


def _combine_two(a: float, b: float, policy: str) -> float:
    aa = float(a) if (a is not None and math.isfinite(float(a))) else float("nan")
    bb = float(b) if (b is not None and math.isfinite(float(b))) else float("nan")

    vals = [x for x in (aa, bb) if math.isfinite(x)]
    if not vals:
        return float("nan")

    pol = (policy or "min").lower()
    if pol == "min":
        return float(min(vals))
    if pol == "max":
        return float(max(vals))
    if pol == "mean":
        return float(np.mean(vals))
    raise ValueError(f"Unknown seg-conf policy: {policy} (expected min|mean|max)")


def weighted_mean(vals, weights) -> float:
    v = np.asarray(vals, dtype=float)
    w = np.asarray(weights, dtype=float)
    ok = np.isfinite(v) & np.isfinite(w) & (w > 0)
    if not ok.any():
        return float("nan")
    return float(np.sum(v[ok] * w[ok]) / np.sum(w[ok]))


def _sum_int(vals) -> pd._libs.missing.NAType | int:
    s = 0
    any_ok = False
    for x in vals:
        if x is None or x is pd.NA:
            continue
        if isinstance(x, float) and (not math.isfinite(x)):
            continue
        try:
            s += int(x)
            any_ok = True
        except Exception:
            continue
    return int(s) if any_ok else pd.NA


def _mean_float(vals) -> float:
    vv = []
    for x in vals:
        if x is None or x is pd.NA:
            continue
        if isinstance(x, float) and (not math.isfinite(x)):
            continue
        try:
            vv.append(float(x))
        except Exception:
            pass
    if not vv:
        return float("nan")
    return float(np.mean(vv))


def _as_float(x) -> float:
    try:
        xx = float(x)
        return xx
    except Exception:
        return float("nan")


# ----------------------------
# Name de-dup helper (fixes your double-first-name symptom)
# ----------------------------

def _dedup_names_preserve_order(names: List[str]) -> List[str]:
    seen = set()
    out = []
    for nm in names:
        nm = str(nm)
        if nm not in seen:
            out.append(nm)
            seen.add(nm)
    return out



def _collect_block_names(block: pd.DataFrame) -> list[str]:
    """
    Collect an ordered, de-duplicated list of segment names represented by rows in `block`.

    Priority:
      1) If fused_segment_names exists and is populated, expand it (comma-separated)
      2) Else fall back to the row's `name`

    Robust to NA/NaN, floats, and already-fused rows.
    """
    names: list[str] = []

    def _add_token(tok: str):
        tok = tok.strip()
        if tok and tok != "NA":
            names.append(tok)

    if "fused_segment_names" in block.columns:
        for v in block["fused_segment_names"].tolist():
            if v is None or v is pd.NA:
                continue
            # pandas NaN is float
            if isinstance(v, float):
                if not math.isfinite(v):
                    continue
                v = str(v)
            s = str(v).strip()
            if not s or s == "NA":
                continue
            for tok in s.split(","):
                _add_token(tok)

    # If we didn't get anything, fall back to `name`
    if not names and "name" in block.columns:
        for v in block["name"].tolist():
            if v is None or v is pd.NA:
                continue
            if isinstance(v, float):
                if not math.isfinite(v):
                    continue
                v = str(v)
            _add_token(str(v))

    # De-duplicate while preserving order
    seen = set()
    uniq: list[str] = []
    for nm in names:
        if nm not in seen:
            uniq.append(nm)
            seen.add(nm)
    return uniq



# ----------------------------
# Fusion engine (force or eligibility-driven)
# ----------------------------

def _weights_for_block(block: pd.DataFrame) -> np.ndarray:
    # weight priority: n_probes if present and >0, else seg_len_bp, else span length
    wp = pd.to_numeric(block.get("n_probes", np.nan), errors="coerce").to_numpy(dtype=float)
    wlen = pd.to_numeric(block.get("seg_len_bp", np.nan), errors="coerce").to_numpy(dtype=float)

    span = (
        pd.to_numeric(block["end"], errors="coerce") - pd.to_numeric(block["start"], errors="coerce")
    ).to_numpy(dtype=float)

    w = np.where(np.isfinite(wp) & (wp > 0), wp, np.nan)
    w = np.where(np.isfinite(w) & (w > 0), w, np.where(np.isfinite(wlen) & (wlen > 0), wlen, np.nan))
    w = np.where(np.isfinite(w) & (w > 0), w, np.where(np.isfinite(span) & (span > 0), span, 1.0))
    w = np.where(np.isfinite(w) & (w > 0), w, 1.0)
    return w


def _block_direction(block: pd.DataFrame, direction_policy: str, fused_ratio: float) -> str:
    pol = (direction_policy or "strict").lower()
    if pol not in ("strict", "weighted", "ignore"):
        raise ValueError(f"Unknown --direction-policy: {direction_policy}")

    # if strict: require identical direction labels across block (gain/loss)
    if pol == "strict":
        dirs = [str(x) for x in block.get("direction", "NA").astype(str).tolist()]
        dirs = [d for d in dirs if d not in ("NA", "neutral", "nan")]
        if len(dirs) == 0:
            return call_direction(fused_ratio)
        if len(set(dirs)) != 1:
            raise RuntimeError(f"direction-policy strict: mixed directions in block: {sorted(set(dirs))}")
        return dirs[0]

    # weighted/ignore: recompute from fused ratio
    return call_direction(fused_ratio)


def _recompute_call_fields(
    row: Dict,
    weak_ratio: float,
    strong_ratio: float,
    weak_z: float,
    strong_z: float,
    keep_mode: str,
    neutral_band: bool,
) -> Dict:
    ratio = _as_float(row.get("ratio", float("nan")))
    seg_conf_z = _as_float(row.get("seg_conf_z", float("nan")))
    direction = str(row.get("direction", "NA"))

    rtier = ratio_tier(ratio, weak_ratio, strong_ratio)
    ztier = z_tier(seg_conf_z, weak_z, strong_z)

    ratio_keep = (rtier in ("weak", "strong"))
    z_keep = (ztier in ("weak", "strong"))

    km = (keep_mode or "ratio_only").lower()
    if km == "ratio_only":
        keep_suggested = ratio_keep
    elif km == "z_only":
        keep_suggested = z_keep
    elif km == "both":
        keep_suggested = ratio_keep and z_keep
    else:
        raise ValueError(f"Unknown keep-mode: {keep_mode}")

    # neutral-band safety: if fused ratio is inside neutral band, force "no effect"
    if neutral_band and math.isfinite(ratio):
        if (1.0 / weak_ratio) < ratio < weak_ratio:
            rtier = "low"
            ratio_keep = False
            keep_suggested = False
            direction = call_direction(ratio)

    if not ratio_keep:
        call = "no_CNV_effect"
    else:
        call = f"CNV_{direction}_{rtier}Effect_{ztier}Support"

    row.update({
        "ratio_tier": rtier,
        "support_tier": ztier,
        "keep_suggested": bool(keep_suggested),
        "call": call,
        "direction": direction,
        "weak_ratio": float(weak_ratio),
        "strong_ratio": float(strong_ratio),
        "weak_z": float(weak_z),
        "strong_z": float(strong_z),
        "keep_mode": str(keep_mode),
    })
    return row


def fuse_one_block(
    block: pd.DataFrame,
    seg_conf_policy: str,
    direction_policy: str,
    weak_ratio: float,
    strong_ratio: float,
    weak_z: float,
    strong_z: float,
    keep_mode: str,
    neutral_band: bool,
) -> Tuple[Dict, List[int]]:
    """
    Fuse a *user-provided* block (already filtered/ordered).
    Returns (fused_row_dict, list_of_original_row_indices_in_input_df).
    """

    if block.empty:
        raise RuntimeError("Cannot fuse empty block")

    chroms = set(block["chrom"].astype(str).tolist())
    if len(chroms) != 1:
        raise RuntimeError(f"Force-fuse requires one chromosome per fused block; got: {sorted(chroms)}")
    chrom = list(chroms)[0]

    block = block.sort_values(["start", "end"]).copy()
    w = _weights_for_block(block)

    new_start = int(pd.to_numeric(block["start"], errors="coerce").min())
    new_end = int(pd.to_numeric(block["end"], errors="coerce").max())

    ratio_vals = pd.to_numeric(block.get("ratio", np.nan), errors="coerce").to_numpy(dtype=float)
    y_vals = pd.to_numeric(block.get("y", np.nan), errors="coerce").to_numpy(dtype=float)
    cn_vals = pd.to_numeric(block.get("cn", np.nan), errors="coerce").to_numpy(dtype=float)

    fused_ratio = weighted_mean(ratio_vals, w)
    fused_y = weighted_mean(y_vals, w)
    fused_cn = weighted_mean(cn_vals, w)

    # n_probes sum if possible
    nprobes_vals = pd.to_numeric(block.get("n_probes", np.nan), errors="coerce").to_numpy(dtype=float)
    fused_nprobes = float(np.nansum(nprobes_vals)) if np.isfinite(nprobes_vals).any() else float("nan")

    # boundary endpoints
    left_name = str(block.iloc[0].get("left_boundary_name", "NA"))
    right_name = str(block.iloc[-1].get("right_boundary_name", block.iloc[-1].get("name", "NA")))
    left_z = _as_float(block.iloc[0].get("left_boundary_z", float("nan")))
    right_z = _as_float(block.iloc[-1].get("right_boundary_z", float("nan")))

    fused_seg_conf_z = _combine_two(left_z, right_z, seg_conf_policy)

    # direction policy
    direction = _block_direction(block, direction_policy, fused_ratio)

    # recompute length + support counts
    fused_seg_len_bp = int(new_end - new_start)

    # seg_primary_reads: sum (counts are additive)
    fused_seg_primary_reads = _sum_int(block.get("seg_primary_reads", pd.Series([pd.NA]*len(block))).tolist())

    # For outer boundary support, use endpoint per-segment support columns if present.
    outer_left_support = _as_float(block.iloc[0].get("left_support_reads", float("nan")))
    outer_right_support = _as_float(block.iloc[-1].get("right_support_reads", float("nan")))

    # If split/disco are present as endpoint read counts, keep those too.
    outer_left_split = _as_float(block.iloc[0].get("left_split_reads", float("nan")))
    outer_left_disco = _as_float(block.iloc[0].get("left_disco_reads", float("nan")))
    outer_right_split = _as_float(block.iloc[-1].get("right_split_reads", float("nan")))
    outer_right_disco = _as_float(block.iloc[-1].get("right_disco_reads", float("nan")))

    # Internal boundaries: take right_support_reads of each segment except the last
    internal_support = []
    internal_pos = []
    if "right_support_reads" in block.columns:
        for k in range(len(block) - 1):
            bpos = int(pd.to_numeric(block.iloc[k].get("end", np.nan), errors="coerce"))
            sval = _as_float(block.iloc[k].get("right_support_reads", float("nan")))
            if math.isfinite(sval):
                internal_support.append(float(sval))
                internal_pos.append(bpos)

    mid = (new_start + new_end) / 2.0
    left_int = [s for s, p in zip(internal_support, internal_pos) if p <= mid]
    right_int = [s for s, p in zip(internal_support, internal_pos) if p > mid]

    internal_left_mean = float(np.mean(left_int)) if left_int else float("nan")
    internal_right_mean = float(np.mean(right_int)) if right_int else float("nan")

    # Build fused row from first row for compatibility, then overwrite everything important
    base = block.iloc[0].to_dict()

    # fused_segment_names (order-preserving, de-duplicated)
    fused_names = _collect_block_names(block)
    fused_names_str = ",".join([x for x in fused_names if x not in ("", "NA", "nan")])

    base.update({
        "chrom": str(chrom),
        "start": int(new_start),
        "end": int(new_end),

        # Keep deterministic name, but you can also make a new name if you prefer
        "name": str(block.iloc[0].get("name", "")),

        "y": float(fused_y),
        "ratio": float(fused_ratio),
        "cn": float(fused_cn),

        "n_probes": (int(fused_nprobes) if math.isfinite(fused_nprobes) else pd.NA),
        "direction": direction,

        "left_boundary_name": left_name,
        "right_boundary_name": right_name,
        "left_boundary_z": (float(left_z) if math.isfinite(left_z) else float("nan")),
        "right_boundary_z": (float(right_z) if math.isfinite(right_z) else float("nan")),
        "seg_conf_z": float(fused_seg_conf_z) if math.isfinite(fused_seg_conf_z) else float("nan"),

        # Support + length
        "seg_len_bp": int(fused_seg_len_bp),
        "seg_primary_reads": fused_seg_primary_reads,

        # Endpoint reads copied from endpoints (do NOT sum)
        "left_split_reads": (int(outer_left_split) if math.isfinite(outer_left_split) else pd.NA),
        "left_disco_reads": (int(outer_left_disco) if math.isfinite(outer_left_disco) else pd.NA),
        "left_support_reads": (int(outer_left_support) if math.isfinite(outer_left_support) else pd.NA),

        "right_split_reads": (int(outer_right_split) if math.isfinite(outer_right_split) else pd.NA),
        "right_disco_reads": (int(outer_right_disco) if math.isfinite(outer_right_disco) else pd.NA),
        "right_support_reads": (int(outer_right_support) if math.isfinite(outer_right_support) else pd.NA),

        # Fusion bookkeeping
        "fused": True,
        "fused_n_segments": int(len(block)),
        "fused_internal_boundaries": int(max(0, len(block) - 1)),
        "fused_segment_names": fused_names_str,

        # Preserve boundary evidence summaries
        "outer_left_support_reads": (float(outer_left_support) if math.isfinite(outer_left_support) else float("nan")),
        "outer_right_support_reads": (float(outer_right_support) if math.isfinite(outer_right_support) else float("nan")),
        "internal_left_boundary_support_mean": (float(internal_left_mean) if math.isfinite(internal_left_mean) else float("nan")),
        "internal_right_boundary_support_mean": (float(internal_right_mean) if math.isfinite(internal_right_mean) else float("nan")),
    })

    # Recompute tiers/keep/call based on *fused* ratio & fused seg_conf_z
    base = _recompute_call_fields(
        base,
        weak_ratio=weak_ratio,
        strong_ratio=strong_ratio,
        weak_z=weak_z,
        strong_z=strong_z,
        keep_mode=keep_mode,
        neutral_band=neutral_band,
    )

    # Mark call as fused (without breaking "no_CNV_effect")
    if isinstance(base.get("call", ""), str) and base["call"] and base["call"] != "no_CNV_effect":
        if not base["call"].endswith("_FUSED"):
            base["call"] = base["call"] + "_FUSED"

    # Return fused row, plus original indices to drop
    orig_idx = block.index.tolist()
    return base, orig_idx


# ----------------------------
# Eligibility logic for "re-fuser"
# ----------------------------

def _truthy_bool(x) -> bool:
    if x is True:
        return True
    if x is False:
        return False
    if x is None or x is pd.NA:
        return False
    s = str(x).strip().lower()
    if s in ("1", "true", "t", "yes", "y"):
        return True
    return False


def _eligible_mask(
    df: pd.DataFrame,
    require_keep_suggested: bool,
    min_ratio_tier: str,
    min_support_tier: str,
) -> pd.Series:
    # min tiers: low < weak < strong
    tier_rank = {"low": 0, "weak": 1, "strong": 2, "NA": -1, "nan": -1}

    rt = df.get("ratio_tier", "low").astype(str)
    st = df.get("support_tier", "low").astype(str)

    def _rank(series: pd.Series) -> pd.Series:
        return series.map(lambda z: tier_rank.get(str(z), -1)).astype(int)

    r_ok = _rank(rt) >= tier_rank.get(min_ratio_tier, 0)
    s_ok = _rank(st) >= tier_rank.get(min_support_tier, 0)

    keep_ok = pd.Series([True] * len(df), index=df.index)
    if require_keep_suggested and ("keep_suggested" in df.columns):
        keep_ok = df["keep_suggested"].map(_truthy_bool)

    return (r_ok & s_ok & keep_ok)


def refuse_pass(
    df: pd.DataFrame,
    max_gap: int,
    direction_policy: str,
    seg_conf_policy: str,
    weak_ratio: float,
    strong_ratio: float,
    weak_z: float,
    strong_z: float,
    keep_mode: str,
    neutral_band: bool,
    require_keep_suggested: bool,
    min_ratio_tier: str,
    min_support_tier: str,
) -> Tuple[pd.DataFrame, List[int]]:
    """
    Performs a second-pass fusion across eligible adjacent rows on each chrom.
    Returns (fused_rows_df, list_of_original_indices_fused_away).
    """
    df = df.copy()

    # normalize types
    df["chrom"] = df["chrom"].astype(str)
    df["start"] = pd.to_numeric(df["start"], errors="coerce").astype(int)
    df["end"] = pd.to_numeric(df["end"], errors="coerce").astype(int)

    eligible = _eligible_mask(df, require_keep_suggested, min_ratio_tier, min_support_tier)

    fused_rows: List[Dict] = []
    fused_away: List[int] = []

    for chrom, sub in df.groupby("chrom", sort=False):
        sub = sub.sort_values(["start", "end"]).copy()
        sub["_orig_idx"] = sub.index

        # safe alignment of eligibility (no KeyError on weird indices)
        sub["_eligible"] = (
            eligible.reindex(sub["_orig_idx"])
                    .fillna(False)
                    .to_numpy(dtype=bool)
        )

        i = 0
        while i < len(sub):
            if not bool(sub.iloc[i]["_eligible"]):
                i += 1
                continue

            # Start a candidate block (contiguous eligible with gap<=max_gap)
            block_rows = [i]
            j = i + 1

            while j < len(sub):
                if not bool(sub.iloc[j]["_eligible"]):
                    break

                gap = int(sub.iloc[j]["start"]) - int(sub.iloc[j - 1]["end"])
                if gap > int(max_gap):
                    break

                block_rows.append(j)
                j += 1

            block = sub.iloc[block_rows].copy()

            # If strict direction-policy is requested, this is enforced inside fuse_one_block
            if len(block) > 1:
                fused_row, orig_idx = fuse_one_block(
                    block=block.drop(columns=["_orig_idx", "_eligible"], errors="ignore"),
                    seg_conf_policy=seg_conf_policy,
                    direction_policy=direction_policy,
                    weak_ratio=weak_ratio,
                    strong_ratio=strong_ratio,
                    weak_z=weak_z,
                    strong_z=strong_z,
                    keep_mode=keep_mode,
                    neutral_band=neutral_band,
                )
                fused_rows.append(fused_row)
                fused_away.extend(orig_idx)

            i = j

    fused_df = pd.DataFrame(fused_rows) if fused_rows else pd.DataFrame(columns=df.columns)
    return fused_df, fused_away


# ----------------------------
# Force-fuse grouping
# ----------------------------

def _parse_force_group_specs(specs: List[str]) -> List[List[str]]:
    """
    Each --force-group accepts a comma-separated list of names.
    Example:
      --force-group "A,B,C" --force-group "D,E"
    """
    out = []
    for s in specs:
        s = (s or "").strip()
        if not s:
            continue
        names = [x.strip() for x in s.split(",") if x.strip()]
        if names:
            out.append(names)
    return out


def _ensure_same_chrom(df_block: pd.DataFrame, allow_multi_chrom: bool) -> None:
    chroms = set(df_block["chrom"].astype(str).tolist())
    if len(chroms) == 1:
        return
    if allow_multi_chrom:
        return
    raise RuntimeError(f"Force-fuse block spans multiple chromosomes: {sorted(chroms)}")


def force_fuse(
    df: pd.DataFrame,
    groups: List[List[str]],
    by_col: str,
    allow_multi_chrom: bool,
    seg_conf_policy: str,
    direction_policy: str,
    weak_ratio: float,
    strong_ratio: float,
    weak_z: float,
    strong_z: float,
    keep_mode: str,
    neutral_band: bool,
) -> Tuple[pd.DataFrame, List[int]]:
    """
    Force-fuse either:
      - explicit --force-group name lists (matches df['name'])
      - or group by a column (e.g. manual_group)
    Returns (fused_rows_df, list_of_original_indices_fused_away).
    """
    df = df.copy()
    df["chrom"] = df["chrom"].astype(str)
    df["start"] = pd.to_numeric(df["start"], errors="coerce").astype(int)
    df["end"] = pd.to_numeric(df["end"], errors="coerce").astype(int)

    fused_rows: List[Dict] = []
    fused_away: List[int] = []

    if groups:
        for names in groups:
            block = df[df["name"].astype(str).isin([str(x) for x in names])].copy()
            if block.empty:
                raise RuntimeError(f"--force-group matched no rows: {names}")
            _ensure_same_chrom(block, allow_multi_chrom)
            fused_row, orig_idx = fuse_one_block(
                block=block,
                seg_conf_policy=seg_conf_policy,
                direction_policy=direction_policy,
                weak_ratio=weak_ratio,
                strong_ratio=strong_ratio,
                weak_z=weak_z,
                strong_z=strong_z,
                keep_mode=keep_mode,
                neutral_band=neutral_band,
            )
            fused_rows.append(fused_row)
            fused_away.extend(orig_idx)

        fused_df = pd.DataFrame(fused_rows)
        return fused_df, fused_away

    # Otherwise, group by column
    if not by_col:
        raise RuntimeError("force-fuse requires either --force-group or --force-by-col")

    if by_col not in df.columns:
        raise RuntimeError(f"--force-by-col {by_col} not found in input columns")

    for gval, block in df.groupby(by_col, sort=False):
        # skip NA/blank groups
        if gval is None or gval is pd.NA:
            continue
        if isinstance(gval, float) and (not math.isfinite(gval)):
            continue
        if str(gval).strip() in ("", "NA", "nan"):
            continue

        block = block.copy()
        _ensure_same_chrom(block, allow_multi_chrom)

        # only fuse if block has >1 row
        if len(block) <= 1:
            continue

        fused_row, orig_idx = fuse_one_block(
            block=block,
            seg_conf_policy=seg_conf_policy,
            direction_policy=direction_policy,
            weak_ratio=weak_ratio,
            strong_ratio=strong_ratio,
            weak_z=weak_z,
            strong_z=strong_z,
            keep_mode=keep_mode,
            neutral_band=neutral_band,
        )
        # record group in output (helpful)
        fused_row[by_col] = gval

        fused_rows.append(fused_row)
        fused_away.extend(orig_idx)

    fused_df = pd.DataFrame(fused_rows) if fused_rows else pd.DataFrame(columns=df.columns)
    return fused_df, fused_away


# ----------------------------
# Output shaping
# ----------------------------

def _drop_redundant_fused_names(df: pd.DataFrame) -> pd.DataFrame:
    """
    Fix pathological fused_segment_names like:
      A,A,B,C
    by de-duplicating in-place.
    """
    if "fused_segment_names" not in df.columns:
        return df
    out = df.copy()
    fixed = []
    for s in out["fused_segment_names"].astype(str).tolist():
        if not s or s in ("NA", "nan"):
            fixed.append(s)
            continue
        names = [x for x in s.split(",") if x]
        names = _dedup_names_preserve_order(names)
        fixed.append(",".join(names))
    out["fused_segment_names"] = fixed
    return out


def _sort_stable(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["chrom"] = out["chrom"].astype(str)
    out["start"] = pd.to_numeric(out["start"], errors="coerce").astype(int)
    out["end"] = pd.to_numeric(out["end"], errors="coerce").astype(int)

    if "fused" in out.columns:
        out["_fused_sort"] = out["fused"].map(_truthy_bool).astype(int)
        out = out.sort_values(["chrom", "start", "end", "_fused_sort"], kind="mergesort").drop(columns=["_fused_sort"])
    else:
        out = out.sort_values(["chrom", "start", "end"], kind="mergesort")

    return out


def _ensure_fused_bool(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    if "fused" not in out.columns:
        out["fused"] = False
    else:
        out["fused"] = out["fused"].map(_truthy_bool)
    return out


# ----------------------------
# Main
# ----------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Post-hoc fuser for CLOverCNV segments.calls.tsv(.gz): force-fuse and re-fuse with tunable knobs."
    )

    ap.add_argument("--in-tsv", required=True, help="Input segments.calls.tsv(.gz ok)")
    ap.add_argument("--out-tsv", required=True, help="Output TSV(.gz ok)")

    ap.add_argument("--mode", choices=["force", "refuse"], required=True,
                    help="force: fuse user-designated groups; refuse: second-pass eligibility-driven fusion")

    # direction mixing behavior
    ap.add_argument("--direction-policy", choices=["strict", "weighted", "ignore"], default="strict",
                    help="strict: require same direction; weighted/ignore: allow mixed and recompute from fused ratio")

    # seg_conf_z combination behavior for fused rows
    ap.add_argument("--seg-conf-policy", choices=["min", "mean", "max"], default="mean",
                    help="How to compute fused seg_conf_z from outer left/right boundary z")

    # tier thresholds + keep-mode to recompute call fields
    ap.add_argument("--weak-ratio", type=float, default=1.25)
    ap.add_argument("--strong-ratio", type=float, default=1.50)
    ap.add_argument("--weak-z", type=float, default=1.0)
    ap.add_argument("--strong-z", type=float, default=2.0)
    ap.add_argument("--keep-mode", choices=["ratio_only", "z_only", "both"], default="both")

    ap.add_argument("--neutral-band", action="store_true",
                    help="If fused_ratio is within (1/weak_ratio, weak_ratio), force no_CNV_effect")

    # output shaping
    ap.add_argument("--keep-originals", action="store_true",
                    help="Keep original rows that were fused (default drops them)")
    ap.add_argument("--dedup-fused-names", action="store_true",
                    help="De-duplicate fused_segment_names in output (recommended)")
    ap.add_argument("--emit-fused-only", action="store_true",
                    help="Output only fused rows (drops all non-fused rows)")

    # --- force mode options ---
    ap.add_argument("--force-group", action="append", default=[],
                    help="Comma-separated segment 'name' list to fuse into one row; repeatable")
    ap.add_argument("--force-by-col", default="",
                    help="Fuse each group defined by this column (e.g. manual_group)")
    ap.add_argument("--allow-multi-chrom", action="store_true",
                    help="Allow force-fuse groups to span multiple chromosomes (NOT recommended)")

    # --- refuse mode options ---
    ap.add_argument("--max-gap", type=int, default=0, help="Max bp gap between adjacent rows to be fused")
    ap.add_argument("--require-keep-suggested", action="store_true",
                    help="Only fuse rows where keep_suggested is TRUE")
    ap.add_argument("--min-ratio-tier", choices=["low", "weak", "strong"], default="weak",
                    help="Minimum ratio_tier required for eligibility")
    ap.add_argument("--min-support-tier", choices=["low", "weak", "strong"], default="low",
                    help="Minimum support_tier required for eligibility")

    args = ap.parse_args()

    # sanity
    if args.weak_ratio <= 1.0 or args.strong_ratio <= 1.0:
        raise SystemExit("--weak-ratio and --strong-ratio must be > 1.0")
    if args.strong_ratio < args.weak_ratio:
        raise SystemExit("--strong-ratio must be >= --weak-ratio")
    if args.strong_z < args.weak_z:
        raise SystemExit("--strong-z must be >= --weak-z")

    df = read_tsv_any(args.in_tsv, header="infer")

    required = {"chrom", "start", "end", "name"}
    missing = required - set(df.columns)
    if missing:
        raise SystemExit(f"Input is missing required columns: {sorted(missing)}")

    df = _ensure_fused_bool(df)

    fused_df = pd.DataFrame()
    fused_away: List[int] = []

    if args.mode == "force":
        groups = _parse_force_group_specs(args.force_group)
        by_col = (args.force_by_col or "").strip()

        fused_df, fused_away = force_fuse(
            df=df,
            groups=groups,
            by_col=by_col,
            allow_multi_chrom=bool(args.allow_multi_chrom),
            seg_conf_policy=args.seg_conf_policy,
            direction_policy=args.direction_policy,
            weak_ratio=args.weak_ratio,
            strong_ratio=args.strong_ratio,
            weak_z=args.weak_z,
            strong_z=args.strong_z,
            keep_mode=args.keep_mode,
            neutral_band=bool(args.neutral_band),
        )

    else:
        fused_df, fused_away = refuse_pass(
            df=df,
            max_gap=int(args.max_gap),
            direction_policy=args.direction_policy,
            seg_conf_policy=args.seg_conf_policy,
            weak_ratio=args.weak_ratio,
            strong_ratio=args.strong_ratio,
            weak_z=args.weak_z,
            strong_z=args.strong_z,
            keep_mode=args.keep_mode,
            neutral_band=bool(args.neutral_band),
            require_keep_suggested=bool(args.require_keep_suggested),
            min_ratio_tier=str(args.min_ratio_tier),
            min_support_tier=str(args.min_support_tier),
        )

    # combine outputs
    out = df.copy()

    # Drop originals that were fused away unless requested
    if fused_away and (not args.keep_originals):
        out = out.drop(index=list(set(fused_away)), errors="ignore")

    # Append fused rows
    if not fused_df.empty:
        fused_df = _ensure_fused_bool(fused_df)
        out = pd.concat([out, fused_df], ignore_index=True, sort=False)

    # Optional output filtering
    if args.emit_fused_only:
        out = out[out["fused"].map(_truthy_bool)].copy()

    # Fix fused_segment_names redundancy if requested
    if args.dedup_fused_names:
        out = _drop_redundant_fused_names(out)

    out = _sort_stable(out)

    _atomic_write_no_follow(write_tsv_any, out, args.out_tsv)


if __name__ == "__main__":
    main()
