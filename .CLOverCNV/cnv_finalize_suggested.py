#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import math
import sys
from dataclasses import dataclass
from typing import List, Tuple

import numpy as np
import pandas as pd


# ============================================================
# I/O
# ============================================================

def is_gzip_path(path: str) -> bool:
    try:
        with open(path, "rb") as fh:
            return fh.read(2) == b"\x1f\x8b"
    except OSError:
        return False


def open_maybe_gzip_read(path: str):
    if is_gzip_path(path):
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return open(path, "rt", encoding="utf-8", newline="")


def open_maybe_gzip_write(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "wt", encoding="utf-8", newline="")
    return open(path, "wt", encoding="utf-8", newline="")


def read_tsv(path: str) -> pd.DataFrame:
    with open_maybe_gzip_read(path) as fh:
        return pd.read_csv(
            fh,
            sep="\t",
            dtype=str,
            keep_default_na=True,
            na_values=["", "NA", "NaN", "nan"],
        )


def write_tsv(df: pd.DataFrame, path: str) -> None:
    with open_maybe_gzip_write(path) as fh:
        df.to_csv(fh, sep="\t", index=False, na_rep="NA")


# ============================================================
# Utils
# ============================================================

NUMERIC_COLUMNS_HINT = [
    "start", "end", "y", "ratio", "n_probes", "cn_base", "cn",
    "left_boundary_z", "right_boundary_z", "seg_conf_z",
    "weak_ratio", "strong_ratio", "weak_z", "strong_z",
    "seg_len_bp", "seg_primary_reads",
    "left_split_reads", "left_disco_reads", "left_support_reads",
    "right_split_reads", "right_disco_reads", "right_support_reads",
    "fused_n_segments", "fused_internal_boundaries",
    "internal_left_boundary_support_mean", "internal_right_boundary_support_mean",
    "fused_left_boundaries_n", "fused_right_boundaries_n",
    "outer_left_split_reads", "outer_left_disco_reads", "outer_left_support_reads",
    "outer_right_split_reads", "outer_right_disco_reads", "outer_right_support_reads",
    "ratio_z",
]


def to_float(x) -> float:
    if x is None:
        return float("nan")
    if isinstance(x, (float, int, np.floating, np.integer)):
        return float(x)
    s = str(x).strip()
    if s == "" or s.upper() == "NA":
        return float("nan")
    try:
        return float(s)
    except ValueError:
        return float("nan")


def to_bool(x) -> bool:
    if isinstance(x, bool):
        return x
    if x is None:
        return False
    if isinstance(x, float) and math.isnan(x):
        return False
    s = str(x).strip().lower()
    return s in {"true", "t", "1", "yes", "y"}


def bool_str(x: bool) -> str:
    return "TRUE" if x else "FALSE"


def coerce_numeric_columns(df: pd.DataFrame) -> pd.DataFrame:
    for col in df.columns:
        if col in NUMERIC_COLUMNS_HINT:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def get_numeric(row: pd.Series, col: str) -> float:
    if col not in row.index:
        return float("nan")
    return to_float(row[col])


def set_na(df: pd.DataFrame, idx: int, cols: List[str]) -> None:
    for col in cols:
        if col in df.columns:
            df.at[idx, col] = np.nan


def minmax_ratio(a: float, b: float) -> float:
    if math.isnan(a) or math.isnan(b):
        return float("nan")
    mx = max(a, b)
    mn = min(a, b)
    if mx <= 0:
        return 0.0
    return mn / mx


def harmonic_mean(a: float, b: float) -> float:
    if math.isnan(a) or math.isnan(b) or a <= 0 or b <= 0:
        return 0.0
    return 2.0 * a * b / (a + b)


def dedup_keep_order(items: List[str]) -> List[str]:
    out: List[str] = []
    seen = set()
    for x in items:
        if not x or x in seen:
            continue
        seen.add(x)
        out.append(x)
    return out


# ============================================================
# Support extraction
# ============================================================

def support_count_for_side(row: pd.Series, side: str) -> float:
    fused = to_bool(row.get("fused", False))

    if fused:
        v = get_numeric(row, f"outer_{side}_support_reads")
        if not math.isnan(v):
            return v

        split_v = get_numeric(row, f"outer_{side}_split_reads")
        disco_v = get_numeric(row, f"outer_{side}_disco_reads")
        if not math.isnan(split_v) or not math.isnan(disco_v):
            return (0.0 if math.isnan(split_v) else split_v) + (0.0 if math.isnan(disco_v) else disco_v)

    v = get_numeric(row, f"{side}_support_reads")
    if not math.isnan(v):
        return v

    split_v = get_numeric(row, f"{side}_split_reads")
    disco_v = get_numeric(row, f"{side}_disco_reads")
    if not math.isnan(split_v) or not math.isnan(disco_v):
        return (0.0 if math.isnan(split_v) else split_v) + (0.0 if math.isnan(disco_v) else disco_v)

    return float("nan")


# ============================================================
# Dataclasses
# ============================================================

@dataclass
class SideInfo:
    z_raw: float
    z_eff: float
    z_missing: bool
    support_count: float


@dataclass
class SegmentEval:
    left: SideInfo
    right: SideInfo

    keep_ratio: bool
    keep_z: bool
    keep_balance: bool
    keep_neighbor: bool

    rescue_applies: bool
    rescue_side: str
    balance_source: str
    balance_ratio: float
    balance_score: float

    imbalance_override: bool
    reasons_override: List[str]

    shared_boundary_suspect: bool
    neighbor_competitor: str

    reasons_fail: List[str]
    reasons_rescue: List[str]
    reasons_note: List[str]
    reason_summary: str


# ============================================================
# Core logic
# ============================================================

def build_side_info(row: pd.Series, side: str) -> SideInfo:
    z_raw = get_numeric(row, f"{side}_boundary_z")
    z_missing = math.isnan(z_raw)
    z_eff = 0.0 if z_missing else max(z_raw, 0.0)
    support_count = support_count_for_side(row, side)
    return SideInfo(
        z_raw=z_raw,
        z_eff=z_eff,
        z_missing=z_missing,
        support_count=support_count,
    )


def choose_balance_pair(left: SideInfo, right: SideInfo) -> Tuple[float, float, str]:
    if (not left.z_missing) and (not right.z_missing):
        return left.z_eff, right.z_eff, "z"
    return left.support_count, right.support_count, "count"


def compute_keep_ratio(row: pd.Series) -> bool:
    ratio_tier = str(row.get("ratio_tier", "NA")).strip().lower()
    weak_ratio = to_float(row.get("weak_ratio"))
    ratio = to_float(row.get("ratio"))

    if ratio_tier in {"weak", "strong"}:
        return True

    if not math.isnan(ratio) and not math.isnan(weak_ratio):
        if ratio >= weak_ratio:
            return True
        if weak_ratio > 0 and ratio <= (1.0 / weak_ratio):
            return True

    return False


def compute_keep_z(left: SideInfo, right: SideInfo, weak_z: float, rescue_z: float) -> Tuple[bool, bool, str, List[str], List[str]]:
    fail_reasons: List[str] = []
    rescue_reasons: List[str] = []
    rescue_applies = False
    rescue_side = "none"

    left_below = (not left.z_missing) and (left.z_eff < weak_z)
    right_below = (not right.z_missing) and (right.z_eff < weak_z)

    if left_below and right.z_eff >= rescue_z:
        rescue_applies = True
        rescue_side = "right"
        rescue_reasons.append("left_boundary_below_weak_z_rescued_by_extreme_right_z")

    if right_below and left.z_eff >= rescue_z:
        rescue_applies = True
        rescue_side = "left"
        rescue_reasons.append("right_boundary_below_weak_z_rescued_by_extreme_left_z")

    if rescue_applies:
        return True, True, rescue_side, fail_reasons, rescue_reasons

    present_vals = []
    if not left.z_missing:
        present_vals.append(left.z_eff)
    if not right.z_missing:
        present_vals.append(right.z_eff)

    if not present_vals:
        fail_reasons.append("both_boundary_z_NA")
        return False, False, rescue_side, fail_reasons, rescue_reasons

    if any(v < weak_z for v in present_vals):
        if not any(v >= weak_z for v in present_vals):
            fail_reasons.append("no_boundary_meets_weak_z")
        else:
            fail_reasons.append("one_boundary_below_weak_z")
        return False, False, rescue_side, fail_reasons, rescue_reasons

    return True, False, rescue_side, fail_reasons, rescue_reasons


def compute_keep_balance(left: SideInfo, right: SideInfo, balance_threshold: float, rescue_applies: bool) -> Tuple[bool, float, float, str, List[str]]:
    fail_reasons: List[str] = []

    a, b, source = choose_balance_pair(left, right)
    ratio = minmax_ratio(a, b)
    score = harmonic_mean(0.0 if math.isnan(a) else a, 0.0 if math.isnan(b) else b)

    if rescue_applies:
        return True, ratio, score, source, fail_reasons

    if math.isnan(ratio):
        fail_reasons.append("balance_not_computable")
        return False, ratio, score, source, fail_reasons

    if ratio < balance_threshold:
        fail_reasons.append("boundary_balance_below_threshold")
        return False, ratio, score, source, fail_reasons

    return True, ratio, score, source, fail_reasons


def compute_imbalance_override(row: pd.Series, left: SideInfo, right: SideInfo, strong_z: float) -> Tuple[bool, List[str]]:
    reasons: List[str] = []

    ratio_tier = str(row.get("ratio_tier", "NA")).strip().lower()
    left_strong = (not left.z_missing) and (left.z_eff >= strong_z)
    right_strong = (not right.z_missing) and (right.z_eff >= strong_z)

    if ratio_tier == "strong":
        reasons.append("imbalance_overridden_by_strong_ratio")
        return True, reasons

    if left_strong and right_strong:
        reasons.append("imbalance_overridden_by_both_boundaries_strong")
        return True, reasons

    return False, reasons


def initial_eval(row: pd.Series, weak_z: float, rescue_z: float, balance_threshold: float, strong_z: float) -> SegmentEval:
    left = build_side_info(row, "left")
    right = build_side_info(row, "right")

    keep_ratio = compute_keep_ratio(row)
    keep_z, rescue_applies, rescue_side, z_fail_reasons, rescue_reasons = compute_keep_z(
        left, right, weak_z=weak_z, rescue_z=rescue_z
    )
    keep_balance, balance_ratio, balance_score, balance_source, bal_fail_reasons = compute_keep_balance(
        left, right, balance_threshold=balance_threshold, rescue_applies=rescue_applies
    )
    imbalance_override, override_reasons = compute_imbalance_override(
        row, left, right, strong_z=strong_z
    )

    if imbalance_override:
        keep_balance = True
        bal_fail_reasons = [x for x in bal_fail_reasons if x != "boundary_balance_below_threshold"]

    reasons_fail: List[str] = []
    reasons_rescue: List[str] = []
    reasons_note: List[str] = []

    if not keep_ratio:
        reasons_fail.append("ratio_below_threshold")

    reasons_fail.extend(z_fail_reasons)
    reasons_fail.extend(bal_fail_reasons)
    reasons_rescue.extend(rescue_reasons)

    if left.z_missing:
        reasons_note.append("left_boundary_z_NA")
    if right.z_missing:
        reasons_note.append("right_boundary_z_NA")
    if balance_source == "count":
        reasons_note.append("balance_used_support_counts")
    if rescue_applies:
        reasons_note.append(f"rescued_by_{rescue_side}_boundary")

    return SegmentEval(
        left=left,
        right=right,
        keep_ratio=keep_ratio,
        keep_z=keep_z,
        keep_balance=keep_balance,
        keep_neighbor=True,
        rescue_applies=rescue_applies,
        rescue_side=rescue_side,
        balance_source=balance_source,
        balance_ratio=balance_ratio,
        balance_score=balance_score,
        imbalance_override=imbalance_override,
        reasons_override=override_reasons,
        shared_boundary_suspect=False,
        neighbor_competitor="NA",
        reasons_fail=reasons_fail,
        reasons_rescue=reasons_rescue,
        reasons_note=reasons_note,
        reason_summary="",
    )


# ============================================================
# Neighbor competition
# ============================================================

def boundary_shared_strength(a: pd.Series, b: pd.Series) -> float:
    arz = get_numeric(a, "right_boundary_z")
    blz = get_numeric(b, "left_boundary_z")
    arz = 0.0 if math.isnan(arz) else max(arz, 0.0)
    blz = 0.0 if math.isnan(blz) else max(blz, 0.0)
    return max(arz, blz)


def outer_metric(ev: SegmentEval, side: str) -> float:
    s = ev.left if side == "left" else ev.right
    if ev.balance_source == "z":
        return s.z_eff
    return s.support_count


def apply_neighbor_competition(df: pd.DataFrame, evals: List[SegmentEval], weak_z: float) -> None:
    n = len(df)
    for i in range(n - 1):
        a = df.iloc[i]
        b = df.iloc[i + 1]

        if str(a.get("chrom")) != str(b.get("chrom")):
            continue

        shared = boundary_shared_strength(a, b)
        if shared < weak_z:
            continue

        eva = evals[i]
        evb = evals[i + 1]

        if not (eva.keep_ratio and eva.keep_z):
            continue
        if not (evb.keep_ratio and evb.keep_z):
            continue

        a_outer = outer_metric(eva, "left")
        b_outer = outer_metric(evb, "right")

        a_outer_bad = math.isnan(a_outer) or a_outer <= 0
        b_outer_bad = math.isnan(b_outer) or b_outer <= 0

        if a_outer_bad and (not b_outer_bad):
            eva.keep_neighbor = False
            eva.shared_boundary_suspect = True
            eva.neighbor_competitor = str(b.get("name", "NA"))
            eva.reasons_fail.append("shared_boundary_suspect")
            continue

        if b_outer_bad and (not a_outer_bad):
            evb.keep_neighbor = False
            evb.shared_boundary_suspect = True
            evb.neighbor_competitor = str(a.get("name", "NA"))
            evb.reasons_fail.append("shared_boundary_suspect")
            continue

        if (not eva.keep_balance) and (not eva.imbalance_override) and evb.keep_balance and (evb.balance_score > eva.balance_score * 1.5):
            eva.keep_neighbor = False
            eva.shared_boundary_suspect = True
            eva.neighbor_competitor = str(b.get("name", "NA"))
            eva.reasons_fail.append("shared_boundary_suspect")
            continue

        if (not evb.keep_balance) and (not evb.imbalance_override) and eva.keep_balance and (eva.balance_score > evb.balance_score * 1.5):
            evb.keep_neighbor = False
            evb.shared_boundary_suspect = True
            evb.neighbor_competitor = str(a.get("name", "NA"))
            evb.reasons_fail.append("shared_boundary_suspect")
            continue


def finalize_reasons(evals: List[SegmentEval]) -> None:
    for ev in evals:
        fail_parts = dedup_keep_order(ev.reasons_fail)
        rescue_parts = dedup_keep_order(ev.reasons_rescue)
        override_parts = dedup_keep_order(ev.reasons_override)
        note_parts = dedup_keep_order(ev.reasons_note)

        parts: List[str] = []

        if fail_parts:
            parts.append("fail=" + ",".join(fail_parts))
        if rescue_parts:
            parts.append("rescue=" + ",".join(rescue_parts))
        if override_parts:
            parts.append("override=" + ",".join(override_parts))
        if note_parts:
            parts.append("note=" + ",".join(note_parts))
        if ev.neighbor_competitor != "NA":
            parts.append(f"neighbor={ev.neighbor_competitor}")

        if not fail_parts and not rescue_parts and not override_parts and not note_parts and ev.neighbor_competitor == "NA":
            ev.reason_summary = "pass"
        elif not fail_parts:
            ev.reason_summary = "; ".join(parts)
        else:
            ev.reason_summary = "; ".join(parts)


def final_keep(ev: SegmentEval) -> bool:
    if not ev.keep_ratio:
        return False
    if not ev.keep_z:
        return False
    if not ev.keep_neighbor:
        return False
    if ev.keep_balance:
        return True
    if ev.imbalance_override:
        return True
    return False


# ============================================================
# Boundary nulling for neighbors of kept CNVs
# ============================================================

def kept_flags_from_evals(evals: List[SegmentEval]) -> List[bool]:
    return [final_keep(ev) for ev in evals]


def null_neighbor_shared_boundary(df: pd.DataFrame, idx: int, side: str) -> None:
    if side == "left":
        cols = [
            "left_boundary_name",
            "left_boundary_z",
            "left_split_reads",
            "left_disco_reads",
            "left_support_reads",
            "outer_left_split_reads",
            "outer_left_disco_reads",
            "outer_left_support_reads",
        ]
    elif side == "right":
        cols = [
            "right_boundary_name",
            "right_boundary_z",
            "right_split_reads",
            "right_disco_reads",
            "right_support_reads",
            "outer_right_split_reads",
            "outer_right_disco_reads",
            "outer_right_support_reads",
        ]
    else:
        return

    set_na(df, idx, cols)


def recompute_single_eval(df: pd.DataFrame, idx: int, weak_z: float, rescue_z: float, balance_threshold: float, strong_z: float) -> SegmentEval:
    row = df.iloc[idx]
    ev = initial_eval(
        row,
        weak_z=weak_z,
        rescue_z=rescue_z,
        balance_threshold=balance_threshold,
        strong_z=strong_z,
    )
    return ev


def apply_kept_neighbor_boundary_nulling(
    df: pd.DataFrame,
    evals: List[SegmentEval],
    weak_z: float,
    rescue_z: float,
    balance_threshold: float,
    strong_z: float,
) -> List[SegmentEval]:
    """
    For each kept segment, null the shared boundary evidence on immediate
    non-kept same-chromosome neighbors, then recompute those neighbor rows.
    """
    n = len(df)
    keeps = kept_flags_from_evals(evals)
    touched = set()

    for i in range(n):
        if not keeps[i]:
            continue

        chrom_i = str(df.iloc[i].get("chrom"))

        j = i - 1
        if j >= 0 and str(df.iloc[j].get("chrom")) == chrom_i and not keeps[j]:
            null_neighbor_shared_boundary(df, j, "right")
            touched.add(j)

        j = i + 1
        if j < n and str(df.iloc[j].get("chrom")) == chrom_i and not keeps[j]:
            null_neighbor_shared_boundary(df, j, "left")
            touched.add(j)

    for idx in sorted(touched):
        new_ev = recompute_single_eval(
            df=df,
            idx=idx,
            weak_z=weak_z,
            rescue_z=rescue_z,
            balance_threshold=balance_threshold,
            strong_z=strong_z,
        )
        new_ev.reasons_note.append("shared_boundary_nullified_by_neighbor")
        evals[idx] = new_ev

    apply_neighbor_competition(df, evals, weak_z=weak_z)
    finalize_reasons(evals)

    return evals


# ============================================================
# Output shaping
# ============================================================

def add_internal_backup_columns(df: pd.DataFrame) -> pd.DataFrame:
    if "_call_original" not in df.columns:
        df["_call_original"] = df["call"] if "call" in df.columns else "NA"
    if "_keep_original" not in df.columns:
        if "keep_suggested" in df.columns:
            df["_keep_original"] = df["keep_suggested"]
        else:
            df["_keep_original"] = "FALSE"

    drop_cols = [c for c in ["keep_suggested", "call"] if c in df.columns]
    if drop_cols:
        df = df.drop(columns=drop_cols)

    return df


def sort_segments(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["_orig_order"] = np.arange(len(df))

    for col in ["start", "end"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    sort_cols = [c for c in ["chrom", "start", "end", "name"] if c in df.columns]
    if sort_cols:
        df = df.sort_values(sort_cols, kind="mergesort").reset_index(drop=True)
    return df


def attach_columns(df: pd.DataFrame, evals: List[SegmentEval], sample: str, verbose: bool) -> pd.DataFrame:
    df["sample_name"] = sample
    df["keep_ratio"] = [bool_str(ev.keep_ratio) for ev in evals]
    df["keep_z"] = [bool_str(ev.keep_z) for ev in evals]
    df["keep_balance"] = [bool_str(ev.keep_balance) for ev in evals]
    df["keep_neighbor"] = [bool_str(ev.keep_neighbor) for ev in evals]
    df["balance_ratio"] = [ev.balance_ratio for ev in evals]
    df["imbalance_override"] = [bool_str(ev.imbalance_override) for ev in evals]
    df["heuristic_reason"] = [ev.reason_summary for ev in evals]

    if verbose:
        df["balance_source"] = [ev.balance_source for ev in evals]
        df["rescue_applies"] = [bool_str(ev.rescue_applies) for ev in evals]
        df["shared_boundary_suspect"] = [bool_str(ev.shared_boundary_suspect) for ev in evals]
        df["neighbor_competitor"] = [ev.neighbor_competitor for ev in evals]
        df["balance_score"] = [ev.balance_score for ev in evals]
        df["left_z_raw_used"] = [ev.left.z_raw for ev in evals]
        df["right_z_raw_used"] = [ev.right.z_raw for ev in evals]
        df["left_z_eff_used"] = [ev.left.z_eff for ev in evals]
        df["right_z_eff_used"] = [ev.right.z_eff for ev in evals]
        df["left_support_count_used"] = [ev.left.support_count for ev in evals]
        df["right_support_count_used"] = [ev.right.support_count for ev in evals]

    return df


def update_keep_and_call(df: pd.DataFrame, evals: List[SegmentEval]) -> pd.DataFrame:
    out_keep: List[str] = []
    out_call: List[str] = []

    for idx, row in df.iterrows():
        ev = evals[idx]
        keep = final_keep(ev)
        out_keep.append(bool_str(keep))

        if keep:
            if ev.rescue_applies:
                out_call.append("CNV_effect_rescued")
            elif ev.imbalance_override:
                out_call.append("CNV_effect_imbalance_override")
            else:
                out_call.append(str(row.get("_call_original", "CNV_effect")))
        else:
            if ev.shared_boundary_suspect:
                out_call.append("shared_boundary_suspect")
            elif not ev.keep_ratio:
                out_call.append("ratio_below_threshold")
            elif not ev.keep_z:
                out_call.append("boundary_support_fail")
            elif not ev.keep_balance and not ev.imbalance_override:
                out_call.append("boundary_balance_fail")
            elif not ev.keep_neighbor:
                out_call.append("neighbor_competition_fail")
            else:
                out_call.append("no_CNV_effect")

    df["keep_suggested"] = out_keep
    df["call"] = out_call
    return df


def trim_columns(df: pd.DataFrame, verbose: bool) -> pd.DataFrame:
    df = df.drop(columns=[c for c in ["_call_original", "_keep_original", "_orig_order"] if c in df.columns], errors="ignore")

    if verbose:
        return df

    drop_if_present = [
        "weak_ratio", "strong_ratio", "weak_z", "strong_z", "keep_mode",
        "left_split_reads", "left_disco_reads", "left_support_reads",
        "right_split_reads", "right_disco_reads", "right_support_reads",
        "fused", "fused_internal_boundaries", "fused_segment_names",
        "internal_left_boundary_support_mean", "internal_right_boundary_support_mean",
        "fused_left_boundaries_n", "fused_right_boundaries_n",
        "outer_left_split_reads", "outer_left_disco_reads", "outer_left_support_reads",
        "outer_right_split_reads", "outer_right_disco_reads", "outer_right_support_reads",
        "balance_source", "rescue_applies", "shared_boundary_suspect", "neighbor_competitor",
        "balance_score",
        "left_z_raw_used", "right_z_raw_used",
        "left_z_eff_used", "right_z_eff_used",
        "left_support_count_used", "right_support_count_used",
    ]

    df = df.drop(columns=[c for c in drop_if_present if c in df.columns], errors="ignore")

    preferred_tail = [
        "sample_name",
        "keep_ratio",
        "keep_z",
        "keep_balance",
        "keep_neighbor",
        "balance_ratio",
        "imbalance_override",
        "heuristic_reason",
        "keep_suggested",
        "call",
    ]

    existing_tail = [c for c in preferred_tail if c in df.columns]
    base_cols = [c for c in df.columns if c not in existing_tail]
    return df[base_cols + existing_tail]


# ============================================================
# Main
# ============================================================

def process(
    in_tsv: str,
    out_tsv: str,
    sample: str,
    balance: float,
    rescue_z: float,
    weak_z: float,
    strong_z: float,
    keep_suggested_only: bool,
    verbose: bool,
    debug: bool,
) -> None:
    df = read_tsv(in_tsv)
    if df.empty:
        write_tsv(df, out_tsv)
        return

    df = add_internal_backup_columns(df)
    df = coerce_numeric_columns(df)
    df = sort_segments(df)

    evals = [
        initial_eval(row, weak_z=weak_z, rescue_z=rescue_z, balance_threshold=balance, strong_z=strong_z)
        for _, row in df.iterrows()
    ]

    apply_neighbor_competition(df, evals, weak_z=weak_z)
    finalize_reasons(evals)

    evals = apply_kept_neighbor_boundary_nulling(
        df=df,
        evals=evals,
        weak_z=weak_z,
        rescue_z=rescue_z,
        balance_threshold=balance,
        strong_z=strong_z,
    )

    df = attach_columns(df, evals, sample=sample, verbose=verbose)
    df = update_keep_and_call(df, evals)

    if keep_suggested_only:
        df = df[df["keep_suggested"].astype(str).str.upper() == "TRUE"].copy()

    if "_orig_order" in df.columns:
        df = df.sort_values("_orig_order", kind="mergesort")

    df = trim_columns(df, verbose=verbose)
    write_tsv(df, out_tsv)

    if debug:
        kept = int((df["keep_suggested"].astype(str).str.upper() == "TRUE").sum()) if "keep_suggested" in df.columns else 0
        print(f"[cnv_finalize_suggested] wrote: {out_tsv}", file=sys.stderr)
        print(f"[cnv_finalize_suggested] rows_out: {len(df)}", file=sys.stderr)
        print(f"[cnv_finalize_suggested] kept_out: {kept}", file=sys.stderr)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=(
            "Post-finalization suggested-call refinement for CLOverCNV. "
            "Reads cnv_finalize_segments output, applies ratio/z/balance/neighbor logic, "
            "updates keep_suggested and call, and adds sample_name plus diagnostic columns."
        )
    )
    p.add_argument("--in-tsv", required=True, help="Input TSV from cnv_finalize_segments (plain or gzipped).")
    p.add_argument("--out-tsv", required=True, help="Output TSV path (plain or .gz).")
    p.add_argument("--sample", required=True, help="Sample name to populate into sample_name.")
    p.add_argument("--balance", type=float, required=True, help="Balance ratio threshold, typically min(sideA,sideB)/max(sideA,sideB).")
    p.add_argument("--rescue-z", type=float, required=True, help="Extreme one-sided rescue z threshold.")
    p.add_argument("--weak-z", type=float, required=True, help="Weak z threshold from upstream finalize.")
    p.add_argument("--strong-z", type=float, required=True, help="Strong z threshold from upstream finalize.")
    p.add_argument("--keep-suggested-only", action="store_true", help="Emit only rows with updated keep_suggested == TRUE.")
    p.add_argument("--verbose", action="store_true", help="Keep all diagnostic/helper columns.")
    p.add_argument("--debug", action="store_true", help="Verbose debug logging to stderr.")
    return p


def main() -> None:
    args = build_parser().parse_args()

    if not (0.0 <= args.balance <= 1.0):
        raise SystemExit("--balance must be between 0 and 1")
    if args.rescue_z < 0:
        raise SystemExit("--rescue-z must be >= 0")
    if args.weak_z < 0 or args.strong_z < 0:
        raise SystemExit("--weak-z and --strong-z must be >= 0")
    if args.strong_z < args.weak_z:
        raise SystemExit("--strong-z must be >= --weak-z")

    process(
        in_tsv=args.in_tsv,
        out_tsv=args.out_tsv,
        sample=args.sample,
        balance=args.balance,
        rescue_z=args.rescue_z,
        weak_z=args.weak_z,
        strong_z=args.strong_z,
        keep_suggested_only=args.keep_suggested_only,
        verbose=args.verbose,
        debug=args.debug,
    )


if __name__ == "__main__":
    main()
