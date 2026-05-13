#!/usr/bin/env python3
import argparse
import gzip
import math
import os
import tempfile
from typing import Optional, List, Dict, Tuple, Set

import numpy as np
import pandas as pd
import pysam


# ============================================================
# IO helpers
# ============================================================

def _is_gzip_by_magic(path: str) -> bool:
    try:
        with open(path, "rb") as fh:
            return fh.read(2) == b"\x1f\x8b"
    except Exception:
        return False


def read_tsv_any(path: str, header="infer") -> pd.DataFrame:
    """
    Read TSV. If file ends with .gz but isn't actually gzipped, detect by magic bytes.
    """
    is_gz = _is_gzip_by_magic(path)
    if is_gz:
        return pd.read_csv(path, sep="\t", header=header, compression="gzip")
    return pd.read_csv(path, sep="\t", header=header)


def write_tsv_any(df: pd.DataFrame, path: str) -> None:
    """
    Write TSV; if path ends with .gz, gzip it. Always write NA explicitly for missing.
    """
    if path.endswith(".gz"):
        with gzip.open(path, "wt") as f:
            df.to_csv(f, sep="\t", index=False, na_rep="NA")
    else:
        df.to_csv(path, sep="\t", index=False, na_rep="NA")


def write_bedgraph_any(df_bg: pd.DataFrame, path: str) -> None:
    """
    Write bedGraph (no header). Always write NA explicitly for missing.
    """
    if path.endswith(".gz"):
        with gzip.open(path, "wt") as f:
            df_bg.to_csv(f, sep="\t", header=False, index=False, na_rep="NA")
    else:
        df_bg.to_csv(path, sep="\t", header=False, index=False, na_rep="NA")


def _atomic_write_no_follow(write_fn, df: pd.DataFrame, out_path: str) -> None:
    """
    Write output without following an existing symlink at out_path.
    Writes to temp in same dir, then os.replace (replaces symlink itself if present).
    """
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


# ============================================================
# Segment parsing
# ============================================================

def _maybe_infer_ratio_col(df: pd.DataFrame, y_col: int) -> Optional[int]:
    """
    If there is a column immediately after y_col that looks like ratio (positive, ~2^y),
    use it; else return None (we'll compute ratio).
    """
    if y_col + 1 >= df.shape[1]:
        return None

    y = pd.to_numeric(df.iloc[:, y_col], errors="coerce").to_numpy(dtype=float)
    r = pd.to_numeric(df.iloc[:, y_col + 1], errors="coerce").to_numpy(dtype=float)

    finite = np.isfinite(y) & np.isfinite(r)
    if finite.sum() < 50:
        return None

    yy = y[finite]
    rr = r[finite]

    if np.mean(rr > 0) < 0.95:
        return None

    with np.errstate(divide="ignore", invalid="ignore"):
        err = np.nanmedian(np.abs(np.log2(rr) - yy))
    if not np.isfinite(err):
        return None

    if err <= 0.25:
        return y_col + 1

    return None


def _infer_y_col(df: pd.DataFrame, candidate_cols: List[int], n_check: int = 200) -> int:
    """
    Infer which column is y (log2 ratio) from a headerless segments BED-like table.
    Prefer columns with both negative and positive values.
    """
    n = min(n_check, len(df))
    if n == 0:
        raise RuntimeError("Empty segments file; cannot infer y column.")

    best_col = None
    best_score = -1e18

    for c in candidate_cols:
        vals = pd.to_numeric(df.iloc[:n, c], errors="coerce").to_numpy(dtype=float)
        finite = np.isfinite(vals)
        if finite.sum() < max(10, int(0.5 * n)):
            continue

        v = vals[finite]

        has_neg = np.any(v < 0)
        has_pos = np.any(v > 0)
        sign_bonus = 2.0 if (has_neg and has_pos) else (-2.0)

        uniq = len(np.unique(np.round(v, 6)))
        uniq_bonus = min(1.0, uniq / 20.0)

        frac_in = np.mean((v >= -20.0) & (v <= 20.0))
        med_abs = np.median(np.abs(v))

        score = (2.0 * frac_in) + (1.0 / (1.0 + med_abs)) + sign_bonus + uniq_bonus

        if score > best_score:
            best_score = score
            best_col = c

    if best_col is None:
        raise RuntimeError(
            "Could not infer y (log2 ratio) column. "
            "Your segments.bed column layout is unusual-please inspect columns."
        )
    return int(best_col)


def read_segments_bed(path: str) -> pd.DataFrame:
    """
    segments.bed(.gz) expected columns (0-based):
      0 chrom
      1 start
      2 end
      3 name
      6 y (log2 ratio)                 (if >=7 cols)
      7 ratio (optional; computed if missing)
      9 n_probes (optional; if present)
    """
    df = read_tsv_any(path, header=None)

    if df.shape[1] < 4:
        raise RuntimeError("segments.bed must have >= 4 columns (chrom,start,end,name).")

    out = pd.DataFrame({
        "chrom": df.iloc[:, 0].astype(str),
        "start": pd.to_numeric(df.iloc[:, 1], errors="coerce"),
        "end": pd.to_numeric(df.iloc[:, 2], errors="coerce"),
        "name": df.iloc[:, 3].astype(str),
    })

    if out["start"].isna().any() or out["end"].isna().any():
        raise RuntimeError("segments.bed has non-numeric start/end in columns 1/2-unexpected format.")

    out["start"] = out["start"].astype(int)
    out["end"] = out["end"].astype(int)

    if df.shape[1] >= 7:
        y_col = 6
        out["y"] = pd.to_numeric(df.iloc[:, y_col], errors="coerce").astype(float)

        if df.shape[1] >= 8:
            out["ratio"] = pd.to_numeric(df.iloc[:, 7], errors="coerce").astype(float)
        else:
            out["ratio"] = np.power(2.0, out["y"].to_numpy(dtype=float))

        if df.shape[1] >= 10:
            out["n_probes"] = pd.to_numeric(df.iloc[:, 9], errors="coerce").astype("Int64")
        else:
            out["n_probes"] = pd.Series([pd.NA] * len(out), dtype="Int64")

        return out

    candidate_cols = list(range(4, df.shape[1]))
    y_col = _infer_y_col(df, candidate_cols)
    out["y"] = pd.to_numeric(df.iloc[:, y_col], errors="coerce").astype(float)

    ratio_col = _maybe_infer_ratio_col(df, y_col)
    if ratio_col is not None:
        out["ratio"] = pd.to_numeric(df.iloc[:, ratio_col], errors="coerce").astype(float)
    else:
        out["ratio"] = np.power(2.0, out["y"].to_numpy(dtype=float))

    out["n_probes"] = pd.Series([pd.NA] * len(out), dtype="Int64")
    return out


# ============================================================
# Boundary support TSV
# ============================================================

def read_boundary_support(path: str) -> pd.DataFrame:
    """
    boundary_support.tsv expected at least:
      chrom, name, combined_z_chr
    """
    df = read_tsv_any(path, header="infer")
    required = {"chrom", "name", "combined_z_chr"}
    missing = required - set(df.columns)
    if missing:
        raise RuntimeError(f"boundary_support missing required columns: {sorted(missing)}")
    df["chrom"] = df["chrom"].astype(str)
    df["name"] = df["name"].astype(str)
    df["combined_z_chr"] = pd.to_numeric(df["combined_z_chr"], errors="coerce").astype(float)
    return df


# ============================================================
# Circular contigs helpers
# ============================================================

def _read_lines_list(path: str) -> List[str]:
    out: List[str] = []
    if not path:
        return out
    with open(path, "r") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            out.append(s)
    return out


def _norm_contig(s: str) -> str:
    s = (s or "").strip()
    if s.lower().startswith("chr"):
        return s[3:]
    return s


def parse_circular_contigs(csv: str, path: str) -> Set[str]:
    """
    Store BOTH raw and normalized names to maximize matching:
      - raw: as provided
      - norm: stripped + optional chr-prefix removed
    """
    names: Set[str] = set()

    def _add(name: str):
        name = (name or "").strip()
        if not name:
            return
        names.add(name)
        names.add(_norm_contig(name))

    if csv:
        for part in csv.split(","):
            _add(part)

    if path:
        for n in _read_lines_list(path):
            _add(n)

    return names


def is_circular_contig(chrom: str, circular_contigs: Set[str]) -> bool:
    """
    Match against both raw chrom and normalized chrom.
    """
    if not chrom:
        return False
    chrom_raw = str(chrom).strip()
    chrom_norm = _norm_contig(chrom_raw)
    return (chrom_raw in circular_contigs) or (chrom_norm in circular_contigs)


def _with_chr(chrom: str) -> str:
    c = (chrom or "").strip()
    if not c:
        return c
    return c if c.lower().startswith("chr") else ("chr" + c)


def _strip_chr(chrom: str) -> str:
    c = (chrom or "").strip()
    if c.lower().startswith("chr"):
        return c[3:]
    return c


def resolve_bam_contig_name(bam: Optional[pysam.AlignmentFile], chrom: str) -> Optional[str]:
    """
    Return a contig name that exists in bam.references, trying common variants:
      - exact
      - add/remove chr
      - normalized match (_norm_contig equality)
    """
    if bam is None or not chrom:
        return None

    chrom = str(chrom).strip()
    if not chrom:
        return None

    if chrom in bam.references:
        return chrom

    c1 = _with_chr(chrom)
    if c1 in bam.references:
        return c1

    c2 = _strip_chr(chrom)
    if c2 in bam.references:
        return c2

    want = _norm_contig(chrom)
    for r in bam.references:
        if _norm_contig(r) == want:
            return r

    return None


def build_contig_length_map(bam: Optional[pysam.AlignmentFile]) -> Dict[str, int]:
    """
    Build a length map that supports:
      - exact name keys
      - normalized keys (_norm_contig)
    If multiple refs map to same normalized key, keep the MAX length.
    """
    out: Dict[str, int] = {}
    if bam is None:
        return out
    try:
        for r, L in zip(bam.references, bam.lengths):
            if not r:
                continue
            L = int(L)
            if L <= 0:
                continue
            out[r] = L
            nr = _norm_contig(r)
            if nr:
                out[nr] = max(out.get(nr, 0), L)
    except Exception:
        pass
    return out


def get_contig_len(chrom: str, contig_lengths: Dict[str, int]) -> Optional[int]:
    """
    Robust contig length lookup using exact, norm, chr add/remove.
    """
    if not chrom:
        return None
    c = str(chrom).strip()
    if not c:
        return None

    L = contig_lengths.get(c, 0)
    if L and int(L) > 0:
        return int(L)

    n = _norm_contig(c)
    L = contig_lengths.get(n, 0)
    if L and int(L) > 0:
        return int(L)

    L = contig_lengths.get(_with_chr(c), 0)
    if L and int(L) > 0:
        return int(L)
    L = contig_lengths.get(_strip_chr(c), 0)
    if L and int(L) > 0:
        return int(L)

    L = contig_lengths.get(_norm_contig(_with_chr(c)), 0)
    if L and int(L) > 0:
        return int(L)
    L = contig_lengths.get(_norm_contig(_strip_chr(c)), 0)
    if L and int(L) > 0:
        return int(L)

    return None


def _clip_window(L: Optional[int], start: int, end: int) -> Tuple[int, int]:
    if start < 0:
        start = 0
    if end < start:
        end = start
    if L is not None and L > 0:
        if start > L:
            start = L
        if end > L:
            end = L
        if end < start:
            end = start
    return int(start), int(end)


def _wrap_windows_0based_halfopen(L: int, flank: int) -> List[Tuple[int, int]]:
    """
    For circular boundary at end<->start, return two windows:
      [0, flank) and [L-flank, L)
    """
    if L <= 0 or flank <= 0:
        return []
    a0, a1 = 0, min(L, flank)
    b0, b1 = max(0, L - flank), L
    wins = [(a0, a1), (b0, b1)]
    uniq: List[Tuple[int, int]] = []
    for w in wins:
        if w not in uniq and w[1] > w[0]:
            uniq.append(w)
    return uniq


# ============================================================
# BAM counting helpers
# ============================================================

def _open_bam(path: str) -> Optional[pysam.AlignmentFile]:
    if not path:
        return None
    if not os.path.exists(path):
        raise RuntimeError(f"BAM not found: {path}")
    return pysam.AlignmentFile(path, "rb")


def count_reads_in_interval(
    bam: pysam.AlignmentFile,
    chrom: str,
    start: int,
    end: int,
    read1_only: bool = True,
    drop_dup: bool = False,
    include_supp: bool = False,
    include_secondary: bool = False,
) -> int:
    """
    Counts ALIGNMENT RECORDS overlapping [start,end).
    Kept for seg_primary_reads.
    """
    if start < 0:
        start = 0
    if end < start:
        end = start

    bam_chrom = resolve_bam_contig_name(bam, chrom)
    if bam_chrom is None:
        return 0

    n = 0
    for r in bam.fetch(bam_chrom, start, end):
        if r.is_unmapped:
            continue
        if (not include_secondary) and r.is_secondary:
            continue
        if (not include_supp) and r.is_supplementary:
            continue
        if drop_dup and r.is_duplicate:
            continue

        if read1_only:
            if r.is_paired and (not r.is_read1):
                continue
        n += 1
    return n


def _count_unique_qnames_in_windows(
    bam: pysam.AlignmentFile,
    chrom: str,
    windows: List[Tuple[int, int]],
    read1_only: bool = True,
    drop_dup: bool = False,
    include_supp: bool = False,
    include_secondary: bool = False,
) -> int:
    """
    Counts UNIQUE QNAMEs overlapping the union of [start,end) windows.
    De-duplicates across windows by QNAME.
    """
    if not windows:
        return 0

    bam_chrom = resolve_bam_contig_name(bam, chrom)
    if bam_chrom is None:
        return 0

    seen = set()
    for (start, end) in windows:
        if start < 0:
            start = 0
        if end < start:
            end = start

        for r in bam.fetch(bam_chrom, start, end):
            if r.is_unmapped:
                continue
            if (not include_secondary) and r.is_secondary:
                continue
            if (not include_supp) and r.is_supplementary:
                continue
            if drop_dup and r.is_duplicate:
                continue

            if read1_only:
                if r.is_paired and (not r.is_read1):
                    continue

            seen.add(r.query_name)

    return len(seen)


def count_unique_qnames_for_boundary(
    bam: pysam.AlignmentFile,
    chrom: str,
    boundary_pos_0based: Optional[int],
    anchor_pos_0based: Optional[int],
    flank: int,
    circular_contigs: Set[str],
    contig_lengths: Dict[str, int],
    read1_only: bool = True,
    drop_dup: bool = False,
    include_supp: bool = False,
    include_secondary: bool = False,
    debug: bool = False,
    debug_tag: str = "",
) -> Optional[int]:
    """
    Return unique-QNAME count for a boundary.

    Rules:
      - If boundary_pos_0based is not None: use window centered on boundary: [pos-flank, pos+flank)
      - If boundary_pos_0based is None (terminal boundary):
          * If chrom is circular: wrap windows ([0,flank), [L-flank,L))
          * Else (linear): use window centered on anchor_pos_0based (segment start/end):
                [anchor-flank, anchor+flank)
            If anchor is None, return None.
    """
    bam_chrom = resolve_bam_contig_name(bam, chrom)
    if bam_chrom is None:
        if debug:
            print(f"[cnv_finalize_segments][DEBUG]{debug_tag} chrom={chrom} -> bam_chrom=None (counts=0)")
        return 0

    L_opt = get_contig_len(bam_chrom, contig_lengths)
    if debug:
        print(f"[cnv_finalize_segments][DEBUG]{debug_tag} chrom={chrom} bam_chrom={bam_chrom} L_opt={L_opt} "
              f"boundary_pos_0={boundary_pos_0based} anchor_pos_0={anchor_pos_0based} flank={flank}")

    if boundary_pos_0based is None:
        if is_circular_contig(chrom, circular_contigs):
            if L_opt is None or int(L_opt) <= 0:
                if debug:
                    print(f"[cnv_finalize_segments][DEBUG]{debug_tag} TERMINAL+CIRCULAR but no contig length; returning None")
                return None

            wins = _wrap_windows_0based_halfopen(int(L_opt), flank)
            if debug:
                print(f"[cnv_finalize_segments][DEBUG]{debug_tag} TERMINAL+CIRCULAR wrap windows (0-based half-open): {wins}")
            return _count_unique_qnames_in_windows(
                bam, bam_chrom, wins,
                read1_only=read1_only, drop_dup=drop_dup,
                include_supp=include_supp, include_secondary=include_secondary
            )

        if anchor_pos_0based is None:
            if debug:
                print(f"[cnv_finalize_segments][DEBUG]{debug_tag} TERMINAL+LINEAR but anchor=None; returning None")
            return None

        pos = int(anchor_pos_0based)
        w0, w1 = _clip_window(L_opt, pos - flank, pos + flank)
        wins = [(w0, w1)]
        if debug:
            print(f"[cnv_finalize_segments][DEBUG]{debug_tag} TERMINAL+LINEAR windows (0-based half-open): {wins}")
        return _count_unique_qnames_in_windows(
            bam, bam_chrom, wins,
            read1_only=read1_only, drop_dup=drop_dup,
            include_supp=include_supp, include_secondary=include_secondary
        )

    pos = int(boundary_pos_0based)
    w0, w1 = _clip_window(L_opt, pos - flank, pos + flank)
    wins = [(w0, w1)]
    if debug:
        print(f"[cnv_finalize_segments][DEBUG]{debug_tag} INTERNAL windows (0-based half-open): {wins}")
    return _count_unique_qnames_in_windows(
        bam, bam_chrom, wins,
        read1_only=read1_only, drop_dup=drop_dup,
        include_supp=include_supp, include_secondary=include_secondary
    )


# ============================================================
# Generic CNV helpers
# ============================================================

def min2(a, b):
    vals = []
    try:
        if a is not None and math.isfinite(float(a)):
            vals.append(float(a))
    except Exception:
        pass
    try:
        if b is not None and math.isfinite(float(b)):
            vals.append(float(b))
    except Exception:
        pass
    if not vals:
        return float("nan")
    return float(min(vals))


def z_tier(z: float, weak_z: float, strong_z: float) -> str:
    try:
        z = float(z)
    except Exception:
        return "NA"

    if not math.isfinite(z):
        return "NA"

    try:
        weak_z = float(weak_z)
        strong_z = float(strong_z)
    except Exception:
        return "NA"

    if z >= strong_z:
        return "strong"
    if z >= weak_z:
        return "weak"
    return "low"


def ratio_tier(r: float, weak_ratio: float, strong_ratio: float) -> str:
    try:
        r = float(r)
    except Exception:
        return "NA"

    if not math.isfinite(r) or r <= 0:
        return "NA"

    try:
        weak_ratio = float(weak_ratio)
        strong_ratio = float(strong_ratio)
    except Exception:
        return "NA"

    if not math.isfinite(weak_ratio) or not math.isfinite(strong_ratio):
        return "NA"
    if weak_ratio <= 0 or strong_ratio <= 0:
        return "NA"

    if r >= strong_ratio or r <= (1.0 / strong_ratio):
        return "strong"
    if r >= weak_ratio or r <= (1.0 / weak_ratio):
        return "weak"
    return "low"


def call_direction(r: float) -> str:
    try:
        r = float(r)
    except Exception:
        return "NA"

    if not math.isfinite(r) or r <= 0:
        return "NA"
    if r > 1.0:
        return "gain"
    if r < 1.0:
        return "loss"
    return "neutral"


def _sum_int(vals):
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
            pass
    return int(s) if any_ok else pd.NA


def _ensure_fusion_schema(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()

    if "fused" not in out.columns:
        out["fused"] = False
    out["fused"] = out["fused"].fillna(False).astype(bool)

    want_int = [
        "fused_n_segments",
        "fused_internal_boundaries",
        "fused_left_boundaries_n",
        "fused_right_boundaries_n",
        "outer_left_split_reads",
        "outer_left_disco_reads",
        "outer_left_support_reads",
        "outer_right_split_reads",
        "outer_right_disco_reads",
        "outer_right_support_reads",
    ]
    want_float = [
        "internal_left_boundary_support_mean",
        "internal_right_boundary_support_mean",
    ]
    want_str = [
        "fused_segment_names",
    ]

    for c in want_int:
        if c not in out.columns:
            out[c] = pd.NA
    for c in want_float:
        if c not in out.columns:
            out[c] = np.nan
    for c in want_str:
        if c not in out.columns:
            out[c] = pd.NA

    for c in want_int:
        out[c] = pd.to_numeric(out[c], errors="coerce").astype("Int64")
    for c in want_float:
        out[c] = pd.to_numeric(out[c], errors="coerce").astype(float)

    return out


# ============================================================
# Fusion helpers
# ============================================================

def _finite_float(x, default=np.nan):
    try:
        x = float(x)
        return x if math.isfinite(x) else float(default)
    except Exception:
        return float(default)


def _seg_span(row) -> int:
    try:
        return max(0, int(row["end"]) - int(row["start"]))
    except Exception:
        return 0


def _seg_cn(row) -> float:
    cn = _finite_float(row.get("cn", np.nan), default=np.nan)
    if math.isfinite(cn) and cn > 0:
        return cn
    rr = _finite_float(row.get("ratio", np.nan), default=np.nan)
    return rr if math.isfinite(rr) and rr > 0 else 1.0


def _seg_ratio_tier(row) -> str:
    return str(row.get("ratio_tier", "NA"))


def _seg_direction(row) -> str:
    return str(row.get("direction", "NA"))


def _seg_support_tier(row) -> str:
    return str(row.get("support_tier", "NA"))


def _seg_passes_boundary_z(row, weak_z_value: float) -> bool:
    lz = _finite_float(row.get("left_boundary_z", np.nan), default=np.nan)
    rz = _finite_float(row.get("right_boundary_z", np.nan), default=np.nan)

    if math.isfinite(lz) and lz >= weak_z_value:
        return True
    if math.isfinite(rz) and rz >= weak_z_value:
        return True

    ratio_z = _finite_float(row.get("ratio_z", np.nan), default=np.nan)
    if not math.isfinite(ratio_z):
        ratio_z = _finite_float(row.get("seg_conf_z", np.nan), default=np.nan)

    return bool(math.isfinite(ratio_z) and ratio_z >= weak_z_value)


def _block_passes_ratio_z(df_block: pd.DataFrame, weak_z_value: float) -> bool:
    """
    Pass if ANY member of the block has a boundary meeting/exceeding weak_z.
    """
    if df_block is None or df_block.empty:
        return False
    for _, row in df_block.iterrows():
        if _seg_passes_boundary_z(row, weak_z_value):
            return True
    return False


def _is_bridgeable_opposite_insert(
    left_row: pd.Series,
    middle_row: pd.Series,
    right_row: pd.Series,
    max_gap_value: int,
    bridge_max_cap_mult: float = 5.0,
    forbid_strong_opposite: bool = True,
    forbid_supported_opposite: bool = False,
) -> bool:
    """
    Allow:
        gain ... loss ... gain
    or
        loss ... gain ... loss

    only if the middle opposite-direction segment is short enough.

    Threshold:
        allowable_insert = max_gap * mean(cn_left, cn_right)

    capped at:
        max_gap * bridge_max_cap_mult
    """
    d_left = _seg_direction(left_row)
    d_mid = _seg_direction(middle_row)
    d_right = _seg_direction(right_row)

    if d_left not in {"gain", "loss"}:
        return False
    if d_right != d_left:
        return False
    if d_mid not in {"gain", "loss"}:
        return False
    if d_mid == d_left:
        return False

    if forbid_strong_opposite and _seg_ratio_tier(middle_row) == "strong":
        return False

    if forbid_supported_opposite and _seg_support_tier(middle_row) == "strong":
        return False

    gap_left = int(middle_row["start"]) - int(left_row["end"])
    gap_right = int(right_row["start"]) - int(middle_row["end"])

    if gap_left > int(max_gap_value):
        return False
    if gap_right > int(max_gap_value):
        return False

    cn_left = _seg_cn(left_row)
    cn_right = _seg_cn(right_row)
    flank_cn_mean = np.nanmean([cn_left, cn_right])

    if not math.isfinite(flank_cn_mean) or flank_cn_mean <= 0:
        flank_cn_mean = 1.0

    allowable_insert = float(max_gap_value) * float(flank_cn_mean)
    allowable_insert = min(
        allowable_insert,
        float(max_gap_value) * float(bridge_max_cap_mult)
    )

    insert_size = _seg_span(middle_row)
    return bool(insert_size <= allowable_insert)


def fuse_segments_preserve_support(df: pd.DataFrame, max_gap: int = 0) -> pd.DataFrame:
    if df.empty:
        out0 = df.copy()
        out0 = _ensure_fusion_schema(out0)
        return out0

    out = df.copy()

    out["start"] = pd.to_numeric(out["start"], errors="coerce").astype(int)
    out["end"] = pd.to_numeric(out["end"], errors="coerce").astype(int)
    out["ratio"] = pd.to_numeric(out["ratio"], errors="coerce")
    out["y"] = pd.to_numeric(out["y"], errors="coerce")
    out["cn"] = pd.to_numeric(out.get("cn", np.nan), errors="coerce")

    eligible = (
        (out.get("ratio_tier", "low") != "low") &
        (out.get("direction", "NA").isin(["gain", "loss"]))
    )

    def _weights(df_block: pd.DataFrame) -> np.ndarray:
        wp = pd.to_numeric(df_block.get("n_probes", np.nan), errors="coerce").to_numpy(dtype=float)
        span = (
            pd.to_numeric(df_block["end"], errors="coerce") -
            pd.to_numeric(df_block["start"], errors="coerce")
        ).to_numpy(dtype=float)
        w = np.where(np.isfinite(wp) & (wp > 0), wp, span)
        w = np.where(np.isfinite(w) & (w > 0), w, 1.0)
        return w

    def _wavg(vals: np.ndarray, w: np.ndarray) -> float:
        m = np.isfinite(vals) & np.isfinite(w) & (w > 0)
        if m.sum() == 0:
            return float("nan")
        return float(np.sum(vals[m] * w[m]) / np.sum(w[m]))

    def _wmean(vals: List[float], wts: List[float]) -> float:
        vv = np.asarray(vals, dtype=float)
        ww = np.asarray(wts, dtype=float)
        m = np.isfinite(vv) & np.isfinite(ww) & (ww > 0)
        if m.sum() == 0:
            return float("nan")
        return float(np.sum(vv[m] * ww[m]) / np.sum(ww[m]))

    BRIDGE_MAX_CAP_MULT = 5.0
    FORBID_STRONG_OPPOSITE = True
    FORBID_SUPPORTED_OPPOSITE = False

    fused_rows = []
    fused_member_idx = set()

    for chrom, sub in out.groupby("chrom", sort=False):
        sub = sub.sort_values(["start", "end"]).copy()
        sub["_orig_idx"] = sub.index
        sub = sub.reset_index(drop=True)

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

            seed = sub.iloc[i]
            direction = _seg_direction(seed)
            if direction not in {"gain", "loss"}:
                i += 1
                continue

            block_idx = [i]
            j = i + 1
            used_bridge = False

            while j < len(sub):
                cur = sub.iloc[block_idx[-1]]
                nxt = sub.iloc[j]

                if not bool(nxt["_eligible"]):
                    break

                nxt_dir = _seg_direction(nxt)

                if nxt_dir == direction:
                    gap = int(nxt["start"]) - int(cur["end"])
                    if gap > int(max_gap):
                        break
                    block_idx.append(j)
                    j += 1
                    continue

                if (not used_bridge) and (j + 1 < len(sub)):
                    mid = sub.iloc[j]
                    right = sub.iloc[j + 1]

                    if (
                        bool(mid["_eligible"]) and
                        bool(right["_eligible"]) and
                        _is_bridgeable_opposite_insert(
                            left_row=cur,
                            middle_row=mid,
                            right_row=right,
                            max_gap_value=int(max_gap),
                            bridge_max_cap_mult=BRIDGE_MAX_CAP_MULT,
                            forbid_strong_opposite=FORBID_STRONG_OPPOSITE,
                            forbid_supported_opposite=FORBID_SUPPORTED_OPPOSITE,
                        )
                    ):
                        block_idx.extend([j, j + 1])
                        used_bridge = True
                        j += 2
                        continue

                break

            block = sub.iloc[block_idx].copy()

            base = block.iloc[0].to_dict()
            weak_ratio = float(base.get("weak_ratio", np.nan))
            strong_ratio = float(base.get("strong_ratio", np.nan))
            weak_z = float(base.get("weak_z", np.nan))
            strong_z = float(base.get("strong_z", np.nan))
            keep_mode = str(base.get("keep_mode", "ratio_only"))

            if len(block) > 1 and _block_passes_ratio_z(block, weak_z_value=weak_z):
                fused_member_idx.update(block["_orig_idx"].tolist())
                w = _weights(block)

                new_start = int(block["start"].min())
                new_end = int(block["end"].max())

                ratio_vals = pd.to_numeric(block["ratio"], errors="coerce").to_numpy(dtype=float)
                y_vals = pd.to_numeric(block["y"], errors="coerce").to_numpy(dtype=float)
                cn_vals = pd.to_numeric(block.get("cn", np.nan), errors="coerce").to_numpy(dtype=float)
                nprobes_vals = pd.to_numeric(block.get("n_probes", np.nan), errors="coerce").to_numpy(dtype=float)

                fused_ratio = _wavg(ratio_vals, w)
                fused_y = _wavg(y_vals, w)
                fused_cn = _wavg(cn_vals, w)
                fused_nprobes = float(np.nansum(nprobes_vals)) if np.isfinite(nprobes_vals).any() else float("nan")

                left_name = str(block.iloc[0].get("left_boundary_name", "NA"))
                right_name = str(block.iloc[-1].get("right_boundary_name", block.iloc[-1].get("name", "NA")))

                left_z = pd.to_numeric(block.iloc[0].get("left_boundary_z", np.nan), errors="coerce")
                right_z = pd.to_numeric(block.iloc[-1].get("right_boundary_z", np.nan), errors="coerce")

                conf_method = str(base.get("conf_method", base.get("confidence_method", "min"))).lower()
                if conf_method not in {"min", "mean", "max"}:
                    conf_method = "min"

                left_z_val = float(left_z) if math.isfinite(float(left_z)) else float("nan")
                right_z_val = float(right_z) if math.isfinite(float(right_z)) else float("nan")
                z_pair = [z for z in [left_z_val, right_z_val] if math.isfinite(z)]

                if len(z_pair) == 0:
                    fused_seg_conf_z = float("nan")
                elif conf_method == "mean":
                    fused_seg_conf_z = float(np.mean(z_pair))
                elif conf_method == "max":
                    fused_seg_conf_z = float(np.max(z_pair))
                else:
                    fused_seg_conf_z = float(np.min(z_pair))

                mid = (new_start + new_end) / 2.0
                segw = _weights(block)

                left_vals, left_wts = [], []
                right_vals, right_wts = [], []

                for k in range(len(block) - 1):
                    bpos = int(block.iloc[k]["end"])

                    a = pd.to_numeric(block.iloc[k].get("right_support_reads", np.nan), errors="coerce")
                    b2 = pd.to_numeric(block.iloc[k + 1].get("left_support_reads", np.nan), errors="coerce")
                    aval = float(a) if math.isfinite(float(a)) else float("nan")
                    bval = float(b2) if math.isfinite(float(b2)) else float("nan")

                    if math.isfinite(aval) and math.isfinite(bval):
                        sval = 0.5 * (aval + bval)
                    elif math.isfinite(aval):
                        sval = aval
                    elif math.isfinite(bval):
                        sval = bval
                    else:
                        continue

                    w_k = 0.5 * (float(segw[k]) + float(segw[k + 1]))
                    if not math.isfinite(w_k) or w_k <= 0:
                        w_k = 1.0

                    if bpos <= mid:
                        left_vals.append(float(sval))
                        left_wts.append(float(w_k))
                    else:
                        right_vals.append(float(sval))
                        right_wts.append(float(w_k))

                fused_internal_left_support_mean = _wmean(left_vals, left_wts) if left_vals else float("nan")
                fused_internal_right_support_mean = _wmean(right_vals, right_wts) if right_vals else float("nan")

                fused_left_boundaries_n = int(len(left_vals))
                fused_right_boundaries_n = int(len(right_vals))

                names = [str(x) for x in block["name"].astype(str).tolist()]
                seen = set()
                names_uniq = []
                for nm in names:
                    if nm not in seen:
                        names_uniq.append(nm)
                        seen.add(nm)
                fused_names_str = ",".join(names_uniq)

                fused_rtier = ratio_tier(float(fused_ratio), weak_ratio, strong_ratio)
                fused_ztier = z_tier(float(fused_seg_conf_z), weak_z, strong_z)

                ratio_keep = fused_rtier in ["weak", "strong"]
                z_keep = fused_ztier in ["weak", "strong"]

                if keep_mode == "ratio_only":
                    fused_keep = ratio_keep
                elif keep_mode == "z_only":
                    fused_keep = z_keep
                else:
                    fused_keep = ratio_keep and z_keep

                if not ratio_keep:
                    fused_call = "no_CNV_effect"
                else:
                    fused_call = f"CNV_{direction}_{fused_rtier}Effect_{fused_ztier}Support"
                fused_call = fused_call + "_FUSED"

                fused_ratio_z = float(fused_seg_conf_z) if math.isfinite(float(fused_seg_conf_z)) else float("nan")
                fused_seg_len_bp = int(new_end) - int(new_start)

                fused_seg_primary_reads = _sum_int(block["seg_primary_reads"].tolist()) if "seg_primary_reads" in block.columns else pd.NA

                fused_left_split_reads = block.iloc[0].get("left_split_reads", pd.NA) if "left_split_reads" in block.columns else pd.NA
                fused_left_disco_reads = block.iloc[0].get("left_disco_reads", pd.NA) if "left_disco_reads" in block.columns else pd.NA
                fused_right_split_reads = block.iloc[-1].get("right_split_reads", pd.NA) if "right_split_reads" in block.columns else pd.NA
                fused_right_disco_reads = block.iloc[-1].get("right_disco_reads", pd.NA) if "right_disco_reads" in block.columns else pd.NA

                fused_left_support_reads = pd.NA
                fused_right_support_reads = pd.NA
                try:
                    if pd.notna(fused_left_split_reads) and pd.notna(fused_left_disco_reads):
                        fused_left_support_reads = int(fused_left_split_reads) + int(fused_left_disco_reads)
                except Exception:
                    fused_left_support_reads = pd.NA
                try:
                    if pd.notna(fused_right_split_reads) and pd.notna(fused_right_disco_reads):
                        fused_right_support_reads = int(fused_right_split_reads) + int(fused_right_disco_reads)
                except Exception:
                    fused_right_support_reads = pd.NA

                fused_row = dict(base)
                fused_row.update({
                    "chrom": str(chrom),
                    "start": int(new_start),
                    "end": int(new_end),
                    "name": str(block.iloc[0].get("name", "")),

                    "y": float(fused_y),
                    "ratio": float(fused_ratio),
                    "cn": float(fused_cn),
                    "n_probes": int(fused_nprobes) if math.isfinite(fused_nprobes) else pd.NA,
                    "direction": direction,

                    "left_boundary_name": left_name,
                    "right_boundary_name": right_name,
                    "left_boundary_z": left_z_val,
                    "right_boundary_z": right_z_val,
                    "seg_conf_z": float(fused_seg_conf_z),

                    "ratio_tier": fused_rtier,
                    "support_tier": fused_ztier,
                    "keep_suggested": bool(fused_keep),
                    "call": fused_call,

                    "ratio_z": float(fused_ratio_z) if math.isfinite(float(fused_ratio_z)) else float("nan"),


                    "seg_len_bp": fused_seg_len_bp,
                    "seg_primary_reads": fused_seg_primary_reads,

                    "left_split_reads": fused_left_split_reads,
                    "left_disco_reads": fused_left_disco_reads,
                    "left_support_reads": fused_left_support_reads,

                    "right_split_reads": fused_right_split_reads,
                    "right_disco_reads": fused_right_disco_reads,
                    "right_support_reads": fused_right_support_reads,

                    "fused": True,
                    "fused_n_segments": int(len(names_uniq)),
                    "fused_internal_boundaries": int(max(0, len(names_uniq) - 1)),
                    "fused_segment_names": fused_names_str,

                    "internal_left_boundary_support_mean": float(fused_internal_left_support_mean) if math.isfinite(float(fused_internal_left_support_mean)) else float("nan"),
                    "internal_right_boundary_support_mean": float(fused_internal_right_support_mean) if math.isfinite(float(fused_internal_right_support_mean)) else float("nan"),
                    "fused_left_boundaries_n": int(fused_left_boundaries_n),
                    "fused_right_boundaries_n": int(fused_right_boundaries_n),

                    "outer_left_split_reads": fused_left_split_reads,
                    "outer_left_disco_reads": fused_left_disco_reads,
                    "outer_left_support_reads": fused_left_support_reads,

                    "outer_right_split_reads": fused_right_split_reads,
                    "outer_right_disco_reads": fused_right_disco_reads,
                    "outer_right_support_reads": fused_right_support_reads,
                })

                fused_row.pop("_orig_idx", None)
                fused_row.pop("_eligible", None)
                fused_rows.append(fused_row)

            i = j

    if not fused_rows:
        out2 = out.copy()
        out2["fused"] = False
        out2 = _ensure_fusion_schema(out2)
        return out2

    fused_df = pd.DataFrame(fused_rows)

    dedupe_keys = [c for c in ["chrom", "start", "end", "direction", "fused_segment_names"] if c in fused_df.columns]
    if dedupe_keys:
        fused_df = fused_df.drop_duplicates(subset=dedupe_keys, keep="last").reset_index(drop=True)

    out2 = out.copy()
    if "fused" not in out2.columns:
        out2["fused"] = False

    out2["_fused_member"] = out2.index.map(lambda idx: idx in fused_member_idx)
    originals_kept = out2.loc[~out2["_fused_member"]].drop(columns=["_fused_member"])

    combined = pd.concat([originals_kept, fused_df], ignore_index=True, sort=False)

    combined["fused"] = combined.get("fused", False)
    combined["fused"] = combined["fused"].fillna(False).astype(bool)

    combined["start"] = pd.to_numeric(combined["start"], errors="coerce").astype(int)
    combined["end"] = pd.to_numeric(combined["end"], errors="coerce").astype(int)

    combined["fused_sort"] = combined["fused"].astype(int)
    combined = combined.sort_values(["chrom", "start", "end", "fused_sort"], kind="mergesort") \
                       .drop(columns=["fused_sort"])

    combined = _ensure_fusion_schema(combined)
    return combined


# ============================================================
# Main
# ============================================================

def main():
    ap = argparse.ArgumentParser(
        description="Finalize CNV segments by combining segment ratio (effect size tiers) with boundary z support (confidence tiers), "
                    "with optional BAM-based read support and optional fusion."
    )
    ap.add_argument("--segments-bed", required=True, help="segments.bed(.gz ok) from Step 2")
    ap.add_argument("--boundary-support-tsv", required=True, help="boundary_support.tsv(.gz ok) from Step 3")
    ap.add_argument("--out-tsv", required=True, help="Output TSV (plain .tsv recommended)")

    ap.add_argument("--out-bedgraph-y", default="", help="Optional bedGraph of y per segment (.gz ok)")
    ap.add_argument("--out-bedgraph-ratio", default="", help="Optional bedGraph of ratio per segment (.gz ok)")

    ap.add_argument("--weak-ratio", type=float, default=1.25)
    ap.add_argument("--strong-ratio", type=float, default=1.5)

    ap.add_argument("--weak-z", type=float, default=1.0)
    ap.add_argument("--strong-z", type=float, default=2.0)
    ap.add_argument(
        "--confidence-method",
        choices=["min", "mean", "max"],
        default="mean",
        help="How to combine left/right boundary z into seg_conf_z for per-segment calls"
    )
    ap.add_argument("--keep-mode", choices=["ratio_only", "z_only", "both"], default="ratio_only")
    ap.add_argument("--cn-base", type=float, default=1.0)

    ap.add_argument("--fuse", action="store_true", help="Fuse adjacent kept CNVs of same direction (allow weak+strong)")
    ap.add_argument("--fuse-max-gap", type=int, default=0, help="Maximum bp gap allowed between fused segments")

    ap.add_argument("--primary-bam", default="", help="Primary/aligned BAM for segment read counting (.bai required)")
    ap.add_argument("--split-bam", default="", help="Splitters BAM for boundary read counting (.bai required)")
    ap.add_argument("--disco-bam", default="", help="Discordant BAM for boundary read counting (.bai required)")
    ap.add_argument("--count-flank", type=int, default=1250, help="Flank bp on each side of boundary for read counting")

    ap.add_argument("--count-all-alignments", action="store_true", default=False,
                    help="Count all alignments (otherwise counts read1 only for paired reads)")
    ap.add_argument("--count-drop-dup", action="store_true", default=False,
                    help="Exclude duplicate-marked reads (FLAG 1024) from counts")
    ap.add_argument("--count-include-supp", action="store_true", default=False,
                    help="Include supplementary alignments in counts")
    ap.add_argument("--count-include-secondary", action="store_true", default=False,
                    help="Include secondary alignments in counts")

    ap.add_argument(
        "--circular-contigs",
        default="",
        help="Comma-separated contig names to treat as circular (wrap support at terminal boundaries).",
    )
    ap.add_argument(
        "--circular-contigs-file",
        default="",
        help="File with one circular contig name per line (# comments ok).",
    )

    ap.add_argument("--debug", action="store_true", default=False,
                    help="Verbose debug prints (recommended while validating circular wrap counting)")

    args = ap.parse_args()

    if args.weak_ratio <= 1.0 or args.strong_ratio <= 1.0:
        raise RuntimeError("--weak-ratio and --strong-ratio must be > 1.0")
    if args.strong_ratio < args.weak_ratio:
        raise RuntimeError("--strong-ratio must be >= --weak-ratio")
    if args.weak_z < 0 or args.strong_z < 0:
        raise RuntimeError("--weak-z/--strong-z must be >= 0")
    if args.strong_z < args.weak_z:
        raise RuntimeError("--strong-z must be >= --weak-z")

    seg = read_segments_bed(args.segments_bed)
    b = read_boundary_support(args.boundary_support_tsv)

    bz = (
        b.groupby(["chrom", "name"], as_index=False)["combined_z_chr"]
        .max()
        .rename(columns={"combined_z_chr": "boundary_z"})
    )
    zmap = {(r.chrom, r.name): float(r.boundary_z) for r in bz.itertuples(index=False)}

    primary_bam = _open_bam(args.primary_bam) if args.primary_bam else None
    split_bam = _open_bam(args.split_bam) if args.split_bam else None
    disco_bam = _open_bam(args.disco_bam) if args.disco_bam else None

    circular_contigs = parse_circular_contigs(args.circular_contigs, args.circular_contigs_file)

    contig_lengths: Dict[str, int] = {}
    for bam in (primary_bam, split_bam, disco_bam):
        contig_lengths = build_contig_length_map(bam)
        if contig_lengths:
            break

    print(f"[cnv_finalize_segments] circular args csv='{args.circular_contigs}' file='{args.circular_contigs_file}'")
    print(f"[cnv_finalize_segments] parsed circular contigs (raw+norm): {sorted(list(circular_contigs))}")

    seg_contigs = set(seg["chrom"].astype(str).tolist())
    seg_contigs_norm = set(_norm_contig(c) for c in seg_contigs)

    bam_contigs = set()
    bam_contigs_norm = set()
    for bam in (primary_bam, split_bam, disco_bam):
        if bam is not None:
            bam_contigs = set(bam.references)
            bam_contigs_norm = set(_norm_contig(c) for c in bam.references)
            break

    if not circular_contigs:
        print("[cnv_finalize_segments] NOTE: no circular contigs provided; wrap logic disabled.")
    else:
        matched = []
        for c in sorted(circular_contigs):
            if (c in seg_contigs) or (_norm_contig(c) in seg_contigs_norm) or (c in bam_contigs) or (_norm_contig(c) in bam_contigs_norm):
                matched.append(c)
        if not matched:
            print("[cnv_finalize_segments] WARN: none of --circular-contigs matched contigs in segments/BAM.")
        else:
            print(f"[cnv_finalize_segments] circular contigs matched (at least one view): {sorted(matched)}")

    print(f"[cnv_finalize_segments] example seg contigs: {sorted(list(seg_contigs))[:10]}")
    if bam_contigs:
        print(f"[cnv_finalize_segments] example BAM contigs: {sorted(list(bam_contigs))[:10]}")
    if args.debug:
        print(f"[cnv_finalize_segments][DEBUG] contig_lengths keys sample: {sorted(list(contig_lengths.keys()))[:15]}")

    read1_only = not bool(args.count_all_alignments)
    drop_dup = bool(args.count_drop_dup)
    include_supp = bool(args.count_include_supp)
    include_secondary = bool(args.count_include_secondary)
    flank = int(args.count_flank)

    out_rows = []

    for chrom, sub in seg.groupby("chrom", sort=False):
        sub = sub.sort_values(["start", "end"]).reset_index(drop=True)

        chrom_is_circ = is_circular_contig(chrom, circular_contigs)

        if args.debug and chrom_is_circ:
            split_res = resolve_bam_contig_name(split_bam, chrom) if split_bam is not None else None
            disco_res = resolve_bam_contig_name(disco_bam, chrom) if disco_bam is not None else None
            primary_res = resolve_bam_contig_name(primary_bam, chrom) if primary_bam is not None else None
            Ls = get_contig_len(split_res or chrom, contig_lengths) if split_res is not None else get_contig_len(chrom, contig_lengths)
            print(f"[cnv_finalize_segments][DEBUG] chrom={chrom} CIRCULAR=True n_segments={len(sub)} "
                  f"split_res={split_res} disco_res={disco_res} primary_res={primary_res} L={Ls}")

        for i in range(len(sub)):
            r = {k: sub.loc[i, k] for k in sub.columns}
            name = str(r["name"])

            is_first = (i == 0)
            is_last = (i == (len(sub) - 1))

            if is_first:
                left_name = "NA"
                left_z = float("nan")
                left_anchor_pos_0 = int(r["start"])
                left_boundary_pos_0 = None
            else:
                left_name = str(sub.loc[i - 1, "name"])
                left_z = zmap.get((chrom, left_name), float("nan"))
                left_anchor_pos_0 = None
                left_boundary_pos_0 = int(r["start"])

            if is_last:
                right_name = "NA"
                right_z = float("nan")
                right_anchor_pos_0 = int(r["end"])
                right_boundary_pos_0 = None
            else:
                right_name = name
                right_z = zmap.get((chrom, right_name), float("nan"))
                right_anchor_pos_0 = None
                right_boundary_pos_0 = int(r["end"])

            if args.debug and chrom_is_circ and (is_first or is_last):
                print(f"[cnv_finalize_segments][DEBUG] chrom={chrom} seg={name} idx={i}/{len(sub)-1} "
                      f"FORCE_ENDS circ=True is_first={is_first} is_last={is_last} "
                      f"Lpos={left_boundary_pos_0} Rpos={right_boundary_pos_0} "
                      f"Lanchor={left_anchor_pos_0} Ranchor={right_anchor_pos_0}")

            z_pair = []
            for z in [left_z, right_z]:
                try:
                    zf = float(z)
                    if math.isfinite(zf):
                        z_pair.append(zf)
                except Exception:
                    pass

            conf_method = str(args.confidence_method).lower()

            if len(z_pair) == 0:
                seg_conf_z = float("nan")
            elif conf_method == "mean":
                seg_conf_z = float(np.mean(z_pair))
            elif conf_method == "max":
                seg_conf_z = float(np.max(z_pair))
            else:
                seg_conf_z = float(np.min(z_pair))


            ratio = float(r["ratio"])
            y = float(r["y"])
            cn = float(args.cn_base) * ratio

            rtier = ratio_tier(ratio, args.weak_ratio, args.strong_ratio)
            ztier = z_tier(seg_conf_z, args.weak_z, args.strong_z)
            direction = call_direction(ratio)

            ratio_keep = (rtier in ["weak", "strong"])
            z_keep = (ztier in ["weak", "strong"])

            if args.keep_mode == "ratio_only":
                keep_suggested = ratio_keep
            elif args.keep_mode == "z_only":
                keep_suggested = z_keep
            else:
                keep_suggested = ratio_keep and z_keep

            if not ratio_keep:
                call = "no_CNV_effect"
            else:
                call = f"CNV_{direction}_{rtier}Effect_{ztier}Support"

            seg_len = int(r["end"]) - int(r["start"])
            seg_primary_reads = pd.NA

            left_split_reads = pd.NA
            right_split_reads = pd.NA
            left_disco_reads = pd.NA
            right_disco_reads = pd.NA

            if primary_bam is not None:
                seg_primary_reads = count_reads_in_interval(
                    primary_bam, chrom, int(r["start"]), int(r["end"]),
                    read1_only=read1_only, drop_dup=drop_dup,
                    include_supp=include_supp, include_secondary=include_secondary
                )

            if split_bam is not None:
                ls = count_unique_qnames_for_boundary(
                    split_bam, chrom,
                    boundary_pos_0based=left_boundary_pos_0,
                    anchor_pos_0based=left_anchor_pos_0,
                    flank=flank,
                    circular_contigs=circular_contigs,
                    contig_lengths=contig_lengths,
                    read1_only=read1_only, drop_dup=drop_dup,
                    include_supp=False, include_secondary=include_secondary,
                    debug=bool(args.debug and chrom_is_circ and is_first),
                    debug_tag=f" split LEFT seg={name}"
                )
                rs = count_unique_qnames_for_boundary(
                    split_bam, chrom,
                    boundary_pos_0based=right_boundary_pos_0,
                    anchor_pos_0based=right_anchor_pos_0,
                    flank=flank,
                    circular_contigs=circular_contigs,
                    contig_lengths=contig_lengths,
                    read1_only=read1_only, drop_dup=drop_dup,
                    include_supp=False, include_secondary=include_secondary,
                    debug=bool(args.debug and chrom_is_circ and is_last),
                    debug_tag=f" split RIGHT seg={name}"
                )
                left_split_reads = pd.NA if ls is None else int(ls)
                right_split_reads = pd.NA if rs is None else int(rs)

            if disco_bam is not None:
                ld = count_unique_qnames_for_boundary(
                    disco_bam, chrom,
                    boundary_pos_0based=left_boundary_pos_0,
                    anchor_pos_0based=left_anchor_pos_0,
                    flank=flank,
                    circular_contigs=circular_contigs,
                    contig_lengths=contig_lengths,
                    read1_only=read1_only, drop_dup=drop_dup,
                    include_supp=False, include_secondary=include_secondary,
                    debug=bool(args.debug and chrom_is_circ and is_first),
                    debug_tag=f" disco LEFT seg={name}"
                )
                rd = count_unique_qnames_for_boundary(
                    disco_bam, chrom,
                    boundary_pos_0based=right_boundary_pos_0,
                    anchor_pos_0based=right_anchor_pos_0,
                    flank=flank,
                    circular_contigs=circular_contigs,
                    contig_lengths=contig_lengths,
                    read1_only=read1_only, drop_dup=drop_dup,
                    include_supp=False, include_secondary=include_secondary,
                    debug=bool(args.debug and chrom_is_circ and is_last),
                    debug_tag=f" disco RIGHT seg={name}"
                )
                left_disco_reads = pd.NA if ld is None else int(ld)
                right_disco_reads = pd.NA if rd is None else int(rd)

            left_support_reads = pd.NA
            right_support_reads = pd.NA
            if pd.notna(left_split_reads) and pd.notna(left_disco_reads):
                left_support_reads = int(left_split_reads) + int(left_disco_reads)
            if pd.notna(right_split_reads) and pd.notna(right_disco_reads):
                right_support_reads = int(right_split_reads) + int(right_disco_reads)

            r.update(dict(
                cn_base=float(args.cn_base),
                cn=cn,
                direction=direction,

                left_boundary_name=left_name,
                right_boundary_name=right_name,
                left_boundary_z=left_z,
                right_boundary_z=right_z,
                seg_conf_z=seg_conf_z,

                ratio_tier=rtier,
                support_tier=ztier,
                weak_ratio=args.weak_ratio,
                strong_ratio=args.strong_ratio,
                weak_z=args.weak_z,
                strong_z=args.strong_z,

                keep_mode=args.keep_mode,
                keep_suggested=bool(keep_suggested),
                call=call,

                seg_len_bp=int(seg_len),
                seg_primary_reads=seg_primary_reads,

                left_split_reads=left_split_reads,
                right_split_reads=right_split_reads,
                left_disco_reads=left_disco_reads,
                right_disco_reads=right_disco_reads,
                left_support_reads=left_support_reads,
                right_support_reads=right_support_reads,
            ))

            out_rows.append(r)

    out = pd.DataFrame(out_rows)

    y_abs = pd.to_numeric(out["y"], errors="coerce").abs()

    def _robust_mag_z(g: pd.Series) -> pd.Series:
        v = pd.to_numeric(g, errors="coerce").to_numpy(dtype=float)
        med = np.nanmedian(v)
        mad = np.nanmedian(np.abs(v - med))
        denom = mad if (np.isfinite(mad) and mad > 0) else np.nanstd(v)
        if not np.isfinite(denom) or denom <= 0:
            return pd.Series([0.0] * len(v), index=g.index, dtype=float)
        z = 0.6745 * (v - med) / denom
        z = np.abs(z)
        z = np.where(np.isfinite(z), z, 0.0)
        return pd.Series(z, index=g.index, dtype=float)

    out["ratio_z"] = y_abs.groupby(out["chrom"], sort=False).apply(_robust_mag_z).reset_index(level=0, drop=True)

    if args.fuse:
        out = fuse_segments_preserve_support(out, max_gap=int(args.fuse_max_gap))
        out = _ensure_fusion_schema(out)
    else:
        out = _ensure_fusion_schema(out)

    col_order = [
        "chrom", "start", "end", "name",
        "y", "ratio", "n_probes",
        "cn_base", "cn", "direction",
        "left_boundary_name", "right_boundary_name",
        "left_boundary_z", "right_boundary_z", "seg_conf_z",
        "ratio_tier", "support_tier",
        "weak_ratio", "strong_ratio", "weak_z", "strong_z",
        "keep_mode", "keep_suggested", "call",

        "seg_len_bp", "seg_primary_reads",
        "left_split_reads", "left_disco_reads", "left_support_reads",
        "right_split_reads", "right_disco_reads", "right_support_reads",

        "fused",
        "fused_n_segments",
        "fused_internal_boundaries",
        "fused_segment_names",
        "internal_left_boundary_support_mean",
        "internal_right_boundary_support_mean",
        "fused_left_boundaries_n",
        "fused_right_boundaries_n",
        "outer_left_split_reads", "outer_left_disco_reads", "outer_left_support_reads",
        "outer_right_split_reads", "outer_right_disco_reads", "outer_right_support_reads",

        "ratio_z",
    ]
    out = out[[c for c in col_order if c in out.columns] + [c for c in out.columns if c not in col_order]]

    _atomic_write_no_follow(write_tsv_any, out, args.out_tsv)

    if args.fuse:
        out_dir = os.path.dirname(os.path.abspath(args.out_tsv)) or "."
        base = os.path.basename(args.out_tsv)

        if base.endswith(".tsv.gz"):
            fused_only_name = base[:-7] + ".fused.tsv.gz"
        elif base.endswith(".tsv"):
            fused_only_name = base[:-4] + ".fused.tsv"
        elif base.endswith(".gz"):
            fused_only_name = base[:-3] + ".fused.tsv.gz"
        else:
            fused_only_name = base + ".fused.tsv"

        fused_only_path = os.path.join(out_dir, fused_only_name)

        out_fused_only = out[out["fused"].astype(bool)].copy()
        if out_fused_only.empty:
            out_fused_only = out.head(0).copy()
        _atomic_write_no_follow(write_tsv_any, out_fused_only, fused_only_path)

    if args.out_bedgraph_y:
        bg = out[["chrom", "start", "end"]].copy()
        bg["value"] = pd.to_numeric(out["y"], errors="coerce").astype(float)
        _atomic_write_no_follow(write_bedgraph_any, bg, args.out_bedgraph_y)

    if args.out_bedgraph_ratio:
        bg = out[["chrom", "start", "end"]].copy()
        bg["value"] = pd.to_numeric(out["ratio"], errors="coerce").astype(float)
        _atomic_write_no_follow(write_bedgraph_any, bg, args.out_bedgraph_ratio)

    for bam in (primary_bam, split_bam, disco_bam):
        try:
            if bam is not None:
                bam.close()
        except Exception:
            pass


if __name__ == "__main__":
    main()
