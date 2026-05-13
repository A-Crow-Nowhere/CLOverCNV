#!/usr/bin/env python3
import argparse
import gzip
import math
import os
import subprocess
import sys
from pathlib import Path
from typing import Dict, Optional, Tuple

import numpy as np
import pandas as pd
import pyBigWig


def _resolve_genome_fasta(genome_arg: str) -> Path:
    """
    Resolve genome FASTA from:
      - explicit FASTA path
      - directory containing FASTA
      - basename + common FASTA suffixes in CWD
      - MAPI genome keys (e.g. "Dd2") via:
          $CLOVER_GENOMES_DIR, else $MAPI_GENOMES_DIR, else $MAPI_ROOT/genomes
    """
    p = Path(genome_arg)

    # 1) Exact file path
    if p.is_file():
        return p

    # 2) Directory: pick genome.fa|genome.fasta|... or first matching *.fa*
    if p.is_dir():
        for nm in [
            "genome.fa", "genome.fasta", "genome.fna",
            "genome.fa.gz", "genome.fasta.gz", "genome.fna.gz"
        ]:
            q = p / nm
            if q.is_file():
                return q
        patterns = ["*.fa", "*.fasta", "*.fna", "*.fa.gz", "*.fasta.gz", "*.fna.gz"]
        cand = []
        for pat in patterns:
            cand.extend(sorted(p.glob(pat)))
        if cand:
            return cand[0]
        raise RuntimeError(f"Genome is a directory but no FASTA found in: {p}")

    # 3) Basename + suffix in current directory
    for suf in [".fa", ".fasta", ".fna", ".fa.gz", ".fasta.gz", ".fna.gz"]:
        q = Path(str(p) + suf)
        if q.is_file():
            return q

    # 4) MAPI key resolution: look under genomes directory conventions
    #    Priority: explicit env vars > MAPI_ROOT/genomes
    genomes_dir = (
        os.environ.get("CLOVER_GENOMES_DIR")
        or os.environ.get("MAPI_GENOMES_DIR")
        or (Path(os.environ["MAPI_ROOT"]) / "genomes" if "MAPI_ROOT" in os.environ else None)
    )

    if genomes_dir:
        gd = Path(genomes_dir)
        # allow genome_arg like "Dd2" or "Dd2.fasta" etc.
        # try direct join first, then suffixes
        direct = gd / genome_arg
        if direct.is_file():
            return direct
        for suf in [".fa", ".fasta", ".fna", ".fa.gz", ".fasta.gz", ".fna.gz"]:
            q = gd / f"{genome_arg}{suf}"
            if q.is_file():
                return q

    raise RuntimeError(
        f"Could not resolve genome FASTA from: {genome_arg}. "
        f"Provide a FASTA path, a directory containing a FASTA, "
        f"or set MAPI_ROOT (or CLOVER_GENOMES_DIR/MAPI_GENOMES_DIR) for genome-key resolution."
    )

def load_chr_sizes_from_genome(genome_arg: str) -> Dict[str, int]:
    genome = _resolve_genome_fasta(genome_arg)
    fai = Path(str(genome) + ".fai")
    if not fai.exists():
        try:
            print(f"[boundary_support] .fai not found; attempting: samtools faidx {genome}", file=sys.stderr)
            subprocess.run(["samtools", "faidx", str(genome)], check=True)
        except Exception as e:
            raise RuntimeError(f"Missing FASTA index {fai} and failed to create it via samtools faidx.") from e
    if not fai.exists():
        raise RuntimeError(f"FASTA index still missing: {fai}")

    sizes: Dict[str, int] = {}
    with open(fai, "r") as f:
        for line in f:
            if not line.strip():
                continue
            chrom, length = line.rstrip("\n").split("\t")[:2]
            sizes[chrom] = int(length)
    if not sizes:
        raise RuntimeError(f"No chromosome sizes parsed from {fai}")
    print(f"[boundary_support] loaded {len(sizes)} chromosomes from {fai}", file=sys.stderr)
    return sizes


def load_chr_sizes_from_file(path: str) -> Dict[str, int]:
    sizes: Dict[str, int] = {}
    with open(path, "r") as f:
        for line in f:
            if not line.strip():
                continue
            chrom, length = line.rstrip("\n").split("\t")[:2]
            sizes[chrom] = int(length)
    if not sizes:
        raise RuntimeError(f"No chromosome sizes parsed from {path}")
    return sizes


def open_bw_or_die(path: str, label: str):
    p = Path(path)
    if not p.exists():
        raise RuntimeError(f"{label} bigWig not found: {path}")
    if p.is_dir():
        raise RuntimeError(f"{label} bigWig path is a directory: {path}")
    if p.stat().st_size == 0:
        raise RuntimeError(f"{label} bigWig is empty (0 bytes): {path}")
    try:
        bw = pyBigWig.open(str(p))
    except RuntimeError as e:
        raise RuntimeError(f"Failed to open {label} bigWig: {path}") from e
    if bw is None:
        raise RuntimeError(f"pyBigWig returned NULL for {label}: {path}")
    chroms = bw.chroms()
    if not chroms:
        bw.close()
        raise RuntimeError(f"{label} bigWig has no chromosomes: {path}")
    return bw


def bw_stat(bw, chrom: str, start: int, end: int, stat: str) -> float:
    """
    stat: 'mean' or 'sum' (sum approximated as mean*len if sum unsupported)
    Returns 0.0 on missing/no-data.
    """
    if bw is None or end <= start:
        return 0.0

    try:
        if stat == "mean":
            v = bw.stats(chrom, start, end, type="mean")
            x = v[0] if v and v[0] is not None else 0.0
            return 0.0 if (x is None or (isinstance(x, float) and math.isnan(x))) else float(x)

        if stat == "sum":
            # pyBigWig supports type='sum' in many builds; if not, fallback to mean*len
            try:
                v = bw.stats(chrom, start, end, type="sum")
                x = v[0] if v and v[0] is not None else None
                if x is None or (isinstance(x, float) and math.isnan(x)):
                    raise RuntimeError("sum returned None/NaN")
                return float(x)
            except Exception:
                v = bw.stats(chrom, start, end, type="mean")
                x = v[0] if v and v[0] is not None else 0.0
                x = 0.0 if (x is None or (isinstance(x, float) and math.isnan(x))) else float(x)
                return x * float(end - start)

        raise ValueError("stat must be mean or sum")
    except RuntimeError:
        return 0.0


def robust_z(x: np.ndarray, eps: float = 1e-9) -> np.ndarray:
    """
    Robust z-score using median and MAD; returns zeros if MAD is ~0.
    """
    if x.size == 0:
        return x
    med = np.nanmedian(x)
    mad = np.nanmedian(np.abs(x - med))
    if not np.isfinite(mad) or mad < eps:
        return np.zeros_like(x, dtype=float)
    return 0.6745 * (x - med) / mad


def main():
    ap = argparse.ArgumentParser(
        description="Compute boundary support from split/disco bigWigs in flanking windows (RAW pileup; no normalization to primary depth)."
    )
    ap.add_argument("--boundaries-bed", required=True, help="Boundaries BED(.gz ok), e.g. <out-prefix>/boundaries.bed.gz from Step 2.")
    ap.add_argument("--genome", default="", help="Genome FASTA (or dir/prefix) to derive chromosome sizes via .fai (recommended).")
    ap.add_argument("--chr-sizes", default="", help="Optional chrom.sizes file (chrom<TAB>len). Overrides --genome if provided.")
    ap.add_argument("--split-bw", required=True, help="Split-read bigWig (raw signal).")
    ap.add_argument("--disco-bw", required=True, help="Discordant-read bigWig (raw signal).")
    ap.add_argument("--flank", type=int, default=1250, help="Flank size (bp). Default 1250.")
    ap.add_argument("--min-flank", type=int, default=50, help="Minimum effective flank length. Default 50.")
    ap.add_argument("--agg", choices=["sum", "mean"], default="sum",
                    help="Aggregation for support. Default sum (raw pileup in window).")
    ap.add_argument("--w-split", type=float, default=1.0, help="Weight for split support.")
    ap.add_argument("--w-disco", type=float, default=1.0, help="Weight for disco support.")
    ap.add_argument("--eps", type=float, default=1e-9, help="Small epsilon to stabilize math.")
    ap.add_argument("--out-tsv", required=True, help="Output TSV(.gz ok).")
    ap.add_argument("--out-bedgraph", default="", help="Optional 1bp bedGraph of combined_z_chr or combined_support.")
    ap.add_argument("--bedgraph-field", choices=["combined_support", "combined_z_chr"], default="combined_z_chr",
                    help="What to emit in bedGraph. Default combined_z_chr.")

    args = ap.parse_args()

    # chrom sizes
    if args.chr_sizes:
        chr_sizes = load_chr_sizes_from_file(args.chr_sizes)
    elif args.genome:
        chr_sizes = load_chr_sizes_from_genome(args.genome)
    else:
        raise RuntimeError("Provide either --genome (recommended) or --chr-sizes.")

    # open bigwigs
    bw_split = open_bw_or_die(args.split_bw, "split")
    bw_disco = open_bw_or_die(args.disco_bw, "disco")

    # read boundaries BED
    # expected: chrom start end name score strand
    if args.boundaries_bed.endswith(".gz"):
        with gzip.open(args.boundaries_bed, "rt") as f:
            bed = pd.read_csv(f, sep="\t", header=None)
    else:
        bed = pd.read_csv(args.boundaries_bed, sep="\t", header=None)

    if bed.shape[1] < 4:
        raise RuntimeError("boundaries BED must have at least 4 columns: chrom start end name")
    bed = bed.rename(columns={0: "chrom", 1: "start", 2: "end", 3: "name"})
    bed["start"] = bed["start"].astype(int)
    bed["end"] = bed["end"].astype(int)

    rows = []
    for _, r in bed.iterrows():
        chrom = str(r["chrom"])
        name = str(r["name"])

        if chrom not in chr_sizes:
            continue

        # boundary position: use BED end (1bp interval); fallback to midpoint
        pos = int(r["end"]) if int(r["end"]) > int(r["start"]) else int((int(r["start"]) + int(r["end"])) // 2)
        clen = int(chr_sizes[chrom])

        Ls = max(0, pos - int(args.flank))
        Le = pos
        Rs = pos
        Re = min(clen, pos + int(args.flank))

        # enforce minimum flank length
        if (Le - Ls) < args.min_flank:
            Ls = max(0, Le - args.min_flank)
        if (Re - Rs) < args.min_flank:
            Re = min(clen, Rs + args.min_flank)

        # RAW support in windows: default SUM = "pileup in window"
        split_left = bw_stat(bw_split, chrom, Ls, Le, args.agg)
        split_right = bw_stat(bw_split, chrom, Rs, Re, args.agg)
        disco_left = bw_stat(bw_disco, chrom, Ls, Le, args.agg)
        disco_right = bw_stat(bw_disco, chrom, Rs, Re, args.agg)

        split_support = max(split_left, split_right)
        disco_support = max(disco_left, disco_right)

        combined_support = (args.w_split * split_support) + (args.w_disco * disco_support)

        rows.append(
            dict(
                chrom=chrom,
                pos=pos,
                name=name,
                flank=args.flank,
                agg=args.agg,

                split_left=split_left,
                split_right=split_right,
                split_support=split_support,

                disco_left=disco_left,
                disco_right=disco_right,
                disco_support=disco_support,

                combined_support=combined_support,

                left_start=Ls,
                left_end=Le,
                right_start=Rs,
                right_end=Re,
            )
        )

    out = pd.DataFrame(rows)
    if out.empty:
        raise RuntimeError("No boundaries processed (chrom mismatch or empty BED).")

    # robust z per chromosome
    out["split_z_chr"] = 0.0
    out["disco_z_chr"] = 0.0
    out["combined_z_chr"] = 0.0
    out["support_tier_chr"] = "D"

    for chrom, idx in out.groupby("chrom").groups.items():
        i = np.array(list(idx), dtype=int)

        sz = robust_z(out.loc[i, "split_support"].to_numpy(dtype=float), eps=args.eps)
        dz = robust_z(out.loc[i, "disco_support"].to_numpy(dtype=float), eps=args.eps)
        cz = robust_z(out.loc[i, "combined_support"].to_numpy(dtype=float), eps=args.eps)

        out.loc[i, "split_z_chr"] = sz
        out.loc[i, "disco_z_chr"] = dz
        out.loc[i, "combined_z_chr"] = cz

        # simple tiers by combined_z_chr
        tier = np.full(i.shape, "D", dtype=object)
        tier[cz >= 0.5] = "C"
        tier[cz >= 1.0] = "B"
        tier[cz >= 2.0] = "A"
        out.loc[i, "support_tier_chr"] = tier

    # write TSV
    if args.out_tsv.endswith(".gz"):
        with gzip.open(args.out_tsv, "wt") as f:
            out.to_csv(f, sep="\t", index=False)
    else:
        out.to_csv(args.out_tsv, sep="\t", index=False)

    # optional bedGraph (1bp bars)
    if args.out_bedgraph:
        bg = out[["chrom"]].copy()
        bg["start"] = (out["pos"].astype(int) - 1).clip(lower=0)
        bg["end"] = out["pos"].astype(int)
        bg["value"] = out[args.bedgraph_field].astype(float)

        if args.out_bedgraph.endswith(".gz"):
            with gzip.open(args.out_bedgraph, "wt") as f:
                bg.to_csv(f, sep="\t", header=False, index=False)
        else:
            bg.to_csv(args.out_bedgraph, sep="\t", header=False, index=False)

    bw_split.close()
    bw_disco.close()

    print(f"[boundary_support] wrote: {args.out_tsv}", file=sys.stderr)
    if args.out_bedgraph:
        print(f"[boundary_support] wrote: {args.out_bedgraph}", file=sys.stderr)


if __name__ == "__main__":
    main()
