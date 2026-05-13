#!/usr/bin/env bash
set -euo pipefail

die() { echo "ERROR: $*" >&2; exit 1; }
say() { echo "[mapi][cov_model_apply] $*" >&2; }

usage() {
  cat <<'EOF'
Usage:
  mapi modules cov_model_apply \
    --dir <bam_dir> \
    --model <model.rds> \
    --genome <GENOME_KEY|/path/to.fa[.gz]> \
    --out-dir <DIR> \
    [--sample-key <SAMPLE>] \
    [--binsize <bp>] \
    [--threads <n>] \
    [--primary-prefer primary|aligned] \
    [--mappability-track <tsv>] \
    [--emit-raw] \
    [--r-modules "<module list>"] \
    [--r-exe <path/to/Rscript>] \
    [--keep-tmp]

What it does:
  - For each sample in --dir (or one sample via --sample-key):
    1) Select a "primary" BAM (prefers *.primary* or *.aligned* depending on --primary-prefer)
    2) Build fixed bins over genome (default 50bp)
    3) Compute AT per bin (from FASTA)
    4) Compute mean depth per bin from PRIMARY BAM (samtools depth; streaming)
    5) Apply model to predict expected depth per bin
    6) Compute correction factor: cfactor = median(expected)/expected
    7) For EVERY BAM belonging to that sample:
         corr_depth_bw (raw * cfactor) [always]
         raw_depth_bw  (raw)           [only if --emit-raw]
    8) Write factor bigWigs (expected + cfactor)

Outputs (BigWig only):
  <out-dir>/<sample>/corr_bw/*.bw
  <out-dir>/<sample>/raw_bw/*.bw        (only if --emit-raw)
  <out-dir>/<sample>/factors/*.bw       (expected + cfactor)

Notes:
  - Mappability is OFF by default. If enabled, user must provide --mappability-track.
  - Depth is mean pileup depth per bin (samtools depth). This avoids bedtools OOM.
EOF
}

# ---------------------------
# Defaults
# ---------------------------
DIR=""
MODEL_RDS=""
OUT_DIR=""
SAMPLE_KEY=""
GENOME=""
BINSIZE=50
THREADS=4

PRIMARY_PREFER="primary"   # primary|aligned

MAP_TRACK=""
NO_MAPPABILITY=1           # OFF by default

EMIT_RAW=0                 # OFF by default; opt-in via --emit-raw

R_MODULES=""
R_EXE="Rscript"

KEEP_TMP=0

# ---------------------------
# Parse args
# ---------------------------
if [[ $# -eq 0 ]]; then usage; exit 0; fi
while [[ $# -gt 0 ]]; do
  case "$1" in
    --dir) DIR="${2:-}"; shift 2 ;;
    --model) MODEL_RDS="${2:-}"; shift 2 ;;
    --out-dir) OUT_DIR="${2:-}"; shift 2 ;;
    --sample-key) SAMPLE_KEY="${2:-}"; shift 2 ;;
    --genome) GENOME="${2:-}"; shift 2 ;;
    --binsize) BINSIZE="${2:-}"; shift 2 ;;
    --threads) THREADS="${2:-}"; shift 2 ;;
    --primary-prefer) PRIMARY_PREFER="${2:-}"; shift 2 ;;
    --mappability-track) MAP_TRACK="${2:-}"; NO_MAPPABILITY=0; shift 2 ;;
    --emit-raw) EMIT_RAW=1; shift 1 ;;
    --r-modules) R_MODULES="${2:-}"; shift 2 ;;
    --r-exe) R_EXE="${2:-}"; shift 2 ;;
    --keep-tmp) KEEP_TMP=1; shift 1 ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown option: $1" ;;
  esac
done

# ---------------------------
# Validate
# ---------------------------
[[ -n "$DIR" ]] || die "--dir is required"
[[ -d "$DIR" ]] || die "Directory not found: $DIR"
[[ -n "$MODEL_RDS" ]] || die "--model is required"
[[ -f "$MODEL_RDS" ]] || die "Model not found: $MODEL_RDS"
[[ -n "$OUT_DIR" ]] || die "--out-dir is required"
[[ -n "$GENOME" ]] || die "--genome is required"
mkdir -p "$OUT_DIR"

if [[ "$PRIMARY_PREFER" != "primary" && "$PRIMARY_PREFER" != "aligned" ]]; then
  die "--primary-prefer must be 'primary' or 'aligned'"
fi

command -v samtools >/dev/null 2>&1 || die "samtools not found"
command -v bedtools >/dev/null 2>&1 || die "bedtools not found (needed for getfasta)"
command -v bedGraphToBigWig >/dev/null 2>&1 || die "bedGraphToBigWig not found"
command -v python3 >/dev/null 2>&1 || die "python3 not found"
command -v "$R_EXE" >/dev/null 2>&1 || die "R executable not found in PATH: $R_EXE"

MAPI_ROOT="${MAPI_ROOT:-$HOME/MalariAPI}"

# Hidden R helper (primary only)
R_HELPER="$MAPI_ROOT/bin/modules/.CLOverCNV/.cov_model_primary_factor.R"
[[ -f "$R_HELPER" ]] || die "Missing helper: $R_HELPER"

# ---------------------------
# Helpers: filename sanitizing
# ---------------------------
safe_label() {
  local x="${1:-}"
  x="$(basename "$x")"
  x="${x// /_}"
  x="${x//\//_}"
  x="${x//\\/_}"
  x="${x//:/_}"
  x="${x//;/_}"
  x="${x//|/_}"
  x="${x//,/_}"
  x="${x//(/_}"
  x="${x//)/_}"
  x="${x//[/_}"
  x="${x//]/_}"
  x="${x//\{/_}"
  x="${x//\}/_}"
  x="${x//+/_}"
  x="${x//=/_}"
  x="${x//#/_}"
  x="${x//@/_}"
  x="${x//\?/_}"
  x="${x//\*/_}"
  x="${x//\"/_}"
  x="${x//\'/_}"
  while [[ "$x" == *"__"* ]]; do
    x="${x//__/_}"
  done
  x="${x#_}"
  x="${x%_}"
  printf "%s" "$x"
}

# ---------------------------
# FASTA resolve
# ---------------------------
say "Resolving FASTA for genome spec '$GENOME' under $MAPI_ROOT/genomes/..."
FASTA=""

if [[ -f "$GENOME" ]]; then
  FASTA="$GENOME"
fi

if [[ -z "$FASTA" ]]; then
  FASTA="$(find "$MAPI_ROOT/genomes" -type f \
    \( -iname "${GENOME}.fa" -o -iname "${GENOME}.fasta" -o -iname "${GENOME}.fna" \
       -o -iname "${GENOME}.fa.gz" -o -iname "${GENOME}.fasta.gz" -o -iname "${GENOME}.fna.gz" \) \
    ! -iname "*.fai" ! -iname "*.gzi" \
    2>/dev/null | head -n 1 || true)"
fi

[[ -n "$FASTA" && -f "$FASTA" ]] || die "Could not resolve FASTA for genome '$GENOME'."

# Use only for filenames
GENOME_LABEL="$(basename "$FASTA")"
GENOME_LABEL="${GENOME_LABEL%.gz}"
GENOME_LABEL="${GENOME_LABEL%.fasta}"
GENOME_LABEL="${GENOME_LABEL%.fa}"
GENOME_LABEL="${GENOME_LABEL%.fna}"
GENOME_LABEL="$(safe_label "$GENOME_LABEL")"

[[ -n "$GENOME_LABEL" ]] || die "Could not derive genome label from FASTA: $FASTA"

FASTA_RESOLVED="$FASTA"
FASTA_TMP_COPY=""
if [[ "$FASTA" =~ \.gz$ ]]; then
  say "FASTA is gzipped; creating decompressed copy for this run."
  FASTA_TMP_COPY="$OUT_DIR/.tmp_${GENOME_LABEL}_fa_$$.fa"
  gzip -cd "$FASTA" > "$FASTA_TMP_COPY"
  FASTA_RESOLVED="$FASTA_TMP_COPY"
fi

FASTA="$FASTA_RESOLVED"
FAI="${FASTA}.fai"
if [[ ! -f "$FAI" ]]; then
  say "Indexing FASTA: samtools faidx $FASTA"
  samtools faidx "$FASTA"
fi

cleanup_fa() {
  if [[ -n "$FASTA_TMP_COPY" && -f "$FASTA_TMP_COPY" && "$KEEP_TMP" -ne 1 ]]; then
    rm -f "$FASTA_TMP_COPY" "${FASTA_TMP_COPY}.fai" 2>/dev/null || true
  fi
}
trap cleanup_fa EXIT

# Chrom sizes
CHROM_SIZES="$OUT_DIR/${GENOME_LABEL}.chrom.sizes"
cut -f1,2 "$FAI" > "$CHROM_SIZES"
[[ -s "$CHROM_SIZES" ]] || die "Failed to write chrom.sizes: $CHROM_SIZES"

# ---------------------------
# Helpers
# ---------------------------
run_r() {
  local script="$1"; shift
  unset R_LIBS_USER 2>/dev/null || true

  if [[ -n "$R_MODULES" ]]; then
    bash -lc "module purge >/dev/null 2>&1 || true; module load $R_MODULES; \"$R_EXE\" \"$script\" $*"
  else
    "$R_EXE" "$script" "$@"
  fi
}

extract_sample_from_filename() {
  local base
  base="$(basename "$1")"
  base="${base%.bam}"
  echo "${base%%.*}"
}

detect_tag_from_filename() {
  local base lc
  base="$(basename "$1")"
  lc="$(echo "$base" | tr '[:upper:]' '[:lower:]')"

  local ordered=(unmapped supplementary discordant splitters secondary primary aligned)
  for t in "${ordered[@]}"; do
    if [[ "$lc" == *"$t"* ]]; then
      echo "$t"
      return 0
    fi
    if [[ "$t" == "splitters" && "$lc" == *"split"* ]]; then
      echo "splitters"
      return 0
    fi
  done

  echo "aligned"
}

bam_mapped_reads() {
  local bam="$1"
  samtools view -c -@ "$THREADS" -F 4 "$bam" 2>/dev/null || echo 0
}

build_bins_bed() {
  local out_bed="$1"
  awk -v OFS='\t' -v B="$BINSIZE" '{for (s=0; s<$2; s+=B){e=s+B; if(e>$2)e=$2; print $1,s,e}}' "$CHROM_SIZES" > "$out_bed"
}

compute_at() {
  local bins_bed="$1"
  local out_at_tsv="$2"
  local tmpdir="$3"

  bedtools getfasta -fi "$FASTA" -bed "$bins_bed" -fo "$tmpdir/bins.fa" >/dev/null

  python3 - <<'PY' "$tmpdir/bins.fa" "$out_at_tsv"
import sys, re
fa, out = sys.argv[1], sys.argv[2]

def at_frac(seq):
    seq = seq.upper()
    if not seq:
        return 0.0
    return (seq.count("A") + seq.count("T")) / len(seq)

with open(fa) as f, open(out, "w") as o:
    o.write("chrom\tstart\tend\tAT\n")
    chrom = start = end = None
    seq = []
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            if chrom is not None:
                s = "".join(seq)
                o.write(f"{chrom}\t{start}\t{end}\t{at_frac(s):.6f}\n")
            seq = []
            m = re.match(r">([^:]+):(\d+)-(\d+)", line)
            if not m:
                raise SystemExit(f"Unexpected FASTA header: {line}")
            chrom = m.group(1)
            start = int(m.group(2))
            end = int(m.group(3))
        else:
            seq.append(line)
    if chrom is not None:
        s = "".join(seq)
        o.write(f"{chrom}\t{start}\t{end}\t{at_frac(s):.6f}\n")
PY
}

resolve_mappability() {
  local bins_bed="$1"
  local out_map_tsv="$2"

  if [[ "$NO_MAPPABILITY" -eq 1 ]]; then
    awk -v OFS='\t' 'BEGIN{print "chrom\tstart\tend\tM"} {print $1,$2,$3,1.0}' "$bins_bed" > "$out_map_tsv"
    return 0
  fi

  [[ -f "$MAP_TRACK" ]] || die "Mappability enabled but track not found: $MAP_TRACK"

  python3 - <<'PY' "$MAP_TRACK" "$out_map_tsv"
import sys
inp, out = sys.argv[1], sys.argv[2]

with open(inp) as f:
    first = f.readline().rstrip("\n")
    parts = first.split("\t")
    has_header = (len(parts) >= 4) and (
        parts[0].lower() in ("chrom", "chr")
        or parts[-1].lower().startswith("map")
        or parts[-1].lower() == "m"
    )

with open(inp) as f, open(out, "w") as o:
    o.write("chrom\tstart\tend\tM\n")
    if has_header:
        next(f)
    for line in f:
        if not line.strip():
            continue
        p = line.rstrip("\n").split("\t")
        if len(p) < 4:
            raise SystemExit("mappability track must have 4 columns: chrom start end mappability")
        o.write(f"{p[0]}\t{p[1]}\t{p[2]}\t{p[3]}\n")
PY
}

compute_mean_depth_bins() {
  local bam="$1"
  local out_bg="$2"
  local tmpdir="$3"

  local perbase="$tmpdir/depth.perbase.tsv"
  samtools depth -@ "$THREADS" -aa "$bam" > "$perbase"

  awk -v OFS='\t' -v B="$BINSIZE" '
    BEGIN{chr=""; bin_s=0; bin_e=B; sum=0; n=0;}
    function flush() {
      if(chr!=""){
        mean=(n>0 ? sum/n : 0);
        print chr, bin_s, bin_e, mean;
      }
      sum=0; n=0;
    }
    function adv(pos0){
      while(pos0 >= bin_e){
        flush();
        bin_s=bin_e; bin_e=bin_s+B;
      }
    }
    {
      c=$1; pos1=$2; d=$3;
      pos0=pos1-1;
      if(chr=="" || c!=chr){
        if(chr!="") flush();
        chr=c; bin_s=0; bin_e=B; sum=0; n=0;
      }
      adv(pos0);
      sum += d;
      n += 1;
    }
    END{ if(chr!="") flush(); }
  ' "$perbase" > "$out_bg"

  python3 - <<'PY' "$out_bg" "$CHROM_SIZES"
import sys, os
bg, cs = sys.argv[1], sys.argv[2]
sizes = {}
with open(cs) as f:
    for line in f:
        if not line.strip():
            continue
        c, l = line.rstrip("\n").split("\t")[:2]
        sizes[c] = int(l)

tmp = bg + ".tmp"
with open(bg) as inp, open(tmp, "w") as out:
    for line in inp:
        c, s, e, v = line.rstrip("\n").split("\t")
        s = int(float(s))
        e = int(float(e))
        L = sizes.get(c)
        if L is None:
            continue
        if s >= L:
            continue
        if e > L:
            e = L
        if e <= s:
            continue
        out.write(f"{c}\t{s}\t{e}\t{v}\n")
os.replace(tmp, bg)
PY

  sort -k1,1 -k2,2n "$out_bg" -o "$out_bg"
}

join_features() {
  local depth_bg="$1"
  local at_tsv="$2"
  local map_tsv="$3"
  local out_features="$4"

  python3 - <<'PY' "$depth_bg" "$at_tsv" "$map_tsv" "$out_features"
import sys
depth_bg, at_tsv, map_tsv, out = sys.argv[1:]

def load_keyed(path, want_cols):
    d = {}
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        for line in f:
            if not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            key = (p[idx["chrom"]], p[idx["start"]], p[idx["end"]])
            d[key] = [p[idx[c]] for c in want_cols]
    return d

depth = {}
with open(depth_bg) as f:
    for line in f:
        if not line.strip():
            continue
        c, s, e, v = line.rstrip("\n").split("\t")
        depth[(c, s, e)] = v

at = load_keyed(at_tsv, ["AT"])
mp = load_keyed(map_tsv, ["M"])

with open(out, "w") as o:
    o.write("chrom\tstart\tend\tdepth\tAT\tM\n")
    miss = 0
    for key, v in depth.items():
        if key not in at or key not in mp:
            miss += 1
            continue
        c, s, e = key
        o.write(f"{c}\t{s}\t{e}\t{v}\t{at[key][0]}\t{mp[key][0]}\n")
    if miss:
        print(f"[warn] skipped {miss} bins missing AT/M", file=sys.stderr)
PY
}

bedgraph_to_bw() {
  local in_bg="$1"
  local out_bw="$2"
  bedGraphToBigWig "$in_bg" "$CHROM_SIZES" "$out_bw"
}

sanitize_bedgraph_numeric() {
  local bg="$1"
  [[ -s "$bg" ]] || return 0
  local bad
  bad="$(awk '$4=="NA" || $4=="" || $4=="nan" || $4=="NaN" || $4=="inf" || $4=="-inf"{c++} END{print c+0}' "$bg")"
  if [[ "$bad" -gt 0 ]]; then
    say "NOTE: Dropping $bad non-numeric rows from bedGraph: $bg"
    local tmp="${bg}.clean"
    awk 'BEGIN{OFS="\t"} $4!="NA" && $4!="" && $4!="nan" && $4!="NaN" && $4!="inf" && $4!="-inf"{print $1,$2,$3,$4}' "$bg" > "$tmp"
    mv "$tmp" "$bg"
  fi
}

sanitize_4col_numeric() {
  local f="$1"
  [[ -s "$f" ]] || return 0
  local bad
  bad="$(awk '$4=="NA" || $4=="" || $4=="nan" || $4=="NaN" || $4=="inf" || $4=="-inf"{c++} END{print c+0}' "$f")"
  if [[ "$bad" -gt 0 ]]; then
    say "NOTE: Dropping $bad non-numeric rows from 4-col file: $f"
    local tmp="${f}.clean"
    awk 'BEGIN{OFS="\t"} $4!="NA" && $4!="" && $4!="nan" && $4!="NaN" && $4!="inf" && $4!="-inf"{print $1,$2,$3,$4}' "$f" > "$tmp"
    mv "$tmp" "$f"
  fi
}

apply_cfactor_to_bg() {
  local raw_bg="$1"
  local cfactor_tsv="$2"
  local out_bg="$3"

  python3 - <<'PY' "$raw_bg" "$cfactor_tsv" "$out_bg"
import sys
raw_bg, cf_tsv, out_bg = sys.argv[1:]

cf = {}
with open(cf_tsv) as f:
    for line in f:
        if not line.strip():
            continue
        c, s, e, x = line.rstrip("\n").split("\t")
        cf[(c, s, e)] = float(x)

with open(raw_bg) as inp, open(out_bg, "w") as out:
    for line in inp:
        if not line.strip():
            continue
        c, s, e, v = line.rstrip("\n").split("\t")
        mult = cf.get((c, s, e), 1.0)
        out.write(f"{c}\t{s}\t{e}\t{float(v) * mult}\n")
PY

  sort -k1,1 -k2,2n "$out_bg" -o "$out_bg"
}

select_primary_bam() {
  local list_file="$1"

  local want1="$PRIMARY_PREFER"
  local want2="aligned"
  [[ "$want1" == "aligned" ]] && want2="primary"

  local bam1 bam2 any
  bam1="$(awk -F'\t' -v t="$want1" '$1==t{print $2; exit}' "$list_file" || true)"
  bam2="$(awk -F'\t' -v t="$want2" '$1==t{print $2; exit}' "$list_file" || true)"
  any="$(awk -F'\t' 'NF>=2{print $2; exit}' "$list_file" || true)"

  if [[ -n "$bam1" && -f "$bam1" ]]; then
    echo "$bam1"
  elif [[ -n "$bam2" && -f "$bam2" ]]; then
    echo "$bam2"
  elif [[ -n "$any" && -f "$any" ]]; then
    echo "$any"
  else
    echo ""
  fi
}

# ---------------------------
# Gather BAMs
# ---------------------------
say "Scanning BAMs in: $DIR"
mapfile -d '' ALL_BAMS < <(find "$DIR" -maxdepth 1 -type f -name "*.bam" -print0)
[[ "${#ALL_BAMS[@]}" -gt 0 ]] || die "No .bam files found in: $DIR"

declare -A SAMPLE_TO_ROWS
for bam in "${ALL_BAMS[@]}"; do
  bn="$(basename "$bam")"

  if [[ -n "$SAMPLE_KEY" ]]; then
    [[ "$bn" == *"$SAMPLE_KEY"* ]] || continue
    sample="$SAMPLE_KEY"
  else
    sample="$(extract_sample_from_filename "$bam")"
  fi

  tag="$(detect_tag_from_filename "$bam")"
  SAMPLE_TO_ROWS["$sample"]+="${tag}"$'\t'"${bam}"$'\n'
done

# ---------------------------
# Process each sample
# ---------------------------
for sample in "${!SAMPLE_TO_ROWS[@]}"; do
  say "Sample: $sample"

  SAMPLE_LABEL="$(safe_label "$sample")"

  sample_dir="$OUT_DIR/$SAMPLE_LABEL"
  corr_dir="$sample_dir/corr_bw"
  raw_dir="$sample_dir/raw_bw"
  fac_dir="$sample_dir/factors"
  tmp_dir="$sample_dir/.tmp"
  mkdir -p "$corr_dir" "$fac_dir" "$tmp_dir"
  if [[ "$EMIT_RAW" -eq 1 ]]; then
    mkdir -p "$raw_dir"
  fi

  list_file="$tmp_dir/bams.tsv"
  printf "%s" "${SAMPLE_TO_ROWS[$sample]}" | awk -F'\t' 'NF>=2{print $1"\t"$2}' > "$list_file"

  primary_bam="$(select_primary_bam "$list_file")"
  [[ -n "$primary_bam" ]] || { say "No BAMs usable for sample=$sample; skipping."; continue; }

  n_mapped="$(bam_mapped_reads "$primary_bam")"
  if [[ "$n_mapped" -eq 0 ]]; then
    say "Primary candidate has 0 mapped reads; skipping sample=$sample: $primary_bam"
    continue
  fi
  say "Primary BAM for factors: $primary_bam"

  bins_bed="$tmp_dir/${GENOME_LABEL}.bins${BINSIZE}.bed"
  build_bins_bed "$bins_bed"

  at_tsv="$tmp_dir/${SAMPLE_LABEL}.at.tsv"
  map_tsv="$tmp_dir/${SAMPLE_LABEL}.map.tsv"
  compute_at "$bins_bed" "$at_tsv" "$tmp_dir"
  resolve_mappability "$bins_bed" "$map_tsv"

  primary_raw_bg="$tmp_dir/${SAMPLE_LABEL}.primary.raw.bins${BINSIZE}.bedgraph"
  compute_mean_depth_bins "$primary_bam" "$primary_raw_bg" "$tmp_dir"

  features_tsv="$tmp_dir/${SAMPLE_LABEL}.primary.features.tsv"
  join_features "$primary_raw_bg" "$at_tsv" "$map_tsv" "$features_tsv"

  expected_tsv="$tmp_dir/${SAMPLE_LABEL}.expected.bins${BINSIZE}.tsv"
  cfactor_tsv="$tmp_dir/${SAMPLE_LABEL}.cfactor.bins${BINSIZE}.tsv"

  say "Computing expected + cfactor via model..."
  run_r "$R_HELPER" \
    --features "$features_tsv" \
    --model "$MODEL_RDS" \
    --out-expected "$expected_tsv" \
    --out-cfactor "$cfactor_tsv"

  sanitize_4col_numeric "$expected_tsv"
  sanitize_4col_numeric "$cfactor_tsv"

  expected_bg="$tmp_dir/${SAMPLE_LABEL}.expected.bedgraph"
  cfactor_bg="$tmp_dir/${SAMPLE_LABEL}.cfactor.bedgraph"
  cp "$expected_tsv" "$expected_bg"
  cp "$cfactor_tsv" "$cfactor_bg"
  sort -k1,1 -k2,2n "$expected_bg" -o "$expected_bg"
  sort -k1,1 -k2,2n "$cfactor_bg" -o "$cfactor_bg"

  sanitize_bedgraph_numeric "$expected_bg"
  sanitize_bedgraph_numeric "$cfactor_bg"

  bedgraph_to_bw "$expected_bg" "$fac_dir/${SAMPLE_LABEL}.expected.bins${BINSIZE}.bw"
  bedgraph_to_bw "$cfactor_bg" "$fac_dir/${SAMPLE_LABEL}.cfactor.bins${BINSIZE}.bw"

  while IFS=$'\t' read -r tag bam; do
    [[ -n "$bam" && -f "$bam" ]] || continue

    TAG_LABEL="$(safe_label "$tag")"

    n_mapped="$(bam_mapped_reads "$bam")"
    if [[ "$n_mapped" -eq 0 ]]; then
      say "Skipping BAM with 0 mapped reads (tag=$tag): $bam"
      continue
    fi

    say "Processing BAM tag=$tag bam=$bam"

    raw_bg="$tmp_dir/${SAMPLE_LABEL}.${TAG_LABEL}.raw.bins${BINSIZE}.bedgraph"
    compute_mean_depth_bins "$bam" "$raw_bg" "$tmp_dir"

    corr_bg="$tmp_dir/${SAMPLE_LABEL}.${TAG_LABEL}.corr.bins${BINSIZE}.bedgraph"
    apply_cfactor_to_bg "$raw_bg" "$cfactor_tsv" "$corr_bg"

    sanitize_bedgraph_numeric "$raw_bg"
    sanitize_bedgraph_numeric "$corr_bg"

    bedgraph_to_bw "$corr_bg" "$corr_dir/${SAMPLE_LABEL}.${TAG_LABEL}.corr.bins${BINSIZE}.bw"
    if [[ "$EMIT_RAW" -eq 1 ]]; then
      bedgraph_to_bw "$raw_bg" "$raw_dir/${SAMPLE_LABEL}.${TAG_LABEL}.raw.bins${BINSIZE}.bw"
    fi

    rm -f "$raw_bg" "$corr_bg" "$tmp_dir/depth.perbase.tsv" 2>/dev/null || true
  done < "$list_file"

  if [[ "$KEEP_TMP" -ne 1 ]]; then
    rm -rf "$tmp_dir"
  else
    say "Keeping tmp for sample: $tmp_dir"
  fi

  say "Sample done: $sample -> $sample_dir"
done

say "All done. Output root: $OUT_DIR"
