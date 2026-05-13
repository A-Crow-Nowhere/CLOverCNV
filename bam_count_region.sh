#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  bam_region_summary.sh \
    --bam-dir /path/to/bams \
    --region chr:start-end \
    --out summary.tsv \
    [--threads 1]

Description:
  Loops through every BAM in a directory and summarizes:
    - reads overlapping the target region
    - mean depth in the target region
    - mapped reads on the full chromosome
    - mean depth across the full chromosome
    - region/chromosome depth ratio

Requirements:
  - samtools
  - BAMs should be indexed (.bai)

Example:
  bam_region_summary.sh \
    --bam-dir ./bams \
    --region PfDd2_11:100000-120000 \
    --out region_summary.tsv
EOF
}

die() {
  echo "ERROR: $*" >&2
  exit 1
}

BAM_DIR=""
REGION=""
OUT=""
THREADS=1

while [[ $# -gt 0 ]]; do
  case "$1" in
    --bam-dir) BAM_DIR="${2:-}"; shift 2 ;;
    --region)  REGION="${2:-}"; shift 2 ;;
    --out)     OUT="${2:-}"; shift 2 ;;
    --threads) THREADS="${2:-1}"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown argument: $1" ;;
  esac
done

[[ -n "$BAM_DIR" ]] || die "--bam-dir is required"
[[ -d "$BAM_DIR" ]] || die "BAM directory not found: $BAM_DIR"
[[ -n "$REGION" ]]  || die "--region is required"
[[ -n "$OUT" ]]     || die "--out is required"

if [[ ! "$REGION" =~ ^([^:]+):([0-9]+)-([0-9]+)$ ]]; then
  die "Region must look like chr:start-end"
fi

CHROM="${BASH_REMATCH[1]}"
START="${BASH_REMATCH[2]}"
END="${BASH_REMATCH[3]}"

if (( START > END )); then
  die "Region start is greater than end"
fi

command -v samtools >/dev/null 2>&1 || die "samtools not found in PATH"

shopt -s nullglob
bams=( "$BAM_DIR"/*.bam )
(( ${#bams[@]} > 0 )) || die "No BAM files found in: $BAM_DIR"

printf "bam_file\tsample\tregion\tchromosome\tregion_read_count\tregion_mean_depth\tchrom_read_count\tchrom_mean_depth\tregion_vs_chrom_depth_ratio\n" > "$OUT"

for bam in "${bams[@]}"; do
  bam_base="$(basename "$bam")"
  sample="${bam_base%.bam}"

  bai1="${bam}.bai"
  bai2="${bam%.bam}.bai"
  if [[ ! -f "$bai1" && ! -f "$bai2" ]]; then
    echo "[INFO] Indexing $bam_base" >&2
    samtools index -@ "$THREADS" "$bam"
  fi

  # Check that chromosome exists in BAM
  if ! samtools idxstats "$bam" | cut -f1 | grep -Fxq "$CHROM"; then
    echo "[WARN] $bam_base does not contain chromosome $CHROM; writing NA row" >&2
    printf "%s\t%s\t%s\t%s\tNA\tNA\tNA\tNA\tNA\n" \
      "$bam_base" "$sample" "$REGION" "$CHROM" >> "$OUT"
    continue
  fi

  # 1) Count alignment records overlapping the region
  region_read_count="$(samtools view -@ "$THREADS" -c "$bam" "$REGION")"

  # 2) Mean depth in region, including zero-covered positions
  region_mean_depth="$(
    samtools depth -aa -r "$REGION" "$bam" | \
      awk 'BEGIN{s=0;n=0} {s+=$3; n++} END{if(n>0) printf "%.6f", s/n; else print "NA"}'
  )"

  # 3) Number of mapped reads/alignment records on full chromosome
  chrom_read_count="$(
    samtools idxstats "$bam" | \
      awk -v c="$CHROM" '$1==c {print $3; found=1} END{if(!found) print "NA"}'
  )"

  # 4) Mean depth across full chromosome, including zero-covered positions
  chrom_mean_depth="$(
    samtools depth -aa -r "$CHROM" "$bam" | \
      awk 'BEGIN{s=0;n=0} {s+=$3; n++} END{if(n>0) printf "%.6f", s/n; else print "NA"}'
  )"

  # 5) Ratio of region mean depth to chromosome mean depth
  ratio="$(
    awk -v r="$region_mean_depth" -v c="$chrom_mean_depth" '
      BEGIN {
        if (r=="NA" || c=="NA" || c==0) print "NA";
        else printf "%.6f", r/c
      }'
  )"

  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$bam_base" "$sample" "$REGION" "$CHROM" \
    "$region_read_count" "$region_mean_depth" \
    "$chrom_read_count" "$chrom_mean_depth" "$ratio" >> "$OUT"
done

echo "[INFO] Wrote: $OUT" >&2
