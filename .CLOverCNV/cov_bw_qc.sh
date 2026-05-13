#!/usr/bin/env bash
set -euo pipefail

# ------------------------------------------------------------
# cov_bw_qc: QC compare raw vs corrected/expected bigWigs per sample
# + optional dedup bigWig generation (requires BAM w/ dup marks)
# ------------------------------------------------------------

print_help() {
  cat <<'EOF'
Usage:
  mapi modules cov_bw_qc --outdir DIR [--out TSV] [--pattern REGEX] [--skip-missing]
                         [--no-dedup] [--bam-root DIR] [--bam-pattern GLOB] [--eps X]

What it does:
  - Scans DIR for sample subdirectories.
  - Finds raw and corrected bigWig directories inside each sample.
  - Pairs *.corr.*.bw with corresponding *.raw.*.bw (same basename except corr/raw token).
  - Uses bedtools unionbedg + awk to compute mean/sd/CV for raw and corrected, and residual stats.
  - Writes a TSV summary with:
      - flag_corr_more_variant: TRUE if corrected/expected is more variant than raw (CV comparison)
      - stability scores:
          * stability_eff_sd = resid_sd / sqrt(bin_size)
          * stability_eff_cv = (resid_sd / sqrt(bin_size)) / raw_mean
  - By default, also attempts to create deduplicated depth bigWigs:
      * dedup_raw  = raw - dup_depth
      * dedup_corr = corr * (dedup_raw / (raw + eps))
    and writes duplicate rows in TSV with de_dup=TRUE.

Arguments:
  --outdir DIR        Parent dir containing /<sample>/... subdirs (required)
  --out TSV           Output TSV path (default: <outdir>/cov_bw_qc.summary.tsv)
  --pattern REGEX     Only evaluate corrected files whose basename matches REGEX (default: ".*\\.bw$")
  --skip-missing      Skip pairs where either raw or corrected is missing (default: error row with status)

Dedup options (default ON):
  --no-dedup          Disable dedup bigWig generation and dedup TSV rows
  --bam-root DIR      Where to search for BAMs (default: sample_dir)
  --bam-pattern GLOB  BAM filename glob (default: "*.bam")
  --eps X             Small constant for ratios (default: 1e-6)

Notes:
  - Requires: bigWigToBedGraph, bigWigInfo, bedGraphToBigWig (UCSC tools) and bedtools, samtools
  - bin_size inferred from filename token ".bins<INT>." (e.g. ".bins50.")
EOF
}

# ------------------------------------------------------------
# Args
# ------------------------------------------------------------
OUTDIR=""
OUT=""
PATTERN='.*\.bw$'
SKIP_MISSING=0
DO_DEDUP=1
BAM_ROOT=""
BAM_PATTERN="*.bam"
EPS="1e-6"

if [[ $# -eq 0 ]]; then
  print_help
  exit 1
fi

while [[ $# -gt 0 ]]; do
  case "$1" in
    --outdir) OUTDIR="${2:-}"; shift 2;;
    --out) OUT="${2:-}"; shift 2;;
    --pattern) PATTERN="${2:-}"; shift 2;;
    --skip-missing) SKIP_MISSING=1; shift 1;;
    --no-dedup) DO_DEDUP=0; shift 1;;
    --bam-root) BAM_ROOT="${2:-}"; shift 2;;
    --bam-pattern) BAM_PATTERN="${2:-}"; shift 2;;
    --eps) EPS="${2:-1e-6}"; shift 2;;
    -h|--help) print_help; exit 0;;
    *)
      echo "[cov_bw_qc] ERROR: Unknown option: $1" >&2
      print_help >&2
      exit 2
      ;;
  esac
done

# ------------------------------------------------------------
# Checks
# ------------------------------------------------------------
if [[ -z "${OUTDIR}" ]]; then
  echo "[cov_bw_qc] ERROR: --outdir is required" >&2
  exit 2
fi
if [[ ! -d "${OUTDIR}" ]]; then
  echo "[cov_bw_qc] ERROR: --outdir does not exist: ${OUTDIR}" >&2
  exit 2
fi

for tool in bigWigToBedGraph bedtools; do
  command -v "$tool" >/dev/null 2>&1 || { echo "[cov_bw_qc] ERROR: $tool not found in PATH" >&2; exit 2; }
done

if [[ "$DO_DEDUP" -eq 1 ]]; then
  for tool in bigWigInfo bedGraphToBigWig samtools; do
    command -v "$tool" >/dev/null 2>&1 || { echo "[cov_bw_qc] ERROR: $tool not found in PATH (needed for dedup)" >&2; exit 2; }
  done
fi

if [[ -z "${OUT}" ]]; then
  OUT="${OUTDIR%/}/cov_bw_qc.summary.tsv"
fi
mkdir -p "$(dirname "$OUT")"

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------
_find_subdir() {
  local base="$1"; shift
  local cand
  for cand in "$@"; do
    if [[ -d "${base}/${cand}" ]]; then
      echo "${base}/${cand}"
      return 0
    fi
  done
  return 1
}

_guess_raw_from_corr() {
  local corr="$1"
  local raw=""
  raw="${corr/.corr./.raw.}";        [[ "$raw" != "$corr" ]] && { echo "$raw"; return 0; }
  raw="${corr/.corrected./.raw.}";   [[ "$raw" != "$corr" ]] && { echo "$raw"; return 0; }
  raw="${corr/.expected./.raw.}";    [[ "$raw" != "$corr" ]] && { echo "$raw"; return 0; }
  echo ""
}

_infer_bin_size() {
  local name="$1"
  awk -v s="$name" 'BEGIN{
    if (match(s, /[.]bins([0-9]+)[.]/, m)) print m[1];
    else print "NA";
  }'
}

# Extract chrom sizes from a bigWig (no external genome.sizes needed)
_bw_to_chrom_sizes() {
  local bw="$1"
  local out="$2"

  # bigWigInfo -chroms output varies across UCSC tool builds.
  # Robust rule: keep lines whose LAST field is an integer size,
  # and take the FIRST field as chrom name.
  bigWigInfo -chroms "$bw" 2>/dev/null \
  | awk '($NF ~ /^[0-9]+$/) {print $1 "\t" $NF}' \
  | sed 's/:$//' \
  > "$out"

  [[ -s "$out" ]] || return 1
}


_union_stats() {
  local raw_bw="$1"
  local corr_bw="$2"

  bedtools unionbedg \
    -i <(bigWigToBedGraph "$raw_bw" /dev/stdout) \
       <(bigWigToBedGraph "$corr_bw" /dev/stdout) \
    -filler 0 \
  | awk '
      {r=$4; c=$5;
       sR+=r; ssR+=r*r;
       sC+=c; ssC+=c*c;
       d=r-c; sD+=d; ssD+=d*d;
       n++}
      END{
        if(n==0){ print "0\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA"; exit }
        mR=sR/n; vR=ssR/n-mR*mR; sdR=(vR>0?sqrt(vR):0); cvR=(mR!=0? sdR/mR : "NA");
        mC=sC/n; vC=ssC/n-mC*mC; sdC=(vC>0?sqrt(vC):0); cvC=(mC!=0? sdC/mC : "NA");
        mD=sD/n; vD=ssD/n-mD*mD; sdD=(vD>0?sqrt(vD):0);
        printf "%d\t%.6g\t%.6g\t%s\t%.6g\t%.6g\t%s\t%.6g\t%.6g\t%.6g\n", n, mR, sdR, cvR, mC, sdC, cvC, mD, sdD, vD
      }'
}

# Find a BAM for a given sample+track (best-effort)
_find_bam_for_track() {
  local sample_dir="$1"
  local track_prefix="$2"   # e.g. "C01_1.primary" or "C7B.primary"
  local search_dir="$3"

  # If user didn't set --bam-root, default to sample_dir
  [[ -z "$search_dir" ]] && search_dir="$sample_dir"
  [[ -d "$search_dir" ]] || { echo ""; return 0; }

  # Prefer BAMs matching the track prefix
  local cand
  shopt -s nullglob
  local matches=("$search_dir"/"${track_prefix}"*.bam)
  if [[ ${#matches[@]} -gt 0 ]]; then
    echo "${matches[0]}"
    return 0
  fi

  # Otherwise, if there's exactly one BAM in the directory, use it
  local all=("$search_dir"/$BAM_PATTERN)
  if [[ ${#all[@]} -eq 1 ]]; then
    echo "${all[0]}"
    return 0
  fi

  # Otherwise none / ambiguous
  echo ""
  return 0
}

# Build dedup bigWigs next to existing bw files.
# Requires a BAM with duplicate marks (flag 1024).
_make_dedup_bigwigs() {
  local raw_bw="$1"
  local corr_bw="$2"
  local bam="$3"
  local eps="$4"
  local out_dedup_raw_bw="$5"
  local out_dedup_corr_bw="$6"

  local tmpd; tmpd="$(mktemp -d)"
  trap 'rm -rf "$tmpd"' RETURN

  local chrom_sizes="$tmpd/chrom.sizes"
  _bw_to_chrom_sizes "$raw_bw" "$chrom_sizes" || {
    echo "[cov_bw_qc] WARN: Failed to derive chrom sizes from $raw_bw" >&2
    return 1
  }

  # Use raw bigWig bedGraph as the canonical bin set (chr,start,end)
  # We'll compute dup-only mean depth per those bins from the BAM using bedtools coverage.
  local bins_bed="$tmpd/bins.bed"
 bigWigToBedGraph "$raw_bw" /dev/stdout \
 | awk 'BEGIN{OFS="\t"}{print $1,$2,$3}' \
 | sort -k1,1 -k2,2n > "$bins_bed"
 

  # Make a dup-only BAM (temporary) so bedtools can read it reliably
  local dup_bam="$tmpd/dup.bam"
  samtools view -bh -f 1024 "$bam" > "$dup_bam"
  samtools index "$dup_bam"

  # dup mean depth per bin: bedtools coverage -mean returns a BED-like line with mean in last column.
  # Columns: chrom start end ... mean
  local dup_bg="$tmpd/dup.bg"
  bedtools coverage -a "$bins_bed" -b "$dup_bam" -mean \
  | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$NF}' > "$dup_bg"

  # dedup_raw = raw - dup (floored)
  local dedup_raw_bg="$tmpd/dedup_raw.bg"
  bedtools unionbedg -i \
    <(bigWigToBedGraph "$raw_bw" /dev/stdout) \
    "$dup_bg" \
    -filler 0 \
  | awk 'BEGIN{OFS="\t"}{v=$4-$5; if(v<0) v=0; print $1,$2,$3,v}' \
  | sort -k1,1 -k2,2n > "$dedup_raw_bg"

  bedGraphToBigWig "$dedup_raw_bg" "$chrom_sizes" "$out_dedup_raw_bw"

  # dedup_corr = corr * (dedup_raw / (raw+eps))
  local dedup_corr_bg="$tmpd/dedup_corr.bg"
  bedtools unionbedg -i \
    <(bigWigToBedGraph "$corr_bw" /dev/stdout) \
    <(bigWigToBedGraph "$raw_bw" /dev/stdout) \
    "$dedup_raw_bg" \
    -filler 0 \
  | awk -v eps="$eps" 'BEGIN{OFS="\t"}
      {corr=$4; raw=$5; ded=$6;
       f = (raw>0 ? ded/(raw+eps) : 0);
       v = corr * f;
       if(v<0) v=0;
       print $1,$2,$3,v
      }' \
  | sort -k1,1 -k2,2n > "$dedup_corr_bg"

  bedGraphToBigWig "$dedup_corr_bg" "$chrom_sizes" "$out_dedup_corr_bw"

  return 0
}

# ------------------------------------------------------------
# Main
# ------------------------------------------------------------

echo -e "sample\ttrack\tfile_type\tbin_size\tde_dup\traw_bw\tcorr_bw\tn\traw_mean\traw_sd\traw_cv\tcorr_mean\tcorr_sd\tcorr_cv\tresid_mean\tresid_sd\tresid_var\tstability_eff_sd\tstability_eff_cv\tflag_corr_more_variant\tstatus" > "$OUT"

shopt -s nullglob

for sample_dir in "${OUTDIR%/}/"*; do
  [[ -d "$sample_dir" ]] || continue
  sample="$(basename "$sample_dir")"

  raw_dir="$(_find_subdir "$sample_dir" raw_bw raw rawBW 2>/dev/null || true)"
  corr_dir="$(_find_subdir "$sample_dir" corr_bw corrected_bw corrected corr expected_bw expected corrBW 2>/dev/null || true)"

  if [[ -z "$raw_dir" || -z "$corr_dir" ]]; then
    echo -e "${sample}\t.\t.\tNA\tFALSE\t.\t.\t0\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tmissing_raw_or_corr_dir" >> "$OUT"
    continue
  fi

  mapfile -t corr_bws < <(find "$corr_dir" -maxdepth 1 -type f -name "*.bw" | sort)
  if [[ "${#corr_bws[@]}" -eq 0 ]]; then
    echo -e "${sample}\t.\t.\tNA\tFALSE\t.\t.\t0\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tno_corr_bw_found" >> "$OUT"
    continue
  fi

  for corr_bw in "${corr_bws[@]}"; do
    bn="$(basename "$corr_bw")"
    [[ "$bn" =~ $PATTERN ]] || continue

    bin_size="$(_infer_bin_size "$bn")"

    # Identify file_type (aligned/primary/supplementary/etc.) by removing sample prefix and tokens
    # Keep it simple: classify by common tokens present in filename
    file_type="."
    case "$bn" in
      *primary*) file_type="primary" ;;
      *aligned*) file_type="aligned" ;;
      *supplementary*) file_type="supplementary" ;;
      *discordant*) file_type="discordant" ;;
      *splitters*) file_type="splitters" ;;
    esac

    # track label (basename without .bw and without corr/corrected/expected)
    track="${bn%.bw}"
    track="${track//.corr/}"
    track="${track//.corrected/}"
    track="${track//.expected/}"

    raw_guess="$(_guess_raw_from_corr "$bn")"
    raw_bw=""
    if [[ -n "$raw_guess" && -f "${raw_dir}/${raw_guess}" ]]; then
      raw_bw="${raw_dir}/${raw_guess}"
    else
      raw_guess2="$(_guess_raw_from_corr "$corr_bw")"
      [[ -n "$raw_guess2" && -f "$raw_guess2" ]] && raw_bw="$raw_guess2"
    fi

    if [[ -z "$raw_bw" || ! -f "$raw_bw" ]]; then
      if [[ "$SKIP_MISSING" -eq 1 ]]; then
        continue
      fi
      echo -e "${sample}\t${track}\t${file_type}\t${bin_size}\tFALSE\t.\t${corr_bw}\t0\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tmissing_raw_pair" >> "$OUT"
      continue
    fi

    # ---- row 1: non-dedup stats ----
    stats="$(_union_stats "$raw_bw" "$corr_bw")"
    raw_mean="$(echo "$stats" | awk -F'\t' '{print $2}')"
    raw_cv="$(echo "$stats"   | awk -F'\t' '{print $4}')"
    corr_cv="$(echo "$stats"  | awk -F'\t' '{print $7}')"
    resid_sd="$(echo "$stats" | awk -F'\t' '{print $9}')"

    flag="NA"
    if [[ "$raw_cv" != "NA" && "$corr_cv" != "NA" ]]; then
      flag="$(awk -v a="$corr_cv" -v b="$raw_cv" 'BEGIN{print (a>b)?"TRUE":"FALSE"}')"
    fi

    stability_eff_sd="NA"
    stability_eff_cv="NA"
    if [[ "$bin_size" != "NA" ]]; then
      stability_eff_sd="$(awk -v rs="$resid_sd" -v bs="$bin_size" 'BEGIN{ if(bs>0 && rs!="NA") printf "%.6g", rs/sqrt(bs); else print "NA"}')"
      stability_eff_cv="$(awk -v es="$stability_eff_sd" -v rm="$raw_mean" 'BEGIN{ if(rm>0 && es!="NA") printf "%.6g", es/rm; else print "NA"}')"
    fi

    echo -e "${sample}\t${track}\t${file_type}\t${bin_size}\tFALSE\t${raw_bw}\t${corr_bw}\t${stats}\t${stability_eff_sd}\t${stability_eff_cv}\t${flag}\tok" >> "$OUT"

    # ---- row 2: dedup outputs + stats (default) ----
    if [[ "$DO_DEDUP" -eq 1 ]]; then
      # Find BAM (best effort). We'll try to match "<sample>.<file_type>" first.
      # Example: C01_1.primary.*.bam
      bam_search_dir="${BAM_ROOT:-$sample_dir}"
      track_prefix="${sample}.${file_type}"
      bam="$(_find_bam_for_track "$sample_dir" "$track_prefix" "$bam_search_dir")"

      if [[ -z "$bam" ]]; then
        # fallback: try just sample prefix
        bam="$(_find_bam_for_track "$sample_dir" "$sample" "$bam_search_dir")"
      fi

      if [[ -z "$bam" ]]; then
        echo -e "${sample}\t${track}\t${file_type}\t${bin_size}\tTRUE\t.\t.\t0\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tmissing_bam_for_dedup" >> "$OUT"
        continue
      fi

      # Put dedup bigWigs alongside existing ones
      dedup_raw_bw="${raw_bw%.bw}.dedup.bw"
      dedup_corr_bw="${corr_bw%.bw}.dedup.bw"

      if [[ ! -f "$dedup_raw_bw" || ! -f "$dedup_corr_bw" ]]; then
        if ! _make_dedup_bigwigs "$raw_bw" "$corr_bw" "$bam" "$EPS" "$dedup_raw_bw" "$dedup_corr_bw"; then
          echo -e "${sample}\t${track}\t${file_type}\t${bin_size}\tTRUE\t.\t.\t0\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tdedup_failed" >> "$OUT"
          continue
        fi
      fi

      dedup_stats="$(_union_stats "$dedup_raw_bw" "$dedup_corr_bw")"

      d_raw_mean="$(echo "$dedup_stats" | awk -F'\t' '{print $2}')"
      d_raw_cv="$(echo "$dedup_stats"   | awk -F'\t' '{print $4}')"
      d_corr_cv="$(echo "$dedup_stats"  | awk -F'\t' '{print $7}')"
      d_resid_sd="$(echo "$dedup_stats" | awk -F'\t' '{print $9}')"

      d_flag="NA"
      if [[ "$d_raw_cv" != "NA" && "$d_corr_cv" != "NA" ]]; then
        d_flag="$(awk -v a="$d_corr_cv" -v b="$d_raw_cv" 'BEGIN{print (a>b)?"TRUE":"FALSE"}')"
      fi

      d_stability_eff_sd="NA"
      d_stability_eff_cv="NA"
      if [[ "$bin_size" != "NA" ]]; then
        d_stability_eff_sd="$(awk -v rs="$d_resid_sd" -v bs="$bin_size" 'BEGIN{ if(bs>0 && rs!="NA") printf "%.6g", rs/sqrt(bs); else print "NA"}')"
        d_stability_eff_cv="$(awk -v es="$d_stability_eff_sd" -v rm="$d_raw_mean" 'BEGIN{ if(rm>0 && es!="NA") printf "%.6g", es/rm; else print "NA"}')"
      fi

      echo -e "${sample}\t${track}\t${file_type}\t${bin_size}\tTRUE\t${dedup_raw_bw}\t${dedup_corr_bw}\t${dedup_stats}\t${d_stability_eff_sd}\t${d_stability_eff_cv}\t${d_flag}\tok" >> "$OUT"
    fi
  done
done

echo "[cov_bw_qc] Wrote: $OUT"
