#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# blast.sh  (MAPI module)
#   + DB builder integration (helpers in .blast/)
#   + BAM->FASTA input support
#   + SIGPIPE-safe BAM read capping
#
# Helpers expected in:  <this_dir>/.blast/
#   - make_contam_db.sh
#   - summarize_folder.py
#   - summarize_top_genus_proportions.py
# ============================================================

say(){ echo "[blast] $*" >&2; }
die(){ say "ERROR: $*"; exit 1; }

usage(){
  cat <<'EOF'
Usage:
  # Build DB only (no out-dir required):
  mapi modules blast --build-db --db /path/to/db/contam_db

  # Build DB if missing, then blast:
  mapi modules blast --build-db-if-missing --db /path/to/db/contam_db \
    (--fasta-dir DIR | --bam-dir DIR) --out-dir DIR

  # Blast only:
  mapi modules blast (--fasta-dir DIR | --bam-dir DIR) --db PREFIX --out-dir DIR

Inputs:
  --fasta-dir DIR         directory containing *.fasta
  --bam-dir DIR           directory containing *.bam (extract reads -> fasta first)
  --fastas DIR            alias for --fasta-dir
  --bams DIR              alias for --bam-dir

Required:
  --db PREFIX             BLAST db prefix path (e.g. /path/db/contam_db)

Required for BLAST runs:
  --out-dir DIR           output directory

DB build modes:
  --build-db              build DB then exit
  --build-db-if-missing   build DB only if PREFIX.{nin,nhr,nsq} missing, then continue
  --db-out-dir DIR        override where DB is built (default dirname(--db))
  --db-prefix NAME        override DB prefix name (default basename(--db))

BLAST params:
  --threads INT           default 8
  --nseq INT              sequences to BLAST per sample (default 5000; debug uses 100)
  --debug                 blast only 100 sequences
  --min-len INT           filter: minimum sequence length (default 100)
  --max-n-frac FLOAT      filter: max N fraction (default 0.5)
  --qcov INT              default 90
  --evalue STR            default 1e-20
  --dust STR              default no
  --max-target-seqs INT   default 10

Memory guard:
  --min-mem-kb INT        default 700000
  --sleep-after-exit INT  default 30 (retry loop sleep)

Genus summarization:
  --no-entrez             do not query NCBI; use sscinames column from BLAST output
  --entrez-email EMAIL    Entrez email (recommended if not using --no-entrez)

BAM extraction:
  --bam-min-mapq INT      default 0
  --bam-include-unmapped  include unmapped reads (default: mapped only)
  --bam-max-reads INT     cap reads extracted from BAM (0 = no cap)
EOF
}

# -----------------------------
# Defaults
# -----------------------------
FASTA_DIR=""
BAM_DIR=""
DB_PATH=""
OUTDIR=""

SKIP_EXISTING=1

THREADS=8
NSEQ=5000
DEBUG=0

MIN_LEN=100
MAX_N_FRAC="0.5"

QCOV=90
EVALUE="1e-20"
DUST="no"
MAX_TARGET_SEQS=10

MIN_MEM_KB=700000
SLEEP_AFTER_EXIT=30

NO_ENTREZ=0
ENTREZ_EMAIL=""

# BAM opts
BAM_MIN_MAPQ=0
BAM_INCLUDE_UNMAPPED=0
BAM_MAX_READS=0

# DB build opts
BUILD_DB=0
BUILD_DB_IF_MISSING=0
DB_OUT_DIR=""
DB_PREFIX=""

# -----------------------------
# Parse args
# -----------------------------
[[ $# -eq 0 ]] && { usage; exit 0; }
while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help) usage; exit 0 ;;

    --fasta-dir|--fastas) FASTA_DIR="${2:-}"; shift 2 ;;
    --bam-dir|--bams) BAM_DIR="${2:-}"; shift 2 ;;

    --db) DB_PATH="${2:-}"; shift 2 ;;
    --out-dir) OUTDIR="${2:-}"; shift 2 ;;

    --skip-existing) SKIP_EXISTING=1; shift ;;
    --no-skip-existing) SKIP_EXISTING=0; shift ;;

    --threads) THREADS="${2:-}"; shift 2 ;;
    --nseq) NSEQ="${2:-}"; shift 2 ;;
    --debug) DEBUG=1; shift ;;

    --min-len) MIN_LEN="${2:-}"; shift 2 ;;
    --max-n-frac) MAX_N_FRAC="${2:-}"; shift 2 ;;

    --qcov) QCOV="${2:-}"; shift 2 ;;
    --evalue) EVALUE="${2:-}"; shift 2 ;;
    --dust) DUST="${2:-}"; shift 2 ;;
    --max-target-seqs) MAX_TARGET_SEQS="${2:-}"; shift 2 ;;

    --min-mem-kb) MIN_MEM_KB="${2:-}"; shift 2 ;;
    --sleep-after-exit) SLEEP_AFTER_EXIT="${2:-}"; shift 2 ;;

    --no-entrez) NO_ENTREZ=1; shift ;;
    --entrez-email) ENTREZ_EMAIL="${2:-}"; shift 2 ;;

    --bam-min-mapq) BAM_MIN_MAPQ="${2:-}"; shift 2 ;;
    --bam-include-unmapped) BAM_INCLUDE_UNMAPPED=1; shift ;;
    --bam-max-reads) BAM_MAX_READS="${2:-}"; shift 2 ;;

    --build-db) BUILD_DB=1; shift ;;
    --build-db-if-missing) BUILD_DB_IF_MISSING=1; shift ;;
    --db-out-dir) DB_OUT_DIR="${2:-}"; shift 2 ;;
    --db-prefix) DB_PREFIX="${2:-}"; shift 2 ;;

    *) die "Unknown option: $1" ;;
  esac
done

# -----------------------------
# Resolve helper paths
# -----------------------------
MOD_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
HELP_DIR="$MOD_DIR/.blast"

TOPGEN_PY="$HELP_DIR/summarize_top_genus_proportions.py"
FOLDER_SUM_PY="$HELP_DIR/summarize_folder.py"
MAKE_DB_SH="$HELP_DIR/make_contam_db.sh"

if [[ "$BUILD_DB" -eq 1 || "$BUILD_DB_IF_MISSING" -eq 1 ]]; then
  [[ -s "$MAKE_DB_SH" ]] || die "Missing helper: $MAKE_DB_SH"
fi
[[ -s "$TOPGEN_PY" ]] || die "Missing helper: $TOPGEN_PY"
[[ -s "$FOLDER_SUM_PY" ]] || die "Missing helper: $FOLDER_SUM_PY"

# -----------------------------
# Validate core args (build-db gating)
# -----------------------------
[[ -n "$DB_PATH" ]] || die "--db is required"

if [[ "$BUILD_DB" -eq 0 ]]; then
  [[ -n "$OUTDIR" ]] || die "--out-dir is required unless --build-db is used"
fi

if [[ "$BUILD_DB" -eq 0 ]]; then
  if [[ -n "$FASTA_DIR" && -n "$BAM_DIR" ]]; then
    die "Provide only one of --fasta-dir or --bam-dir"
  elif [[ -z "$FASTA_DIR" && -z "$BAM_DIR" ]]; then
    die "Provide one of --fasta-dir or --bam-dir"
  fi
  [[ -z "$FASTA_DIR" || -d "$FASTA_DIR" ]] || die "FASTA dir not found: $FASTA_DIR"
  [[ -z "$BAM_DIR" || -d "$BAM_DIR" ]] || die "BAM dir not found: $BAM_DIR"
else
  if [[ -n "$FASTA_DIR" || -n "$BAM_DIR" || -n "$OUTDIR" ]]; then
    die "--build-db cannot be combined with --fasta-dir/--bam-dir/--out-dir (build DB first, then run blast)"
  fi
fi

# -----------------------------
# Tools
# -----------------------------
command -v blastn >/dev/null 2>&1 || die "blastn not found in PATH"
command -v python3 >/dev/null 2>&1 || die "python3 not found in PATH"
command -v awk >/dev/null 2>&1 || die "awk not found"
command -v sort >/dev/null 2>&1 || die "sort not found"
command -v free >/dev/null 2>&1 || true

if [[ "$BUILD_DB" -eq 0 && -n "$BAM_DIR" ]]; then
  command -v samtools >/dev/null 2>&1 || die "samtools not found (required for --bam-dir)"
fi

if [[ "$BUILD_DB" -eq 1 || "$BUILD_DB_IF_MISSING" -eq 1 ]]; then
  command -v wget >/dev/null 2>&1 || die "wget not found (required for DB build)"
  command -v makeblastdb >/dev/null 2>&1 || die "makeblastdb not found (required for DB build)"
fi

# -----------------------------
# DB build integration
# -----------------------------
db_exists() {
  [[ -s "${DB_PATH}.nhr" && -s "${DB_PATH}.nin" && -s "${DB_PATH}.nsq" ]]
}

run_db_build() {
  local outdir prefix
  outdir="${DB_OUT_DIR:-$(dirname "$DB_PATH")}"
  prefix="${DB_PREFIX:-$(basename "$DB_PATH")}"

  say "Building contaminant DB into: $outdir (prefix=$prefix)"
  mkdir -p "$outdir"
  bash "$MAKE_DB_SH" --out-dir "$outdir" --prefix "$prefix"
  DB_PATH="$outdir/$prefix"
  db_exists || die "DB build finished but DB files missing for prefix: $DB_PATH"
}

if [[ "$BUILD_DB" -eq 1 ]]; then
  run_db_build
  say "DB build complete. Exiting due to --build-db."
  exit 0
fi

if [[ "$BUILD_DB_IF_MISSING" -eq 1 ]]; then
  if db_exists; then
    say "DB exists: $DB_PATH (skipping build)"
  else
    run_db_build
  fi
fi

# -----------------------------
# BLAST runs
# -----------------------------
mkdir -p "$OUTDIR"

# -----------------------------
# Filter FASTA
# -----------------------------
filter_fasta_to_file() {
  local in_fa="$1"
  local out_fa="$2"
  local log="$3"

  awk -v out="$out_fa" -v minlen="$MIN_LEN" -v maxn="$MAX_N_FRAC" '
    BEGIN { RS=">"; ORS="" }
    NR>1 {
      header = substr($0, 1, index($0, "\n")-1)
      seq = substr($0, index($0, "\n")+1)
      gsub(/\n/, "", seq)
      total = length(seq)
      ncount = gsub(/[Nn]/, "", seq)

      if (total >= minlen && (total==0 ? 1 : (ncount/total)) <= maxn) {
        print ">" header "\n" seq "\n" >> out
        passed++
      } else {
        failed++
      }
    }
    END {
      printf("? Filtering complete: %d passed, %d failed\n", passed, failed) > "/dev/stderr"
    }' "$in_fa" 2>>"$log"
}

# -----------------------------
# BAM -> reads FASTA (SIGPIPE-safe)
#   Strategy:
#     1) Convert BAM->FASTA fully (no early exit / no broken pipe)
#     2) If BAM_MAX_READS > 0, cap FASTA to first N records afterward
# -----------------------------
cap_fasta_records_inplace() {
  local in_fa="$1"
  local n="$2"
  local tmp="${in_fa}.cap.tmp"

  awk -v n="$n" '
    BEGIN { RS=">"; ORS=""; c=0 }
    NR==1 { next }
    {
      c++
      if (c<=n) print ">"$0
      else exit
    }
  ' "$in_fa" > "$tmp"
  mv "$tmp" "$in_fa"
}

bam_to_reads_fasta() {
  local bam="$1"
  local out_fa="$2"
  local log="$3"

  local view_flags=()
	if [[ "$BAM_INCLUDE_UNMAPPED" -eq 0 ]]; then
	  view_flags+=( -F 4 )
	  say "BAM extraction: mapped-only (-F 4)" | tee -a "$log"
	else
	  say "BAM extraction: including unmapped reads" | tee -a "$log"
	fi

  if [[ "$BAM_MIN_MAPQ" -gt 0 ]]; then
    view_flags+=( -q "$BAM_MIN_MAPQ" )
  fi

  # write full FASTQ then convert to FASTA without truncating the pipe
  say "BAM->FASTA (no SIGPIPE) for $(basename "$bam")" | tee -a "$log"

  # FASTQ -> FASTA with awk (no seqkit dependency)
  samtools fastq -n "${view_flags[@]}" "$bam" 2>>"$log" \
    | awk '
        NR%4==1 { h=substr($0,2); next }
        NR%4==2 { print ">"h "\n" $0 }
      ' > "$out_fa"

  if [[ "$BAM_MAX_READS" -gt 0 ]]; then
    say "Capping FASTA to first $BAM_MAX_READS records" | tee -a "$log"
    cap_fasta_records_inplace "$out_fa" "$BAM_MAX_READS"
  fi
}

# -----------------------------
# Collect inputs
# -----------------------------
collect_inputs() {
  if [[ -n "$FASTA_DIR" ]]; then
    shopt -s nullglob
    for f in "$FASTA_DIR"/*.fasta; do
      printf "%s\t%s\n" "$(basename "$f" .fasta)" "$f"
    done
    shopt -u nullglob
  else
    shopt -s nullglob
    for b in "$BAM_DIR"/*.bam; do
      printf "%s\t%s\n" "$(basename "$b" .bam)" "$b"
    done
    shopt -u nullglob
  fi
}

# -----------------------------
# Core loop
# -----------------------------
process_samples() {
  local any=0

  while IFS=$'\t' read -r sample inpath; do
    [[ -z "${sample:-}" ]] && continue
    any=1

    local outfile="$OUTDIR/${sample}.nt.blast.out"
    local log="$OUTDIR/${sample}.log"
    local filtered="$OUTDIR/${sample}.filtered.fasta"

    if [[ "$SKIP_EXISTING" -eq 1 && -s "$outfile" ]]; then
      echo "$sample already completed, skipping." | tee -a "$log"
      continue
    fi

    echo "$(date) - Starting $sample" | tee "$log"

    local input_fasta=""
    if [[ -n "$FASTA_DIR" ]]; then
      input_fasta="$inpath"
    else
      local reads_fa="$OUTDIR/${sample}.bam_reads.fasta"
      echo "$(date) - Extracting reads FASTA from BAM: $inpath" | tee -a "$log"
      bam_to_reads_fasta "$inpath" "$reads_fa" "$log"
      [[ -s "$reads_fa" ]] || { echo "?? BAM->FASTA empty; skipping." | tee -a "$log"; continue; }
      input_fasta="$reads_fa"
    fi

    echo "$(date) - Filtering $sample" | tee -a "$log"
    rm -f "$filtered"
    filter_fasta_to_file "$input_fasta" "$filtered" "$log"

    if [[ ! -s "$filtered" ]]; then
      echo "?? No reads passed filter for $sample. Skipping." | tee -a "$log"
      continue
    fi

    local total_seqs
    total_seqs=$(grep -c '^>' "$filtered" || true)

    local nblast
    if [[ "$DEBUG" -eq 1 ]]; then
      nblast=100
    else
      nblast="$NSEQ"
    fi

    local blast_input="$filtered"
    if [[ "$nblast" -lt "$total_seqs" ]]; then
      blast_input="$OUTDIR/${sample}.blast_input.fasta"
      awk -v n="$nblast" 'BEGIN{RS=">"; ORS=""} NR==1{next} NR<=n+1{print ">"$0}' \
        "$filtered" > "$blast_input"
      echo "Blasting only first $nblast sequences for $sample" | tee -a "$log"
    fi

    local free_kb
    free_kb=$(awk '/MemAvailable/ {print $2}' /proc/meminfo)
    if [[ "${free_kb:-999999999}" -lt "$MIN_MEM_KB" ]]; then
      echo "?? Low available memory (${free_kb} KB) before BLAST. Exiting early." | tee -a "$log"
      rm -f "$blast_input"
      return 1
    fi

    echo "$(date) - BLASTing $sample" | tee -a "$log"
    blastn -query "$blast_input" \
      -db "$DB_PATH" \
      -outfmt '6 qseqid sseqid pident length mismatch evalue stitle' \
      -max_target_seqs "$MAX_TARGET_SEQS" \
      -num_threads "$THREADS" \
      -dust "$DUST" \
      -qcov_hsp_perc "$QCOV" \
      -evalue "$EVALUE" \
      -out "$outfile.tmp" >>"$log" 2>&1

    local rc=$?
    if [[ $rc -ne 0 ]]; then
      echo "? BLAST failed (exit $rc) for $sample. Skipping." | tee -a "$log"
      rm -f "$outfile.tmp" "$blast_input"
      continue
    fi

    sort -t $'\t' -k1,1 -k6,6g "$outfile.tmp" | awk -F $'\t' '!seen[$1]++' > "$outfile"
    rm -f "$outfile.tmp" "$blast_input"

    echo "$(date) - ? Finished BLAST for $sample" | tee -a "$log"

    if [[ "$NO_ENTREZ" -eq 1 ]]; then
      python3 "$TOPGEN_PY" --no-entrez "$outfile" | tee -a "$log"
    else
      if [[ -n "$ENTREZ_EMAIL" ]]; then
        python3 "$TOPGEN_PY" --entrez-email "$ENTREZ_EMAIL" "$outfile" | tee -a "$log"
      else
        python3 "$TOPGEN_PY" "$outfile" | tee -a "$log"
      fi
    fi

    free_kb=$(awk '/MemAvailable/ {print $2}' /proc/meminfo)
    if [[ "${free_kb:-999999999}" -lt "$MIN_MEM_KB" ]]; then
      echo "?? Low available memory (${free_kb} KB) after $sample. Exiting early." | tee -a "$log"
      return 1
    fi

    sleep 2
  done < <(collect_inputs)

  [[ $any -eq 1 ]] || say "No inputs found."
  return 0
}

say "Starting BLAST loop. debug=$DEBUG"
while true; do
  process_samples
  rc=$?
  if [[ $rc -eq 0 ]]; then
    break
  fi
  say "Low-mem early exit; sleeping $SLEEP_AFTER_EXIT seconds then retrying..."
  sleep "$SLEEP_AFTER_EXIT"
done

say "Writing folder summary"
if ls "$OUTDIR"/*.genus_proportions.tsv >/dev/null 2>&1; then
  python3 "$FOLDER_SUM_PY" --blast-dir "$OUTDIR" --out "$OUTDIR/summary_total.tsv"
else
  say "No genus_proportions.tsv produced; writing empty summary_total.tsv"
  printf "Category\tMean%%\tStdDev%%\tTotal\tN_samples\n" > "$OUTDIR/summary_total.tsv"
fi
say "Done."

