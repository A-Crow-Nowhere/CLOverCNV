#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# cov_model_train.sh
# ============================================================

# -----------------------------
# Logging helpers
# -----------------------------
say() { echo "[cov_model_train] $*"; }
die() { say "ERROR: $*"; exit 1; }
have() { command -v "$1" >/dev/null 2>&1; }

# -----------------------------
# Determine MAPI_ROOT from script location (robust on HPC)
# -----------------------------
if [[ -z "${MAPI_ROOT:-}" ]]; then
  SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
  MAPI_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
fi
export MAPI_ROOT

# -----------------------------
# Defaults
# -----------------------------
R_MODULES="off"   # off|conda|auto|R|"<mods...>"
THREADS=1
KEEP_TMP=0
WINDOW=100
MAPQ=20
DROP_DUP_FOR_TRAIN=0
SEED=42

# R-helper passthrough defaults
MIN_WINDOW_BP=0
EXCLUDE_CONTIGS=()   # repeatable

# -----------------------------
# Usage
# -----------------------------
usage() {
  cat <<EOF
cov_model_train

Required:
  --bam <bam>
  --genome <genome key or fasta path>
  --out-dir <dir>

Optional:
  --window <int>              [${WINDOW}]
  --mapq <int>                [${MAPQ}]
  --threads <int>             [${THREADS}]
  --drop-dup-for-train        Drop BAM duplicates (flag 1024) for depth generation only
  --keep-tmp                  Keep temp dir inside out-dir (debug)

  --min-window-bp <int>       R helper: drop windows smaller than this (default: 0 = off)
  --exclude-contigs <list>    R helper: contigs/chroms to exclude (repeatable)
                              <list> may be comma- or space-separated, e.g.
                                --exclude-contigs "PfDd2_Mt,PfDd2_Api"
                                --exclude-contigs PfDd2_Mt --exclude-contigs PfDd2_Api

  --r-modules <mode>          How to get Rscript on the compute node
                              Modes:
                                off|conda   (default) Do NOT use HPC modules for R.
                                            Use conda R from:
                                              \$MAPI_ROOT/envs/cov_model_train/bin/Rscript
                                            (recommended for non-HPC/local runs)

                                auto|R      Use HPC module system, auto-pick latest R/<ver>.
                                            Then try each prerequisite combo shown by:
                                              module spider R/<ver>
                                            until the module loads.

                                "<mods...>" Explicit modules to load (most reliable on HPC),
                                            e.g.:
                                              --r-modules "gcc/11.4.0 openmpi/4.1.4 R/4.4.1"

  --seed <int>                Hidden/dev: RNG seed for model fitting [${SEED}]
EOF
}

# -----------------------------
# Arg parsing (optparse-free)
# -----------------------------
BAM=""
GENOME=""
OUT_DIR=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --bam) BAM="$2"; shift 2 ;;
    --genome) GENOME="$2"; shift 2 ;;
    --out-dir) OUT_DIR="$2"; shift 2 ;;
    --window) WINDOW="$2"; shift 2 ;;
    --mapq) MAPQ="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --drop-dup-for-train) DROP_DUP_FOR_TRAIN=1; shift ;;
    --keep-tmp) KEEP_TMP=1; shift ;;
    --r-modules) R_MODULES="$2"; shift 2 ;;
    --seed) SEED="$2"; shift 2 ;;

    --min-window-bp) MIN_WINDOW_BP="$2"; shift 2 ;;
    --exclude-contigs) EXCLUDE_CONTIGS+=("$2"); shift 2 ;;

    -h|--help) usage; exit 0 ;;
    *) die "Unknown argument: $1" ;;
  esac
done

[[ -n "$BAM" ]] || die "--bam is required"
[[ -n "$GENOME" ]] || die "--genome is required"
[[ -n "$OUT_DIR" ]] || die "--out-dir is required"
[[ -f "$BAM" ]] || die "BAM not found: $BAM"

# -----------------------------
# Tool checks (non-R)
# -----------------------------
have samtools || die "Missing required tool: samtools"
have bedtools || die "Missing required tool: bedtools"
have mosdepth || die "Missing required tool: mosdepth"
have gzip || die "Missing required tool: gzip"
have zcat || die "Missing required tool: zcat"
have awk || die "Missing required tool: awk"
have sort || die "Missing required tool: sort"

# -----------------------------
# Working directory (create early; FASTA may need staging here)
# -----------------------------
mkdir -p "$OUT_DIR"
WORKDIR="$(mktemp -d "$OUT_DIR/.tmp.cov_model_train.XXXXXX")"
say "Working directory: $WORKDIR"

export TMPDIR="$WORKDIR/Rtmp"
mkdir -p "$TMPDIR"
say "DEBUG: R temp dirs set: TMPDIR=$TMPDIR"

# -----------------------------
# Optional: HPC module loading (smarter)
# -----------------------------
maybe_enable_modules() {
  if command -v module >/dev/null 2>&1; then return 0; fi
  [[ -f /etc/profile.d/modules.sh ]] && source /etc/profile.d/modules.sh && command -v module >/dev/null 2>&1 && return 0
  [[ -f /usr/share/Modules/init/bash ]] && source /usr/share/Modules/init/bash && command -v module >/dev/null 2>&1 && return 0
  return 1
}

# Normalize R_MODULES
R_MODULES="$(echo "${R_MODULES:-}" | awk '{$1=$1;print}')"
case "${R_MODULES,,}" in
  ""|"off"|"conda") R_MODULES="" ;;
  "auto"|"r") R_MODULES="AUTO" ;;
esac

pick_latest_r_module() {
  module spider R 2>&1 \
    | awk '/Versions:/{flag=1;next} flag && $1 ~ /^R\/[0-9]/{print $1}' \
    | sort -V \
    | tail -n 1
}

# Extract all prerequisite "combos" from module spider output.
# Some systems wrap long lines; we join wrapped lines into the prior line.
extract_prereq_combos() {
  local rmod="$1"
  module spider "$rmod" 2>&1 | awk '
    function trim(s){ sub(/^[ \t]+/,"",s); sub(/[ \t]+$/,"",s); return s }
    /You will need to load all module\(s\) on any one of the lines below/ {flag=1; next}
    flag && NF==0 {exit}
    flag {
      line=$0
      # If line is indented and previous exists, treat as continuation
      if (match(line, /^[ \t]+/) && prev != "") {
        prev = prev " " trim(line)
      } else {
        if (prev != "") print trim(prev)
        prev = trim(line)
      }
      next
    }
    END { if (prev != "") print trim(prev) }
  ' | awk '{$1=$1;print}' | sed '/^$/d'
}

try_module_load_silent() {
  local spec="$1"
  # shellcheck disable=SC2086
  module load $spec >/dev/null 2>&1
}

load_r_with_prereqs() {
  local rmod="$1"

  # Try direct
  if try_module_load_silent "$rmod"; then
    return 0
  fi

  # Try each prereq combo line
  mapfile -t combos < <(extract_prereq_combos "$rmod" || true)

  if [[ "${#combos[@]}" -eq 0 ]]; then
    return 1
  fi

  for combo in "${combos[@]}"; do
    say "Trying prerequisites for $rmod: $combo"
    module purge >/dev/null 2>&1 || true

    # load prereqs
    if ! try_module_load_silent "$combo"; then
      continue
    fi

    # now load R
    if try_module_load_silent "$rmod"; then
      return 0
    fi
  done

  return 1
}

if [[ -n "$R_MODULES" ]]; then
  if ! maybe_enable_modules; then
    say "WARN: --r-modules requested, but module system not available; continuing without module loads."
  else
    module purge >/dev/null 2>&1 || true

    if [[ "$R_MODULES" == "AUTO" ]]; then
      R_PICK="$(pick_latest_r_module || true)"
      [[ -z "$R_PICK" ]] && die "Could not find R/<version> via 'module spider R'"

      say "Auto-selected R module: $R_PICK"
      if ! load_r_with_prereqs "$R_PICK"; then
        die "Failed loading $R_PICK (auto). On HPC, explicit is most reliable, e.g. --r-modules \"gcc/... openmpi/... $R_PICK\""
      fi
    else
      say "Loading HPC modules for R: $R_MODULES"
      # shellcheck disable=SC2086
      module load $R_MODULES || die "Failed to load modules: $R_MODULES"
    fi
  fi
fi

# -----------------------------
# Resolve FASTA from --genome
# -----------------------------
say "Resolving FASTA for genome spec '$GENOME' under $MAPI_ROOT/genomes/..."

FASTA=""

# Case 1: user passed a direct path
if [[ -f "$GENOME" ]]; then
  FASTA="$GENOME"
fi

# Case 2: genome key -> search under genomes/
if [[ -z "$FASTA" ]]; then
  FASTA="$(find "$MAPI_ROOT/genomes" -type f \
    \( -iname "${GENOME}.fa" -o -iname "${GENOME}.fasta" -o -iname "${GENOME}.fna" \
       -o -iname "${GENOME}.fa.gz" -o -iname "${GENOME}.fasta.gz" -o -iname "${GENOME}.fna.gz" \) \
    ! -iname "*.fai" ! -iname "*.gzi" \
    2>/dev/null | head -n 1 || true)"
fi

# Case 3: fallback to any fasta if unique
if [[ -z "$FASTA" ]]; then
  mapfile -t FASTAS < <(find "$MAPI_ROOT/genomes" -type f \
    \( -iname "*.fa" -o -iname "*.fasta" -o -iname "*.fna" \
       -o -iname "*.fa.gz" -o -iname "*.fasta.gz" -o -iname "*.fna.gz" \) \
    ! -iname "*.fai" ! -iname "*.gzi" \
    2>/dev/null)

  if [[ "${#FASTAS[@]}" -eq 1 ]]; then
    FASTA="${FASTAS[0]}"
    say "WARN: Only one FASTA found under genomes/; using it: $FASTA"
  fi
fi

[[ -n "$FASTA" ]] || die "Could not resolve FASTA for genome spec '$GENOME' under $MAPI_ROOT/genomes"

# If gzipped, stage uncompressed copy for tools
if [[ "$FASTA" == *.gz ]]; then
  say "Found gzipped FASTA; staging uncompressed copy in WORKDIR..."
  FASTA_GZ="$FASTA"
  FASTA="$WORKDIR/genome.fasta"
  zcat "$FASTA_GZ" > "$FASTA"
fi

say "Resolved FASTA: $FASTA"
[[ -f "$FASTA" ]] || die "Resolved FASTA does not exist: $FASTA"

# -----------------------------
# Resolve Rscript (module -> conda env)
# -----------------------------
resolve_rscript() {
  local cand=""

  cand="$(command -v Rscript 2>/dev/null || true)"
  if [[ -n "$cand" && -x "$cand" ]]; then
    echo "$cand"; return 0
  fi

  if [[ -n "${CONDA_PREFIX:-}" && -x "$CONDA_PREFIX/bin/Rscript" ]]; then
    echo "$CONDA_PREFIX/bin/Rscript"; return 0
  fi

  if [[ -x "$MAPI_ROOT/envs/cov_model_train/bin/Rscript" ]]; then
    echo "$MAPI_ROOT/envs/cov_model_train/bin/Rscript"; return 0
  fi

  return 1
}

RSCRIPT="$(resolve_rscript || true)"
[[ -n "$RSCRIPT" ]] || die "Rscript not found (PATH / CONDA_PREFIX / MAPI env)"
say "DEBUG: Using Rscript: $RSCRIPT"

# -----------------------------
# R execution helpers (NO bare Rscript calls elsewhere)
# -----------------------------
run_r_inline() {
  local expr="$1"
  local prefix="$MAPI_ROOT/envs/cov_model_train"

  if command -v Rscript >/dev/null 2>&1; then
    Rscript --vanilla -e "$expr"
    return $?
  fi

  if command -v conda >/dev/null 2>&1 && [[ -d "$prefix" ]]; then
    conda run -p "$prefix" --no-capture-output Rscript --vanilla -e "$expr"
  else
    "$RSCRIPT" --vanilla -e "$expr"
  fi
}

run_r_file() {
  local script="$1"; shift
  local prefix="$MAPI_ROOT/envs/cov_model_train"

  if command -v Rscript >/dev/null 2>&1; then
    Rscript --vanilla "$script" "$@"
    return $?
  fi

  if command -v conda >/dev/null 2>&1 && [[ -d "$prefix" ]]; then
    conda run -p "$prefix" --no-capture-output Rscript --vanilla "$script" "$@"
  else
    "$RSCRIPT" --vanilla "$script" "$@"
  fi
}

# -----------------------------
# Smoke test + required packages
# -----------------------------
R_SMOKE="$WORKDIR/R_smoke.txt"
if ! run_r_inline 'cat("R_OK\n")' >"$R_SMOKE" 2>&1; then
  say "ERROR: R inline smoke test failed. Output:"
  sed -n '1,200p' "$R_SMOKE" >&2
  die "R cannot run inline code on this node"
fi

R_PKG_CHECK="$WORKDIR/R_pkg_check.txt"
if ! run_r_inline 'pkgs<-c("mgcv","data.table","jsonlite"); miss<-pkgs[!sapply(pkgs, requireNamespace, quietly=TRUE)]; if(length(miss)){cat("MISSING:", paste(miss, collapse=" "), "\n"); quit(status=2)} else {cat("OK\n")}' \
  >"$R_PKG_CHECK" 2>&1; then
  say "ERROR: Required R packages missing for the active R runtime:"
  cat "$R_PKG_CHECK" >&2
  exit 2
fi

# -----------------------------
# Build training windows BED
# -----------------------------
WIN_BED="$WORKDIR/training_windows.bed"
WIN_BED_GZ="$OUT_DIR/training_windows.bed.gz"

FAI="${FASTA}.fai"
if [[ ! -f "$FAI" ]]; then
  say "Indexing FASTA (.fai missing): $FASTA"
  samtools faidx "$FASTA"
fi

say "Building ${WINDOW}bp windows from FASTA index..."
bedtools makewindows -g "$FAI" -w "$WINDOW" > "$WIN_BED"
gzip -c "$WIN_BED" > "$WIN_BED_GZ"

# -----------------------------
# Compute depth per window with mosdepth
# Output required by R helper: chrom start end mean_depth (no header)
# -----------------------------
say "Computing windowed depth with mosdepth..."
MOS_PREFIX="$WORKDIR/mos"

BAM_FOR_DEPTH="$BAM"
if [[ "$DROP_DUP_FOR_TRAIN" -eq 1 ]]; then
  say "Preparing duplicate-dropped BAM for training depth..."
  BAM_FOR_DEPTH="$WORKDIR/train.nodup.bam"
  samtools view -b -F 1024 -@ "$THREADS" "$BAM" > "$BAM_FOR_DEPTH"
  samtools index "$BAM_FOR_DEPTH"
fi

mosdepth --by "$WIN_BED" --mapq "$MAPQ" --threads "$THREADS" --no-per-base "$MOS_PREFIX" "$BAM_FOR_DEPTH"

REG_BED_GZ="${MOS_PREFIX}.regions.bed.gz"
[[ -f "$REG_BED_GZ" ]] || die "mosdepth did not produce regions file: $REG_BED_GZ"

DEPTH_BED="$WORKDIR/depth.mean.bed"
zcat "$REG_BED_GZ" > "$DEPTH_BED"

# -----------------------------
# GC content per window (bedtools nuc)
# -----------------------------
say "Computing GC% per window with bedtools nuc..."
NUC_TSV="$WORKDIR/windows.nuc.tsv"
bedtools nuc -fi "$FASTA" -bed "$WIN_BED" > "$NUC_TSV"

# -----------------------------
# Mappability per window (best-effort)
# -----------------------------
MAP_BED="$WORKDIR/map.mean.bed"
MAP_TRACK=""

for cand in \
  "$MAPI_ROOT/genomes/$GENOME/mappability.${WINDOW}.bed.gz" \
  "$MAPI_ROOT/genomes/$GENOME/mappability.bed.gz" \
  "$MAPI_ROOT/genomes/$GENOME/map.${WINDOW}.bed.gz" \
  "$MAPI_ROOT/genomes/$GENOME/map.bed.gz"
do
  if [[ -f "$cand" ]]; then MAP_TRACK="$cand"; break; fi
done

if [[ -n "$MAP_TRACK" ]]; then
  say "Using mappability track: $MAP_TRACK"
  zcat "$MAP_TRACK" > "$MAP_BED"
else
  say "WARN: No mappability track found for $GENOME; proceeding without mappability (NA)."
  awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"NA"}' "$WIN_BED" > "$MAP_BED"
fi

# -----------------------------
# Locate training helper
# -----------------------------
R_HELPER="$MAPI_ROOT/bin/modules/.CLOverCNV/.cov_model_train_fit.R"
[[ -f "$R_HELPER" ]] || die "Missing R helper: $R_HELPER"

# -----------------------------
# Output paths (canonical)
# -----------------------------
OUT_MODEL="$OUT_DIR/model.rds"
OUT_META="$OUT_DIR/model.meta.json"
OUT_TABLE="$OUT_DIR/training_table.tsv.gz"
OUT_PLOT="$OUT_DIR/diagnostic_fit.png"

# -----------------------------
# Build R helper args for new passthroughs
# -----------------------------
R_HELPER_EXTRA_ARGS=()
if [[ "${MIN_WINDOW_BP:-0}" -gt 0 ]]; then
  R_HELPER_EXTRA_ARGS+=( --min-window-bp "$MIN_WINDOW_BP" )
fi

# pass each exclude-contigs occurrence through (repeatable)
for x in "${EXCLUDE_CONTIGS[@]}"; do
  R_HELPER_EXTRA_ARGS+=( --exclude-contigs "$x" )
done

# -----------------------------
# Run training helper
# -----------------------------
say "Fitting coverage model (R)..."

R_ALL="$WORKDIR/Rscript.all.txt"
if ! run_r_file "$R_HELPER" \
  --depth "$DEPTH_BED" \
  --nuc "$NUC_TSV" \
  --map "$MAP_BED" \
  --window "$WINDOW" \
  --seed "$SEED" \
  --out-model "$OUT_MODEL" \
  --out-meta "$OUT_META" \
  --out-table "$OUT_TABLE" \
  --out-plot "$OUT_PLOT" \
  "${R_HELPER_EXTRA_ARGS[@]}" \
  >"$R_ALL" 2>&1; then
  say "ERROR: R training failed. Output:"
  sed -n '1,200p' "$R_ALL" >&2
  exit 2
fi

say "Wrote: $OUT_MODEL"
say "Wrote: $OUT_META"
say "Wrote: $WIN_BED_GZ"
say "Wrote: $OUT_TABLE"
say "Wrote: $OUT_PLOT"

# -----------------------------
# Cleanup
# -----------------------------
if [[ "$KEEP_TMP" -eq 1 ]]; then
  say "Keeping tmp dir: $WORKDIR"
else
  rm -rf "$WORKDIR"
fi

say "cov_model_train completed successfully"
