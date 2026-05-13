#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# CLOverCNV.sh
#   Wrapper pipeline:
#     Model_Train -> Model_Apply -> Model_Finalize (bw QC)
#     CNV_Probe -> CNV_Segment -> CNV_Boundary -> CNV_Confidence
#     -> CNV_Finalize -> CNV_Finalize2
#
#   Supports:
#     mapi modules CLOverCNV step <AliasOrInternalOrShort> [options...]
#     mapi modules CLOverCNV run  [pipeline options...]
#
# Notes:
#   - Python-based CNV steps run through the CLOverCNV conda env when present.
#   - R-based segmentation can use either conda R or HPC module R.
#   - Wrapper defaults may intentionally differ from downstream script defaults.
#   - CNV_Finalize2 is a downstream post-pass only. It reads the output of the
#     existing finalize step and does not modify upstream segmentation,
#     boundary-support, or confidence logic.
# ============================================================

# -----------------------------
# Logging / errors
# -----------------------------
LOG_TS_FMT='+%Y-%m-%d %H:%M:%S'
msg() { echo "[CLOverCNV] $(date "$LOG_TS_FMT") $*" >&2; }
die() { msg "ERROR: $*"; exit 1; }

# -----------------------------
# Determine MAPI_ROOT (robust)
# -----------------------------
THIS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MAPI_ROOT="${MAPI_ROOT:-"$(cd "$THIS_DIR/../.." && pwd)"}"
export MAPI_ROOT

msg "DEBUG: wrapper path = ${BASH_SOURCE[0]}"
msg "DEBUG: MAPI_ROOT = $MAPI_ROOT"

# Hidden tool dir (modules live here, wrapper stays visible)
HIDDEN_DIR="$MAPI_ROOT/bin/modules/.CLOverCNV"
[[ -d "$HIDDEN_DIR" ]] || die "Missing hidden dir: $HIDDEN_DIR"

# Conda for Python-side execution
CONDA="$MAPI_ROOT/tools/miniconda3/bin/conda"
[[ -x "$CONDA" ]] || die "conda not found/executable at: $CONDA"

# CLOverCNV env path
CNV_ENV="$MAPI_ROOT/envs/CLOverCNV"
[[ -d "$CNV_ENV" ]] || msg "WARN: env dir not found (will still try PATH tools): $CNV_ENV"

# -----------------------------
# Defaults
# -----------------------------
SAMPLE=""
GENOME=""
GFF=""
OUT_ROOT=""
BAM=""
BAM_DIR=""
MODEL_IN=""

# Model train/apply defaults
TRAIN_WINDOW=100
TRAIN_MAPQ=20
TRAIN_THREADS=8
TRAIN_DROP_DUP=1
TRAIN_EXCLUDE_CONTIGS=()

APPLY_BINSIZE=50
APPLY_THREADS=8
APPLY_EMIT_RAW=0
APPLY_PRIMARY_PREFER="primary"   # wrapper-supported values: primary|aligned

# R runtime selection
R_MODULES="off"                  # off|none|auto|conda|"<mods...>"
R_AUTOINSTALL="${R_AUTOINSTALL:-1}"
MAPI_R_LIBS_USER_ROOT="${MAPI_R_LIBS_USER_ROOT:-}"

# CNV defaults (wrapper defaults; may differ from downstream script defaults)
TILE_BP=75
MIN_TILE_BP=50
EXCLUDE_BED=""
LAMBDA=25
MIN_PROBES_PER_SEG=10
FLANK=1250
AGG="sum"                        # mean|sum|max
CONF_METHOD="mean"               # min|mean|max

FINAL_WEAK_RATIO=2.0
FINAL_STRONG_RATIO=2.3
FINAL_WEAK_Z=1.5
FINAL_STRONG_Z=2.0
FINAL_KEEP_MODE="both"           # ratio_only|z_only|both
CN_BASE=1.0

FINAL_FUSE=0
FINAL_FUSE_MAX_GAP=0

FINAL_COUNT_FLANK=1250
FINAL_COUNT_ALL_ALIGN=0
FINAL_COUNT_DROP_DUP=0
FINAL_COUNT_INCLUDE_SUPP=0
FINAL_COUNT_INCLUDE_SECONDARY=0

FINAL_DEBUG=0

# Finalize2 defaults
FINAL2_BALANCE=0.33
FINAL2_RESCUE_Z=10
FINAL2_KEEP_SUGGESTED_ONLY=0
FINAL2_KEEP_WORKING_DIR=0
FINAL2_VERBOSE=0

FORCE=0
KEEP_TMP=0

CIRCULAR_CONTIGS=""
CIRCULAR_CONTIGS_FILE=""

START_AT="cov_model_train"
STOP_AFTER="cnv_finalize_suggested"

# -----------------------------
# Step registry + aliases
# -----------------------------
declare -A STEP_ALIAS=(
  ["Model_Train"]="cov_model_train"
  ["Model_Apply"]="cov_model_apply"
  ["Model_Finalize"]="cov_bw_qc"
  ["CNV_Probe"]="make_probes"
  ["CNV_Segment"]="fused_lasso_segment"
  ["CNV_Boundary"]="boundary_support_from_bigwigs"
  ["CNV_Confidence"]="segment_confidence_from_boundaries"
  ["CNV_Finalize"]="cnv_finalize_segments"
  ["CNV_Finalize2"]="cnv_finalize_suggested"

  # short ergonomic aliases
  ["train"]="cov_model_train"
  ["apply"]="cov_model_apply"
  ["bw_qc"]="cov_bw_qc"
  ["qc"]="cov_bw_qc"
  ["probe"]="make_probes"
  ["segment"]="fused_lasso_segment"
  ["boundary"]="boundary_support_from_bigwigs"
  ["confidence"]="segment_confidence_from_boundaries"
  ["finalize"]="cnv_finalize_segments"
  ["finalize2"]="cnv_finalize_suggested"
  ["suggested"]="cnv_finalize_suggested"
)

steps=(
  cov_model_train
  cov_model_apply
  cov_bw_qc
  make_probes
  fused_lasso_segment
  boundary_support_from_bigwigs
  segment_confidence_from_boundaries
  cnv_finalize_segments
  cnv_finalize_suggested
)

resolve_step() {
  local s="$1"
  if [[ -n "${STEP_ALIAS[$s]:-}" ]]; then
    echo "${STEP_ALIAS[$s]}"
    return 0
  fi
  local x
  for x in "${steps[@]}"; do
    [[ "$x" == "$s" ]] && { echo "$s"; return 0; }
  done
  echo ""
}

idx() {
  local needle="$1"
  local i
  for i in "${!steps[@]}"; do
    [[ "${steps[$i]}" == "$needle" ]] && { echo "$i"; return 0; }
  done
  echo "-1"
}

# -----------------------------
# Usage
# -----------------------------
usage() {
cat <<'EOF'

CLOverCNV - coverage-normalized CDS-based CNV calling pipeline

USAGE
  CLOverCNV run [options]
  CLOverCNV step <step_name> [options]

STEP NAMES
  Long aliases:
    Model_Train
    Model_Apply
    Model_Finalize
    CNV_Probe
    CNV_Segment
    CNV_Boundary
    CNV_Confidence
    CNV_Finalize
    CNV_Finalize2

  Short aliases:
    train | apply | qc | probe | segment | boundary | confidence | finalize | finalize2 | suggested

  Internal names:
    cov_model_train
    cov_model_apply
    cov_bw_qc
    make_probes
    fused_lasso_segment
    boundary_support_from_bigwigs
    segment_confidence_from_boundaries
    cnv_finalize_segments
    cnv_finalize_suggested


GENERAL INPUTS
  --sample STR                Sample name / identifier. Required for run mode and
                              strongly recommended for all step calls. Finalize2 will
                              also write/populate sample_name from this value.

  --genome STR                Genome key or FASTA path (required for run mode;
                              probe step also requires it explicitly)

  --gff FILE                  Genome annotation in GFF/GFF3 format (.gz ok)

  --bam FILE                  Primary BAM for the sample

  --bam-dir DIR               Directory containing BAMs for derivative/support steps

  --model FILE                Pretrained coverage model to reuse instead of training


OUTPUT / PIPELINE CONTROL
  --out-root DIR              Root output directory for all CLOverCNV results

  --start-at STEP             Start at a specific pipeline step
                              Default: cov_model_train

  --stop-after STEP           Stop after a specific pipeline step
                              Default: cnv_finalize_suggested

  --force                     Overwrite / rerun steps even if outputs exist

  --keep-tmp                  Keep temporary/intermediate files created by upstream
                              steps such as model application


---------------------------------------------------------
MODEL TRAINING
---------------------------------------------------------

  --train-window INT          Window size (bp) used to train the coverage correction model.
                              Default: 100

  --train-mapq INT            Minimum MAPQ for reads used in model training.
                              Default: 20

  --train-threads INT         Number of threads used for model training.
                              Default: 8

  --train-drop-dup            Exclude duplicate-marked reads during training.
                              Wrapper default behavior.

  --train-no-drop-dup         Keep duplicate-marked reads during training.

  --train-exclude-contigs STR Exclude a contig from model training.
                              May be supplied multiple times.


---------------------------------------------------------
MODEL APPLICATION
---------------------------------------------------------

  --apply-binsize INT         Bin size (bp) for corrected coverage output.
                              Default: 50

  --apply-threads INT         Number of threads used during model application.
                              Default: 8

  --apply-emit-raw            Also emit raw (uncorrected) coverage tracks.

  --apply-primary-prefer STR  Primary/alignment preference passed to apply step.
                              Wrapper default: primary
                              Intended values: primary | aligned


---------------------------------------------------------
R / HPC ENVIRONMENT
---------------------------------------------------------

  --r-modules STR             HPC module string to load before R-based steps.
                              Typical values:
                                off
                                none
                                auto
                                conda
                                'gcc/11.4.0 openmpi/4.1.4 R/4.3.1'


---------------------------------------------------------
STEP: PROBE GENERATION
(CDS tiling + corrected coverage summarization)
---------------------------------------------------------

This step tiles CDS space into fixed-size probes and computes overlap-weighted mean
coverage for each probe from the corrected primary coverage bigWig.

  --tile-bp INT               Probe tile size in bp within CDS space.
                              Wrapper default: 75

  --min-tile-bp INT           Minimum retained tile size after edge clipping.
                              Wrapper default: 50

  --exclude-bed FILE          Optional BED of regions to exclude from probe generation.

  Notes:
    - Probe values are computed from the corrected primary bigWig.
    - Tiles are built in CDS space, not across all genomic bases.
    - The downstream make_probes.py also has an --eps pseudocount argument,
      but the wrapper currently uses the script default and does not expose it.


---------------------------------------------------------
STEP: SEGMENTATION
(PELT changepoint calling on probe-level signal)
---------------------------------------------------------

This step segments the probe table into regions of approximately constant mean signal.

  --lambda FLOAT              Segmentation penalty parameter.
                              Higher values => FEWER segments
                              Lower values  => MORE segments
                              Wrapper default: 25

  --min-probes-per-seg INT    Minimum probes allowed per segment.
                              Wrapper default: 10

  Notes:
    - CNV boundaries are the transitions between adjacent segments.
    - Segment genomic spans may appear broad if CDS probes are sparse.


---------------------------------------------------------
STEP: BOUNDARY SUPPORT
(split/disco support around breakpoints)
---------------------------------------------------------

This step quantifies breakpoint support in flanking windows around segment boundaries.

  --flank INT                 Boundary flank size in bp on each side.
                              Wrapper default: 1250

  --agg STR                   Aggregation method for support signal in the flank window.
                              Allowed values:
                                mean
                                sum
                                max
                              Wrapper default: sum

                              mean  = average support across window
                              sum   = total support across window
                              max   = peak support within window


---------------------------------------------------------
STEP: CONFIDENCE SCORING
(combine left/right boundary support into segment confidence)
---------------------------------------------------------

This step combines left/right boundary evidence into a per-segment confidence score.

  --confidence-method STR     Combination rule for left/right boundary support.
                              Allowed values:
                                min
                                mean
                                max
                              Wrapper default: mean

                              min   = weaker boundary determines confidence
                              mean  = average of both boundaries
                              max   = stronger boundary determines confidence

  Notes:
    - Wrapper maps this to downstream: --method
    - Wrapper currently fixes downstream --score-from z
    - Additional downstream options such as --missing-policy, --z-col,
      --raw-col, --fallback-pos-match, and --tol-bp are not exposed here


---------------------------------------------------------
STEP: FINAL CNV CLASSIFICATION / FILTERING
(effect size + support tiering, optional fusion and BAM counting)
---------------------------------------------------------

This step combines:
  1. segment ratio / effect size
  2. boundary z support
to classify and optionally fuse CNV calls.

RATIO THRESHOLDS
  --final-weak-ratio FLOAT    Weak CNV ratio threshold.
                              Must be > 1.0
                              Wrapper default: 2.0

  --final-strong-ratio FLOAT  Strong CNV ratio threshold.
                              Must be >= weak-ratio
                              Wrapper default: 2.3

Z THRESHOLDS
  --final-weak-z FLOAT        Weak boundary support threshold.
                              Must be >= 0
                              Wrapper default: 1.5

  --final-strong-z FLOAT      Strong boundary support threshold.
                              Must be >= weak-z
                              Wrapper default: 2.0

KEEP MODE
  --final-keep-mode STR       Rule for keeping CNV calls.
                              Allowed values:
                                ratio_only
                                z_only
                                both
                              Wrapper default: both

                              ratio_only = keep based only on ratio effect size
                              z_only     = keep based only on boundary z support
                              both       = require both ratio and z support

CN BASELINE
  --cn-base FLOAT             Neutral copy-number baseline used to interpret ratios.
                              Wrapper default: 1.0

OPTIONAL FUSION
  --final-fuse                Fuse adjacent retained CNV segments of same direction
                              (gain+gain or loss+loss)

  --final-fuse-max-gap INT    Maximum bp gap allowed between fused segments
                              Wrapper default: 0

OPTIONAL BAM-BASED READ COUNTING
  --final-count-flank INT     Flank size (bp) around boundaries for counting
                              Wrapper default: 1250

  --final-count-all-alignments
                              Count all alignments

  --final-count-drop-dup      Exclude duplicate-marked reads from counts

  --final-count-include-supp  Include supplementary alignments in counts

  --final-count-include-secondary
                              Include secondary alignments in counts

CIRCULAR CONTIG HANDLING
  --circular-contigs STR      Comma-separated contig names to treat as circular

  --circular-contigs-file FILE
                              File with one circular contig name per line

DEBUGGING
  --debug                     Pass verbose debug mode to finalize step


---------------------------------------------------------
STEP: FINALIZE2 / SUGGESTED CALL POST-PASS
(balance, rescue, strength override, and shared-boundary adjudication)
---------------------------------------------------------

This is a downstream last-pass adjudication layer that reads the output of the
existing finalize step and updates/augments the suggested keep/call fields.

Conceptually, it aims to reduce false positive segments that only inherit one
strong shared boundary from a neighboring true CNV, while preserving very strong
high-amplitude CNVs that may be more asymmetric.

General behavior:
  - existing weak/strong z thresholds are reused from the finalize step
  - a single main balance knob controls how balanced the two sides must be
  - one rescue threshold allows extreme one-sided support to survive when the
    opposite boundary is below weak_z
  - rescue bypasses balance
  - a failed balance test can still be overridden for strong CNVs
  - neighboring segments may be compared to identify shared-boundary passengers

Terminal / NA handling:
  - first/last chromosome segments often have NA left/right boundary z because
    they do not have two ordinary interior breakpoints
  - finalize2 falls back to counted support columns for these cases
    when available, especially:
      left_split_reads + left_disco_reads
      right_split_reads + right_disco_reads
      left_support_reads / right_support_reads

Required upstream input:
  - output of CNV_Finalize, typically segments.calls.tsv.gz

Main post-pass knobs:
  --final2-balance FLOAT      Single main balance heuristic knob.
                              Expected interpretation:
                                min(sideA, sideB) / max(sideA, sideB)
                              Higher => stricter balance requirement
                              Lower  => more permissive
                              Default: 0.33

  --final2-rescue-z FLOAT     Extreme one-sided z threshold that can rescue a
                              segment when the opposite side is below weak_z.
                              Default: 10

Output / filtering:
  --keep-suggested-only       After finalize2, emit only rows with updated
                              keep_suggested == TRUE

  --keep-working-dir          Preserve the finalize2 working directory instead
                              of cleaning/staging only final copied outputs

  --verbose                   Keep all diagnostic/helper columns produced by
                              finalize2. Default output is trimmed to a cleaner
                              set of columns.

Notes:
  - finalize2 is intended to update/add fields such as:
      sample_name
      keep_suggested
      call
      balance_ratio
      keep_ratio / keep_z / keep_balance / keep_neighbor
      imbalance_override
      heuristic_reason
  - finalize2 should not alter upstream probe, segment, boundary, or confidence
    calculations


---------------------------------------------------------
EXAMPLES
---------------------------------------------------------

  Full run through finalize2:
    CLOverCNV run \
      --sample C01_1 \
      --genome Clo \
      --gff genome.gff \
      --bam sample.primary.bam \
      --bam-dir bam_dir \
      --out-root results \
      --tile-bp 75 \
      --min-tile-bp 75 \
      --lambda 25 \
      --min-probes-per-seg 5 \
      --flank 1000 \
      --agg mean \
      --confidence-method mean \
      --final-weak-ratio 2.0 \
      --final-strong-ratio 2.3 \
      --final-weak-z 1.5 \
      --final-strong-z 2.0 \
      --final-keep-mode both \
      --final2-balance 0.33 \
      --final2-rescue-z 10

  Run only the post-pass suggested-call refinement:
    CLOverCNV run \
      --sample C01_1 \
      --genome Clo \
      --out-root results \
      --start-at finalize2 \
      --stop-after finalize2 \
      --final2-balance 0.25 \
      --final2-rescue-z 10 \
      --keep-suggested-only

  Run finalize2 and retain all diagnostic columns:
    CLOverCNV step finalize2 \
      --sample C01_1 \
      --out-root results \
      --final2-balance 0.33 \
      --final2-rescue-z 10 \
      --verbose

  Single step:
    CLOverCNV step finalize2 \
      --sample C01_1 \
      --out-root results \
      --final2-balance 0.33 \
      --final2-rescue-z 10 \
      --keep-suggested-only


HELP
  -h, --help                  Show this help message

EOF
}

# -----------------------------
# Value validators
# -----------------------------
validate_choice() {
  local value="$1"; shift
  local label="$1"; shift
  local ok=0
  local x
  for x in "$@"; do
    if [[ "$value" == "$x" ]]; then
      ok=1
      break
    fi
  done
  [[ "$ok" -eq 1 ]] || die "Invalid value for $label: '$value' (allowed: $*)"
}

validate_args() {
  validate_choice "$AGG" "--agg" mean sum max
  validate_choice "$CONF_METHOD" "--confidence-method" min mean max
  validate_choice "$FINAL_KEEP_MODE" "--final-keep-mode" ratio_only z_only both
  validate_choice "$APPLY_PRIMARY_PREFER" "--apply-primary-prefer" primary aligned

  awk -v x="$FINAL_WEAK_RATIO" 'BEGIN{exit !(x>1.0)}' || die "--final-weak-ratio must be > 1.0"
  awk -v x="$FINAL_STRONG_RATIO" -v y="$FINAL_WEAK_RATIO" 'BEGIN{exit !(x>=y)}' || die "--final-strong-ratio must be >= --final-weak-ratio"
  awk -v x="$FINAL_WEAK_Z" 'BEGIN{exit !(x>=0)}' || die "--final-weak-z must be >= 0"
  awk -v x="$FINAL_STRONG_Z" -v y="$FINAL_WEAK_Z" 'BEGIN{exit !(x>=y)}' || die "--final-strong-z must be >= --final-weak-z"
  awk -v x="$FINAL2_BALANCE" 'BEGIN{exit !(x>=0 && x<=1)}' || die "--final2-balance must be between 0 and 1"
  awk -v x="$FINAL2_RESCUE_Z" 'BEGIN{exit !(x>=0)}' || die "--final2-rescue-z must be >= 0"
}

# -----------------------------
# R user-library helpers (HPC-safe)
# -----------------------------
pick_r_user_lib_root() {
  if [[ -n "${MAPI_R_LIBS_USER_ROOT:-}" ]]; then
    echo "$MAPI_R_LIBS_USER_ROOT"
    return 0
  fi

  if [[ -d "/standard" && -n "${HPC_ALLOCATION:-}" ]]; then
    local cand="/standard/${HPC_ALLOCATION}/${USER}/R/goolf"
    if mkdir -p "$cand" >/dev/null 2>&1; then
      echo "$cand"
      return 0
    fi
  fi

  echo "$HOME/R/goolf"
}

ensure_r_userlib_and_pkgs() {
  command -v Rscript >/dev/null 2>&1 || die "Rscript not found (cannot configure R user library)"

  local r_mm root
  r_mm="$(Rscript --vanilla -e 'cat(paste0(R.version$major,".",strsplit(R.version$minor,"\\.")[[1]][1]))')"
  root="$(pick_r_user_lib_root)"

  export R_LIBS_USER="${root}/${r_mm}"
  mkdir -p "$R_LIBS_USER" || die "Failed to create R_LIBS_USER: $R_LIBS_USER"
  msg "DEBUG: R_LIBS_USER=$R_LIBS_USER"

  if [[ "${R_AUTOINSTALL:-1}" -eq 1 ]]; then
    if ! Rscript --vanilla -e '
req <- c("data.table","R.utils","changepoint")
miss <- req[!vapply(req, requireNamespace, logical(1), quietly=TRUE)]
if (length(miss)) {
  message("[R preflight] installing: ", paste(miss, collapse=", "))
  install.packages(miss, repos="https://cloud.r-project.org")
} else {
  message("[R preflight] ok")
}
' >/dev/null 2>&1; then
      return 1
    fi
  fi

  Rscript --vanilla -e '
req <- c("data.table","R.utils","changepoint")
miss <- req[!vapply(req, requireNamespace, logical(1), quietly=TRUE)]
if (length(miss)) {
  cat("MISSING:", paste(miss, collapse=" "), "\n")
  print(.libPaths())
  quit(status=2)
}
' >/dev/null 2>&1 || return 1

  if [[ "${MAPI_DEBUG:-0}" == "1" ]]; then
    Rscript --vanilla -e 'cat("DEBUG .libPaths():\n"); print(.libPaths())' >&2 || true
  fi
}

# -----------------------------
# Module system helpers (for R)
# -----------------------------
maybe_enable_modules() {
  if command -v module >/dev/null 2>&1; then return 0; fi
  [[ -f /etc/profile.d/modules.sh ]] && source /etc/profile.d/modules.sh && command -v module >/dev/null 2>&1 && return 0
  [[ -f /usr/share/Modules/init/bash ]] && source /usr/share/Modules/init/bash && command -v module >/dev/null 2>&1 && return 0
  return 1
}

pick_latest_r_module() {
  module spider R 2>&1 \
    | awk '/Versions:/{flag=1;next} flag && $1 ~ /^R\/[0-9]/{print $1}' \
    | sort -V \
    | tail -n 1
}

extract_prereq_combos() {
  local rmod="$1"
  module spider "$rmod" 2>&1 | awk '
    function trim(s){ sub(/^[ \t]+/,"",s); sub(/[ \t]+$/,"",s); return s }
    /You will need to load all module\(s\) on any one of the lines below/ {flag=1; next}
    flag && NF==0 {exit}
    flag {
      line=$0
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
  if try_module_load_silent "$rmod"; then return 0; fi
  mapfile -t combos < <(extract_prereq_combos "$rmod" || true)
  [[ "${#combos[@]}" -gt 0 ]] || return 1
  local combo
  for combo in "${combos[@]}"; do
    msg "Trying prerequisites for $rmod: $combo"
    module purge >/dev/null 2>&1 || true
    if ! try_module_load_silent "$combo"; then
      continue
    fi
    if try_module_load_silent "$rmod"; then
      return 0
    fi
  done
  return 1
}

run_r_file_hpc() {
  local script="$1"; shift
  [[ -f "$script" ]] || die "R script not found: $script"

  export TMPDIR="$WORK/tmp_R"
  mkdir -p "$TMPDIR"

  local rm
  rm="$(echo "${R_MODULES:-off}" | awk '{$1=$1;print}')"
  case "${rm,,}" in
    ""|"off"|"none") rm="" ;;
    "auto"|"r") rm="AUTO" ;;
    "conda") rm="" ;;
    *) : ;;
  esac

  if [[ -n "$rm" ]]; then
    maybe_enable_modules || die "R modules requested but module system not available on this node"
    module purge >/dev/null 2>&1 || true

    if [[ "$rm" == "AUTO" ]]; then
      local rpick
      rpick="$(pick_latest_r_module || true)"
      [[ -n "$rpick" ]] || die "Could not find R/<version> via 'module spider R'"
      msg "Auto-selected R module: $rpick"
      load_r_with_prereqs "$rpick" || die "Failed loading $rpick (auto). Use explicit --r-modules \"gcc/... openmpi/... $rpick\""
    else
      msg "Loading HPC modules for R: $rm"
      # shellcheck disable=SC2086
      module load $rm || die "Failed to load modules: $rm"
    fi
  fi

  command -v Rscript >/dev/null 2>&1 || die "Rscript not found in PATH (after module handling)"
  msg "DEBUG: Rscript=$(command -v Rscript)"
  msg "DEBUG: TMPDIR=$TMPDIR"

  if ! ensure_r_userlib_and_pkgs; then
    Rscript --vanilla -e 'cat("R pkg preflight failed.\n"); print(.libPaths())' >&2 || true
    die "R runtime missing required packages (data.table/R.utils/changepoint) and auto-install failed. If compute nodes lack CRAN access, pre-install once on a login node into: $R_LIBS_USER"
  fi

  Rscript --vanilla "$script" "$@"
}

# -----------------------------
# Execution helper
# -----------------------------
run_step() {
  local name="$1"; shift
  local logfile="$1"; shift
  mkdir -p "$(dirname "$logfile")"
  msg ">>> $name"
  { "$@"; } 2>&1 | tee "$logfile"
}

# -----------------------------
# Small utilities
# -----------------------------
ensure_env_python() {
  if [[ -d "$CNV_ENV" ]]; then
    echo "$CONDA run -p $CNV_ENV --no-capture-output"
  else
    echo ""
  fi
}

pick_primary_corr_bw() {
  local corr_dir="$1"
  local sample="$2"
  local bins="$3"

  local cand=""
  cand="$(ls -1 "$corr_dir/${sample}.primary.corr.bins${bins}.bw" 2>/dev/null | head -n 1 || true)"
  [[ -n "$cand" ]] && { echo "$cand"; return 0; }
  cand="$(ls -1 "$corr_dir/${sample}.aligned.corr.bins${bins}.bw" 2>/dev/null | head -n 1 || true)"
  [[ -n "$cand" ]] && { echo "$cand"; return 0; }

  cand="$(ls -1 "$corr_dir"/*.bw 2>/dev/null | head -n 1 || true)"
  echo "$cand"
}

pick_primary_bam() {
  local bam_dir="$1"
  local sample="$2"

  local cand=""
  cand="$(ls -1 "$bam_dir/${sample}."*primary*.bam 2>/dev/null | head -n 1 || true)"
  [[ -n "$cand" ]] && { echo "$cand"; return 0; }
  cand="$(ls -1 "$bam_dir/${sample}."*aligned*.bam 2>/dev/null | head -n 1 || true)"
  [[ -n "$cand" ]] && { echo "$cand"; return 0; }

  cand="$(ls -1 "$bam_dir/${sample}."*.bam 2>/dev/null | grep -Ev 'split|disco|discord|supp|unmap|unmapped' | head -n 1 || true)"
  [[ -n "$cand" ]] && { echo "$cand"; return 0; }

  cand="$(ls -1 "$bam_dir"/*.bam 2>/dev/null | head -n 1 || true)"
  echo "$cand"
}

pick_split_bam() {
  local bam_dir="$1"
  local sample="$2"
  local cand=""
  cand="$(ls -1 "$bam_dir/${sample}."*split*.bam 2>/dev/null | head -n 1 || true)"
  [[ -n "$cand" ]] && { echo "$cand"; return 0; }
  cand="$(ls -1 "$bam_dir"/*split*.bam 2>/dev/null | head -n 1 || true)"
  echo "$cand"
}

pick_disco_bam() {
  local bam_dir="$1"
  local sample="$2"
  local cand=""
  cand="$(ls -1 "$bam_dir/${sample}."*disco*.bam 2>/dev/null | head -n 1 || true)"
  [[ -n "$cand" ]] && { echo "$cand"; return 0; }
  cand="$(ls -1 "$bam_dir/${sample}."*discord*.bam 2>/dev/null | head -n 1 || true)"
  [[ -n "$cand" ]] && { echo "$cand"; return 0; }
  cand="$(ls -1 "$bam_dir"/*disco*.bam 2>/dev/null | head -n 1 || true)"
  [[ -n "$cand" ]] && { echo "$cand"; return 0; }
  cand="$(ls -1 "$bam_dir"/*discord*.bam 2>/dev/null | head -n 1 || true)"
  echo "$cand"
}

# -----------------------------
# Parse subcommand
# -----------------------------
if [[ $# -eq 0 || "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  usage
  exit 0
fi

subcmd="$1"
shift

# -----------------------------
# Parse args (shared)
# -----------------------------
parse_common_args() {
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --sample) SAMPLE="${2:-}"; shift 2 ;;
      --genome) GENOME="${2:-}"; shift 2 ;;
      --gff) GFF="${2:-}"; shift 2 ;;
      --out-root) OUT_ROOT="${2:-}"; shift 2 ;;
      --bam) BAM="${2:-}"; shift 2 ;;
      --bam-dir) BAM_DIR="${2:-}"; shift 2 ;;
      --model) MODEL_IN="${2:-}"; shift 2 ;;

      --start-at) START_AT="${2:-}"; shift 2 ;;
      --stop-after) STOP_AFTER="${2:-}"; shift 2 ;;
      --force) FORCE=1; shift ;;
      --keep-tmp) KEEP_TMP=1; shift ;;

      --train-window) TRAIN_WINDOW="${2:-}"; shift 2 ;;
      --train-mapq) TRAIN_MAPQ="${2:-}"; shift 2 ;;
      --train-threads) TRAIN_THREADS="${2:-}"; shift 2 ;;
      --train-drop-dup) TRAIN_DROP_DUP=1; shift ;;
      --train-no-drop-dup) TRAIN_DROP_DUP=0; shift ;;
      --train-exclude-contigs) TRAIN_EXCLUDE_CONTIGS+=("${2:-}"); shift 2 ;;

      --apply-binsize) APPLY_BINSIZE="${2:-}"; shift 2 ;;
      --apply-threads) APPLY_THREADS="${2:-}"; shift 2 ;;
      --apply-emit-raw) APPLY_EMIT_RAW=1; shift ;;
      --apply-primary-prefer) APPLY_PRIMARY_PREFER="${2:-}"; shift 2 ;;

      --r-modules) R_MODULES="${2:-}"; shift 2 ;;

      --tile-bp) TILE_BP="${2:-}"; shift 2 ;;
      --min-tile-bp) MIN_TILE_BP="${2:-}"; shift 2 ;;
      --exclude-bed) EXCLUDE_BED="${2:-}"; shift 2 ;;
      --lambda) LAMBDA="${2:-}"; shift 2 ;;
      --min-probes-per-seg) MIN_PROBES_PER_SEG="${2:-}"; shift 2 ;;
      --flank) FLANK="${2:-}"; shift 2 ;;
      --agg) AGG="${2:-}"; shift 2 ;;
      --confidence-method) CONF_METHOD="${2:-}"; shift 2 ;;

      --final-weak-ratio) FINAL_WEAK_RATIO="${2:-}"; shift 2 ;;
      --final-strong-ratio) FINAL_STRONG_RATIO="${2:-}"; shift 2 ;;
      --final-weak-z) FINAL_WEAK_Z="${2:-}"; shift 2 ;;
      --final-strong-z) FINAL_STRONG_Z="${2:-}"; shift 2 ;;
      --final-keep-mode) FINAL_KEEP_MODE="${2:-}"; shift 2 ;;
      --cn-base) CN_BASE="${2:-}"; shift 2 ;;

      --final-fuse) FINAL_FUSE=1; shift ;;
      --final-fuse-max-gap) FINAL_FUSE_MAX_GAP="${2:-}"; shift 2 ;;

      --final-count-flank) FINAL_COUNT_FLANK="${2:-}"; shift 2 ;;
      --final-count-all-alignments) FINAL_COUNT_ALL_ALIGN=1; shift ;;
      --final-count-drop-dup) FINAL_COUNT_DROP_DUP=1; shift ;;
      --final-count-include-supp) FINAL_COUNT_INCLUDE_SUPP=1; shift ;;
      --final-count-include-secondary) FINAL_COUNT_INCLUDE_SECONDARY=1; shift ;;

      --final2-balance) FINAL2_BALANCE="${2:-}"; shift 2 ;;
      --final2-rescue-z) FINAL2_RESCUE_Z="${2:-}"; shift 2 ;;
      --keep-suggested-only) FINAL2_KEEP_SUGGESTED_ONLY=1; shift ;;
      --keep-working-dir) FINAL2_KEEP_WORKING_DIR=1; shift ;;
      --verbose) FINAL2_VERBOSE=1; shift ;;

      --circular-contigs) CIRCULAR_CONTIGS="${2:-}"; shift 2 ;;
      --circular-contigs-file) CIRCULAR_CONTIGS_FILE="${2:-}"; shift 2 ;;

      --debug) FINAL_DEBUG=1; shift ;;

      -h|--help) usage; exit 0 ;;
      *) die "Unknown argument: $1" ;;
    esac
  done
}

# -----------------------------
# Step implementations
# -----------------------------
step_cov_model_train() {
  [[ -n "$BAM" ]] || die "--bam is required for Model_Train"
  [[ -f "$BAM" ]] || die "BAM not found: $BAM"
  [[ -n "$GENOME" ]] || die "--genome is required"

  local out="$WORK/train"
  mkdir -p "$out"

  local cmd=( bash "$HIDDEN_DIR/cov_model_train.sh"
    --bam "$BAM"
    --genome "$GENOME"
    --out-dir "$out"
    --window "$TRAIN_WINDOW"
    --mapq "$TRAIN_MAPQ"
    --threads "$TRAIN_THREADS"
  )

  if [[ "$TRAIN_DROP_DUP" -eq 1 ]]; then
    cmd+=( --drop-dup-for-train )
  fi

  local x
  for x in "${TRAIN_EXCLUDE_CONTIGS[@]}"; do
    cmd+=( --exclude-contigs "$x" )
  done

  cmd+=( --r-modules "$R_MODULES" )

  run_step "Model_Train" "$LOG/cov_model_train.log" "${cmd[@]}"
}

step_cov_model_apply() {
  [[ -n "$BAM_DIR" ]] || die "--bam-dir is required for Model_Apply"
  [[ -d "$BAM_DIR" ]] || die "BAM dir not found: $BAM_DIR"
  [[ -n "$GENOME" ]] || die "--genome is required"

  local model="$MODEL_IN"
  if [[ -z "$model" ]]; then
    model="$WORK/train/model.rds"
  fi
  [[ -f "$model" ]] || die "Model .rds not found: $model (provide --model or run Model_Train)"

  local out="$WORK/apply"
  mkdir -p "$out"

  local cmd=( bash "$HIDDEN_DIR/cov_model_apply.sh"
    --dir "$BAM_DIR"
    --sample-key "$SAMPLE"
    --model "$model"
    --genome "$GENOME"
    --out-dir "$out"
    --binsize "$APPLY_BINSIZE"
    --threads "$APPLY_THREADS"
    --primary-prefer "$APPLY_PRIMARY_PREFER"
    --r-modules "$R_MODULES"
  )

  if [[ "$APPLY_EMIT_RAW" -eq 1 ]]; then
    cmd+=( --emit-raw )
  fi
  if [[ "$KEEP_TMP" -eq 1 ]]; then
    cmd+=( --keep-tmp )
  fi

  run_step "Model_Apply" "$LOG/cov_model_apply.log" "${cmd[@]}"

  mkdir -p "$SAMPLE_DIR/corr_bw" "$SAMPLE_DIR/factors"
  rsync -a "$out/$SAMPLE/corr_bw/" "$SAMPLE_DIR/corr_bw/" || true
  rsync -a "$out/$SAMPLE/factors/" "$SAMPLE_DIR/factors/" || true

  if [[ -d "$out/$SAMPLE/raw_bw" ]]; then
    mkdir -p "$SAMPLE_DIR/raw_bw"
    rsync -a "$out/$SAMPLE/raw_bw/" "$SAMPLE_DIR/raw_bw/" || true
  fi
}

step_cov_bw_qc() {
  [[ -n "$BAM_DIR" ]] || die "--bam-dir is required for Model_Finalize"

  local out_tsv="$WORK/cov_bw_qc.summary.tsv"
  run_step "Model_Finalize" "$LOG/cov_bw_qc.log" bash "$HIDDEN_DIR/cov_bw_qc.sh" \
    --outdir "$WORK/apply" \
    --out "$out_tsv" \
    --bam-root "$BAM_DIR"

  cp "$out_tsv" "$SAMPLE_DIR/$SAMPLE.cov_bw_qc.summary.tsv"
}

step_make_probes() {
  [[ -n "$GENOME" ]] || die "--genome required"
  [[ -n "$GFF" ]] || die "--gff required"

  local corr_dir="$SAMPLE_DIR/corr_bw"
  [[ -d "$corr_dir" ]] || die "Missing corr_bw dir: $corr_dir (run Model_Apply first)"

  local primary_bw
  primary_bw="$(pick_primary_corr_bw "$corr_dir" "$SAMPLE" "$APPLY_BINSIZE")"
  [[ -n "$primary_bw" && -f "$primary_bw" ]] || die "Could not select primary corrected BW under: $corr_dir"
  msg "Primary BW selected: $primary_bw"

  local out="$WORK/probes"
  mkdir -p "$out"

  local py="$HIDDEN_DIR/make_probes.py"
  [[ -f "$py" ]] || die "Missing make_probes.py: $py"

  local runner
  runner="$(ensure_env_python)"

  local cmd=( python "$py"
    --genome "$GENOME"
    --gff "$GFF"
    --primary-bw "$primary_bw"
    --out-prefix "$out"
    --tile-bp "$TILE_BP"
    --min-tile-bp "$MIN_TILE_BP"
  )
  if [[ -n "$EXCLUDE_BED" ]]; then
    cmd+=( --exclude-bed "$EXCLUDE_BED" )
  fi

  if [[ -n "$runner" ]]; then
    run_step "CNV_Probe" "$LOG/make_probes.log" $runner "${cmd[@]}"
  else
    run_step "CNV_Probe" "$LOG/make_probes.log" "${cmd[@]}"
  fi
}

step_fused_lasso_segment() {
  local probe_tsv="$WORK/probes/probe_table.tsv.gz"
  [[ -f "$probe_tsv" ]] || die "Missing probe table (run CNV_Probe first): $probe_tsv"

  local out="$WORK/segment"
  mkdir -p "$out"

  local r_script="$HIDDEN_DIR/fused_lasso_segment.R"
  [[ -f "$r_script" ]] || die "Missing R script on this node: $r_script"

  msg "Running fused_lasso_segment (PELT) penalty(lambda)=$LAMBDA min_probes=$MIN_PROBES_PER_SEG"
  msg "R modules mode: ${R_MODULES:-off}"
  msg "R TMPDIR: $WORK/tmp_R"

  run_step "CNV_Segment" "$LOG/fused_lasso_segment.log" \
    run_r_file_hpc "$r_script" \
      --probe-tsv "$probe_tsv" \
      --out-prefix "$out" \
      --lambda "$LAMBDA" \
      --min-probes-per-seg "$MIN_PROBES_PER_SEG" \
      --tmpdir "$WORK/tmp_R"

  cp "$out/segments.tsv.gz" "$SAMPLE_DIR/$SAMPLE.segments.tsv.gz"
  cp "$out/segments.bed.gz" "$SAMPLE_DIR/$SAMPLE.segments.bed.gz"
}

step_boundary_support() {
  local boundaries="$WORK/segment/boundaries.bed.gz"
  [[ -f "$boundaries" ]] || die "Missing boundaries (run CNV_Segment first): $boundaries"

  local corr_dir="$SAMPLE_DIR/corr_bw"
  [[ -d "$corr_dir" ]] || die "Missing corr_bw dir: $corr_dir"

  local split_bw disco_bw
  split_bw="$(ls -1 "$corr_dir/${SAMPLE}."*split*"corr.bins${APPLY_BINSIZE}.bw" 2>/dev/null | head -n 1 || true)"
  disco_bw="$(ls -1 "$corr_dir/${SAMPLE}."*disco*"corr.bins${APPLY_BINSIZE}.bw" 2>/dev/null | head -n 1 || true)"

  [[ -n "$split_bw" && -f "$split_bw" ]] || die "Could not find splitters corr BW in: $corr_dir"
  [[ -n "$disco_bw" && -f "$disco_bw" ]] || die "Could not find discordant corr BW in: $corr_dir"

  msg "Split BW selected: $split_bw"
  msg "Disco BW selected: $disco_bw"

  local out="$WORK/boundary_support"
  mkdir -p "$out"

  local py="$HIDDEN_DIR/boundary_support_from_bigwigs.py"
  [[ -f "$py" ]] || die "Missing boundary_support_from_bigwigs.py: $py"

  local runner
  runner="$(ensure_env_python)"

  local out_tsv="$out/boundary_support.tsv.gz"
  local out_bg="$out/boundary_support.z.bedGraph.gz"

  local cmd=( python "$py"
    --boundaries-bed "$boundaries"
    --genome "$GENOME"
    --split-bw "$split_bw"
    --disco-bw "$disco_bw"
    --flank "$FLANK"
    --agg "$AGG"
    --out-tsv "$out_tsv"
    --out-bedgraph "$out_bg"
    --bedgraph-field combined_z_chr
  )

  if [[ -n "$runner" ]]; then
    run_step "CNV_Boundary" "$LOG/boundary_support_from_bigwigs.log" $runner "${cmd[@]}"
  else
    run_step "CNV_Boundary" "$LOG/boundary_support_from_bigwigs.log" "${cmd[@]}"
  fi
}

step_segment_confidence() {
  local seg_bed="$WORK/segment/segments.bed.gz"
  local bsup="$WORK/boundary_support/boundary_support.tsv.gz"
  [[ -f "$seg_bed" ]] || die "Missing segments.bed.gz: $seg_bed"
  [[ -f "$bsup" ]] || die "Missing boundary_support.tsv.gz: $bsup"

  local out="$WORK/confidence"
  mkdir -p "$out"

  local py="$HIDDEN_DIR/segment_confidence_from_boundaries.py"
  [[ -f "$py" ]] || die "Missing segment_confidence_from_boundaries.py: $py"

  local runner
  runner="$(ensure_env_python)"

  local out_tsv="$out/${SAMPLE}_segments.confidence.tsv.gz"
  local out_bed="$out/${SAMPLE}_segments.confidence.bed.gz"

  local cmd=( python "$py"
    --segments-bed "$seg_bed"
    --boundary-support-tsv "$bsup"
    --out-tsv "$out_tsv"
    --out-bed "$out_bed"
    --method "$CONF_METHOD"
    --score-from z
  )

  if [[ -n "$runner" ]]; then
    run_step "CNV_Confidence" "$LOG/segment_confidence_from_boundaries.log" $runner "${cmd[@]}"
  else
    run_step "CNV_Confidence" "$LOG/segment_confidence_from_boundaries.log" "${cmd[@]}"
  fi
}

step_cnv_finalize() {
  local seg_bed="$WORK/segment/segments.bed.gz"
  local bsup="$WORK/boundary_support/boundary_support.tsv.gz"
  [[ -f "$seg_bed" ]] || die "Missing segments.bed.gz: $seg_bed"
  [[ -f "$bsup" ]] || die "Missing boundary_support.tsv.gz: $bsup"

  local out="$WORK/final"
  mkdir -p "$out"

  local py="$HIDDEN_DIR/cnv_finalize_segments.py"
  [[ -f "$py" ]] || die "Missing cnv_finalize_segments.py: $py"

  local runner
  runner="$(ensure_env_python)"

  local out_calls="$out/segments.calls.tsv.gz"
  local out_y="$out/${SAMPLE}_cn.y.bedGraph.gz"
  local out_ratio="$out/${SAMPLE}_cn.ratio.bedGraph.gz"

  local primary_bam split_bam disco_bam
  primary_bam="$BAM"
  if [[ -z "$primary_bam" && -n "$BAM_DIR" && -d "$BAM_DIR" ]]; then
    primary_bam="$(pick_primary_bam "$BAM_DIR" "$SAMPLE")"
  fi

  split_bam=""
  disco_bam=""
  if [[ -n "$BAM_DIR" && -d "$BAM_DIR" ]]; then
    split_bam="$(pick_split_bam "$BAM_DIR" "$SAMPLE")"
    disco_bam="$(pick_disco_bam "$BAM_DIR" "$SAMPLE")"
  fi

  local cmd=( python "$py"
    --segments-bed "$seg_bed"
    --boundary-support-tsv "$bsup"
    --out-tsv "$out_calls"
    --weak-ratio "$FINAL_WEAK_RATIO"
    --strong-ratio "$FINAL_STRONG_RATIO"
    --confidence-method "$CONF_METHOD"
    --weak-z "$FINAL_WEAK_Z"
    --strong-z "$FINAL_STRONG_Z"
    --keep-mode "$FINAL_KEEP_MODE"
    --cn-base "$CN_BASE"
    --out-bedgraph-y "$out_y"
    --out-bedgraph-ratio "$out_ratio"
  )

  msg "[cnv_finalize_segments wrapper] circular args: csv='${CIRCULAR_CONTIGS:-}' file='${CIRCULAR_CONTIGS_FILE:-}'"

  if [[ -n "${CIRCULAR_CONTIGS:-}" ]]; then
    cmd+=( --circular-contigs "${CIRCULAR_CONTIGS}" )
  fi
  if [[ -n "${CIRCULAR_CONTIGS_FILE:-}" ]]; then
    cmd+=( --circular-contigs-file "${CIRCULAR_CONTIGS_FILE}" )
  fi

  if [[ "$FINAL_FUSE" -eq 1 ]]; then
    cmd+=( --fuse --fuse-max-gap "$FINAL_FUSE_MAX_GAP" )
  fi

  if [[ -n "$primary_bam" && -f "$primary_bam" ]]; then
    msg "Finalize counts: primary_bam=$primary_bam"
    cmd+=( --primary-bam "$primary_bam" )
  else
    msg "WARN: finalize counts: primary_bam not found (skipping seg_primary_reads)"
  fi

  if [[ -n "${split_bam:-}" && -f "${split_bam:-}" ]]; then
    msg "Finalize counts: split_bam=$split_bam"
    cmd+=( --split-bam "$split_bam" )
  else
    msg "WARN: finalize counts: split_bam not found (skipping split counts)"
  fi

  if [[ -n "${disco_bam:-}" && -f "${disco_bam:-}" ]]; then
    msg "Finalize counts: disco_bam=$disco_bam"
    cmd+=( --disco-bam "$disco_bam" )
  else
    msg "WARN: finalize counts: disco_bam not found (skipping disco counts)"
  fi

  cmd+=( --count-flank "$FINAL_COUNT_FLANK" )

  if [[ "$FINAL_COUNT_ALL_ALIGN" -eq 1 ]]; then
    cmd+=( --count-all-alignments )
  fi
  if [[ "$FINAL_COUNT_DROP_DUP" -eq 1 ]]; then
    cmd+=( --count-drop-dup )
  fi
  if [[ "$FINAL_COUNT_INCLUDE_SUPP" -eq 1 ]]; then
    cmd+=( --count-include-supp )
  fi
  if [[ "$FINAL_COUNT_INCLUDE_SECONDARY" -eq 1 ]]; then
    cmd+=( --count-include-secondary )
  fi
  if [[ "$FINAL_DEBUG" -eq 1 ]]; then
    cmd+=( --debug )
  fi

  if [[ -n "$runner" ]]; then
    run_step "CNV_Finalize" "$LOG/cnv_finalize_segments.log" $runner "${cmd[@]}"
  else
    run_step "CNV_Finalize" "$LOG/cnv_finalize_segments.log" "${cmd[@]}"
  fi

  local f bn
  for f in "$WORK/final/"*; do
    [[ -f "$f" ]] || continue
    bn="$(basename "$f")"
    cp "$f" "$SAMPLE_DIR/${SAMPLE}.${bn}"
  done
}

step_cnv_finalize_suggested() {
  local in_tsv="$WORK/final/segments.calls.tsv.gz"
  [[ -f "$in_tsv" ]] || die "Missing finalize calls TSV (run CNV_Finalize first): $in_tsv"

  local out="$WORK/final2"
  mkdir -p "$out"

  local py="$HIDDEN_DIR/cnv_finalize_suggested.py"
  [[ -f "$py" ]] || die "Missing cnv_finalize_suggested.py: $py"

  local runner
  runner="$(ensure_env_python)"

  local out_calls="$out/segments.calls.final2.tsv.gz"

  local cmd=( python "$py"
    --in-tsv "$in_tsv"
    --out-tsv "$out_calls"
    --sample "$SAMPLE"
    --balance "$FINAL2_BALANCE"
    --rescue-z "$FINAL2_RESCUE_Z"
    --weak-z "$FINAL_WEAK_Z"
    --strong-z "$FINAL_STRONG_Z"
  )

  if [[ "$FINAL2_KEEP_SUGGESTED_ONLY" -eq 1 ]]; then
    cmd+=( --keep-suggested-only )
  fi
  if [[ "$FINAL2_VERBOSE" -eq 1 ]]; then
    cmd+=( --verbose )
  fi
  if [[ "$FINAL_DEBUG" -eq 1 ]]; then
    cmd+=( --debug )
  fi

  if [[ -n "$runner" ]]; then
    run_step "CNV_Finalize2" "$LOG/cnv_finalize_suggested.log" $runner "${cmd[@]}"
  else
    run_step "CNV_Finalize2" "$LOG/cnv_finalize_suggested.log" "${cmd[@]}"
  fi

  local f bn
  for f in "$WORK/final2/"*; do
    [[ -f "$f" ]] || continue
    bn="$(basename "$f")"
    cp "$f" "$SAMPLE_DIR/${SAMPLE}.${bn}"
  done

  if [[ "$FINAL2_KEEP_WORKING_DIR" -ne 1 ]]; then
    rm -rf "$WORK/final2"
  fi
}

run_internal_step() {
  local internal="$1"
  case "$internal" in
    cov_model_train) step_cov_model_train ;;
    cov_model_apply) step_cov_model_apply ;;
    cov_bw_qc) step_cov_bw_qc ;;
    make_probes) step_make_probes ;;
    fused_lasso_segment) step_fused_lasso_segment ;;
    boundary_support_from_bigwigs) step_boundary_support ;;
    segment_confidence_from_boundaries) step_segment_confidence ;;
    cnv_finalize_segments) step_cnv_finalize ;;
    cnv_finalize_suggested) step_cnv_finalize_suggested ;;
    *) die "Unknown internal step: $internal" ;;
  esac
}

# -----------------------------
# Subcommand: step
# -----------------------------
if [[ "$subcmd" == "step" ]]; then
  [[ $# -ge 1 ]] || die "step requires a step name (alias, short alias, or internal name)"
  step_name="$1"
  shift

  parse_common_args "$@"
  validate_args

  [[ -n "$SAMPLE" ]] || die "--sample required"
  [[ -n "$OUT_ROOT" ]] || die "--out-root required"

  internal="$(resolve_step "$step_name")"
  [[ -n "$internal" ]] || die "Unknown step: $step_name"

  SAMPLE_DIR="$OUT_ROOT/$SAMPLE/$SAMPLE"
  WORK="$SAMPLE_DIR/_work"
  LOG="$SAMPLE_DIR/_log"
  mkdir -p "$SAMPLE_DIR" "$WORK" "$LOG"

  msg "Sample        : $SAMPLE"
  msg "Output root   : $OUT_ROOT"
  msg "Sample dir    : $SAMPLE_DIR"
  msg "Step          : $step_name -> $internal"

  run_internal_step "$internal"
  exit 0
fi

# -----------------------------
# Subcommand: run
# -----------------------------
if [[ "$subcmd" != "run" ]]; then
  die "Unknown subcommand: $subcmd (expected run|step)"
fi

parse_common_args "$@"
validate_args

[[ -n "$SAMPLE" ]] || die "--sample required"
[[ -n "$GENOME" ]] || die "--genome required"
[[ -n "$OUT_ROOT" ]] || die "--out-root required"

START_AT="$(resolve_step "$START_AT")"
STOP_AFTER="$(resolve_step "$STOP_AFTER")"
[[ -n "$START_AT" ]] || die "Unknown --start-at step"
[[ -n "$STOP_AFTER" ]] || die "Unknown --stop-after step"

SAMPLE_DIR="$OUT_ROOT/$SAMPLE/$SAMPLE"
WORK="$SAMPLE_DIR/_work"
LOG="$SAMPLE_DIR/_log"
mkdir -p "$SAMPLE_DIR" "$WORK" "$LOG"

msg "Sample        : $SAMPLE"
msg "Output root   : $OUT_ROOT"
msg "Sample dir    : $SAMPLE_DIR"
msg "Start at      : $START_AT"
msg "Stop after    : $STOP_AFTER"
msg "Circular args : csv='${CIRCULAR_CONTIGS:-}' file='${CIRCULAR_CONTIGS_FILE:-}'"

sidx="$(idx "$START_AT")"
eidx="$(idx "$STOP_AFTER")"
[[ "$sidx" -ge 0 ]] || die "Bad start step: $START_AT"
[[ "$eidx" -ge 0 ]] || die "Bad stop step: $STOP_AFTER"
[[ "$eidx" -ge "$sidx" ]] || die "--stop-after must be >= --start-at in pipeline order"

for i in $(seq "$sidx" "$eidx"); do
  internal="${steps[$i]}"

  if [[ "$FORCE" -ne 1 ]]; then
    case "$internal" in
      cov_model_train)
        [[ -f "$WORK/train/model.rds" ]] && { msg "Skip cov_model_train (exists)"; continue; }
        ;;
      cov_model_apply)
        [[ -d "$SAMPLE_DIR/corr_bw" && -n "$(ls -1 "$SAMPLE_DIR/corr_bw"/*.bw 2>/dev/null | head -n 1 || true)" ]] && { msg "Skip cov_model_apply (corr_bw exists)"; continue; }
        ;;
      cov_bw_qc)
        [[ -f "$SAMPLE_DIR/$SAMPLE.cov_bw_qc.summary.tsv" ]] && { msg "Skip cov_bw_qc (summary exists)"; continue; }
        ;;
      make_probes)
        [[ -f "$WORK/probes/probe_table.tsv.gz" ]] && { msg "Skip make_probes (probe table exists)"; continue; }
        ;;
      fused_lasso_segment)
        [[ -f "$WORK/segment/segments.bed.gz" ]] && { msg "Skip fused_lasso_segment (segments exist)"; continue; }
        ;;
      boundary_support_from_bigwigs)
        [[ -f "$WORK/boundary_support/boundary_support.tsv.gz" ]] && { msg "Skip boundary_support (exists)"; continue; }
        ;;
      segment_confidence_from_boundaries)
        [[ -f "$WORK/confidence/${SAMPLE}_segments.confidence.tsv.gz" ]] && { msg "Skip segment_confidence (exists)"; continue; }
        ;;
      cnv_finalize_segments)
        [[ -f "$WORK/final/segments.calls.tsv.gz" ]] && { msg "Skip cnv_finalize (calls exist)"; continue; }
        ;;
      cnv_finalize_suggested)
        [[ -f "$SAMPLE_DIR/${SAMPLE}.segments.calls.final2.tsv.gz" ]] && { msg "Skip cnv_finalize_suggested (final2 exists)"; continue; }
        ;;
    esac
  fi

  run_internal_step "$internal"
done

msg "CLOverCNV run complete: $SAMPLE_DIR"
