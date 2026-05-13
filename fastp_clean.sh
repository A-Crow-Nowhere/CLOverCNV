#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# MAPI MODULE: fastp_clean
# Adapter trim + quality trim + length filter + pair integrity
# ============================================================

SELF="$(readlink -f "${BASH_SOURCE[0]}")"
BIN_DIR="$(cd "$(dirname "$SELF")" && pwd)"
MAPI_ROOT="$(cd "$BIN_DIR/../.." && pwd)"

usage() {
  cat <<'EOF'
mapi fastp_clean --r1 R1.fq.gz --r2 R2.fq.gz [options]

Required:
  --r1 FILE                 Input R1 fastq(.gz)
  --r2 FILE                 Input R2 fastq(.gz)

Output control:
  --outdir DIR              Output directory (default: ./fastp_clean_out)
  --prefix STR              Output prefix (default: derived from R1 basename)
  --overwrite               Overwrite existing outputs

Core trimming/filtering (maps to fastp):
  --detect-adapters         Let fastp auto-detect adapters for PE (adds --detect_adapter_for_pe)
  --adapter1 SEQ            Adapter sequence for R1
  --adapter2 SEQ            Adapter sequence for R2

  --cut-front               Enable 5' sliding-window trimming (adds --cut_front)
  --cut-tail                Enable 3' sliding-window trimming (adds --cut_tail)
  --cut-window-size INT     Sliding window size (default: 4)
  --cut-mean-quality INT    Mean quality threshold (default: 20)

  --qualified-phred INT     Phred for "qualified" bases (fastp --qualified_quality_phred, default: 15)
  --unqualified-limit INT   Max % unqualified bases (fastp --unqualified_percent_limit, default: 40)

  --length-required INT     Minimum read length after trimming (fastp --length_required, default: 50)
  --length-limit INT        Maximum read length (fastp --length_limit, default: 0 = no limit)
  --n-base-limit INT        Max number of Ns allowed per read (fastp --n_base_limit, default: 5)

Optional "obvious garbage" preset:
  --garbage-filter          Apply a conservative preset:
                              --n-base-limit 5
                              --length-required 50
                              --low-complexity-filter (with threshold 30)


Garbage filter behavior note (use it if you just dont know what you are doing):
  --garbage-filter applies a conservative preset for structural sanity
  (Ns, minimum length, low complexity). When enabled, these preset values
  override user-supplied values for the same options (e.g. --n-base-limit,
  --length-required) unless they are explicitly re-overridden via
  --fastp-extra (last-wins passthru). Quality-based filters
  (--qualified-phred, --unqualified-limit) are NOT modified by
  --garbage-filter.

Low complexity:
  --low-complexity-filter   Enable fastp low complexity filter
  --complexity-threshold N  fastp --complexity_threshold (default: 30)

Passthru:
  --fastp-extra "ARGS..."   Extra args appended verbatim to fastp (last-wins)

Examples:
  mapi fastp_clean --r1 sample_R1.fq.gz --r2 sample_R2.fq.gz --detect-adapters --cut-front --cut-tail
  mapi fastp_clean --r1 R1.fq.gz --r2 R2.fq.gz --garbage-filter --fastp-extra "--thread 16"

Notes:
- This runs fastp in PE mode, so pair synchronization is preserved automatically.
- Outputs:
    <outdir>/<prefix>.clean.R1.fq.gz
    <outdir>/<prefix>.clean.R2.fq.gz
    <outdir>/<prefix>.fastp.json
    <outdir>/<prefix>.fastp.html
EOF
}

die() { echo "ERROR: $*" >&2; exit 1; }

# --------------------------
# Defaults
# --------------------------
R1=""
R2=""
OUTDIR="./fastp_clean_out"
PREFIX=""
OVERWRITE=0

DETECT_ADAPTERS=0
ADAPTER1=""
ADAPTER2=""

CUT_FRONT=0
CUT_TAIL=0
CUT_WINDOW_SIZE=4
CUT_MEAN_QUAL=20

QUALIFIED_PHRED=15
UNQUAL_LIMIT=40

LENGTH_REQUIRED=50
LENGTH_LIMIT=0
N_BASE_LIMIT=5

GARBAGE_FILTER=0
LOW_COMPLEXITY=0
COMPLEXITY_THRESHOLD=30

FASTP_EXTRA=""

# --------------------------
# Arg parse
# --------------------------
if [[ $# -eq 0 ]]; then usage; exit 1; fi

while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help) usage; exit 0 ;;
    --r1) R1="${2:-}"; shift 2 ;;
    --r2) R2="${2:-}"; shift 2 ;;

    --outdir) OUTDIR="${2:-}"; shift 2 ;;
    --prefix) PREFIX="${2:-}"; shift 2 ;;
    --overwrite) OVERWRITE=1; shift ;;

    --detect-adapters) DETECT_ADAPTERS=1; shift ;;
    --adapter1) ADAPTER1="${2:-}"; shift 2 ;;
    --adapter2) ADAPTER2="${2:-}"; shift 2 ;;

    --cut-front) CUT_FRONT=1; shift ;;
    --cut-tail) CUT_TAIL=1; shift ;;
    --cut-window-size) CUT_WINDOW_SIZE="${2:-}"; shift 2 ;;
    --cut-mean-quality) CUT_MEAN_QUAL="${2:-}"; shift 2 ;;

    --qualified-phred) QUALIFIED_PHRED="${2:-}"; shift 2 ;;
    --unqualified-limit) UNQUAL_LIMIT="${2:-}"; shift 2 ;;

    --length-required) LENGTH_REQUIRED="${2:-}"; shift 2 ;;
    --length-limit) LENGTH_LIMIT="${2:-}"; shift 2 ;;
    --n-base-limit) N_BASE_LIMIT="${2:-}"; shift 2 ;;

    --garbage-filter) GARBAGE_FILTER=1; shift ;;
    --low-complexity-filter) LOW_COMPLEXITY=1; shift ;;
    --complexity-threshold) COMPLEXITY_THRESHOLD="${2:-}"; shift 2 ;;

    --fastp-extra) FASTP_EXTRA="${2:-}"; shift 2 ;;

    *) die "Unknown option: $1 (use --help)" ;;
  esac
done

# --------------------------
# Validate
# --------------------------
[[ -n "$R1" ]] || die "Missing --r1"
[[ -n "$R2" ]] || die "Missing --r2"
[[ -f "$R1" ]] || die "R1 not found: $R1"
[[ -f "$R2" ]] || die "R2 not found: $R2"

mkdir -p "$OUTDIR"

if [[ -z "$PREFIX" ]]; then
  b="$(basename "$R1")"
  # strip common suffixes
  b="${b%.fastq.gz}"; b="${b%.fq.gz}"; b="${b%.fastq}"; b="${b%.fq}"
  b="${b%_R1}"; b="${b%_1}"; b="${b%.R1}"; b="${b%.1}"
  PREFIX="$b"
fi

OUT1="$OUTDIR/${PREFIX}.clean.R1.fq.gz"
OUT2="$OUTDIR/${PREFIX}.clean.R2.fq.gz"
JSON="$OUTDIR/${PREFIX}.fastp.json"
HTML="$OUTDIR/${PREFIX}.fastp.html"

if [[ $OVERWRITE -eq 0 ]]; then
  [[ ! -e "$OUT1" ]] || die "Output exists: $OUT1 (use --overwrite)"
  [[ ! -e "$OUT2" ]] || die "Output exists: $OUT2 (use --overwrite)"
  [[ ! -e "$JSON" ]] || die "Output exists: $JSON (use --overwrite)"
  [[ ! -e "$HTML" ]] || die "Output exists: $HTML (use --overwrite)"
fi

# --------------------------
# Garbage filter preset
# --------------------------
if [[ $GARBAGE_FILTER -eq 1 ]]; then
  # Conservative defaults; user flags still override (last-wins via FASTP_EXTRA)
  N_BASE_LIMIT=5
  LENGTH_REQUIRED=50
  LOW_COMPLEXITY=1
  COMPLEXITY_THRESHOLD=30
fi

# --------------------------
# Build fastp command
# --------------------------
cmd=(fastp
  --in1 "$R1" --in2 "$R2"
  --out1 "$OUT1" --out2 "$OUT2"
  --json "$JSON" --html "$HTML"
  --qualified_quality_phred "$QUALIFIED_PHRED"
  --unqualified_percent_limit "$UNQUAL_LIMIT"
  --length_required "$LENGTH_REQUIRED"
  --n_base_limit "$N_BASE_LIMIT"
)

# optional max length
if [[ "$LENGTH_LIMIT" != "0" ]]; then
  cmd+=(--length_limit "$LENGTH_LIMIT")
fi

# adapter handling
if [[ $DETECT_ADAPTERS -eq 1 ]]; then
  cmd+=(--detect_adapter_for_pe)
fi
if [[ -n "$ADAPTER1" ]]; then
  cmd+=(--adapter_sequence "$ADAPTER1")
fi
if [[ -n "$ADAPTER2" ]]; then
  cmd+=(--adapter_sequence_r2 "$ADAPTER2")
fi

# quality trimming
if [[ $CUT_FRONT -eq 1 ]]; then cmd+=(--cut_front); fi
if [[ $CUT_TAIL -eq 1 ]]; then cmd+=(--cut_tail); fi
cmd+=(--cut_window_size "$CUT_WINDOW_SIZE" --cut_mean_quality "$CUT_MEAN_QUAL")

# low complexity
if [[ $LOW_COMPLEXITY -eq 1 ]]; then
  cmd+=(--low_complexity_filter --complexity_threshold "$COMPLEXITY_THRESHOLD")
fi

# passthru (last-wins)
if [[ -n "$FASTP_EXTRA" ]]; then
  # shellcheck disable=SC2206
  extra_arr=($FASTP_EXTRA)
  cmd+=("${extra_arr[@]}")
fi

echo "==> [fastp_clean] Running:"
printf '    %q ' "${cmd[@]}"; echo
"${cmd[@]}"

echo "==> [fastp_clean] Done."
echo "    R1: $OUT1"
echo "    R2: $OUT2"
echo "    JSON: $JSON"
echo "    HTML: $HTML"
