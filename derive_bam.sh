#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# bam_split_types.sh
# Split a coordinate-sorted BAM into derived BAMs:
#   - unmapped (0x4)
#   - supplementary (0x800)   [EXCLUSIVE: excludes splitters]
#   - primary-only            [exclude 0x100+0x800]
#   - secondary-only (0x100)
#   - splitters               [heuristic: SA tag]
#   - discordant              [heuristic: paired, mapped, primary, not proper pair, excludes splitters]
#
# Batch mode:
#   --bam-glob "glob/pattern/*.bam"
#
# Notes:
# - "splitters" and "discordant" here are best-effort proxies from a final BAM.
#   For exact samblaster outputs, generate them during alignment (samblaster step).
# ============================================================

say()  { echo "[bam_split_types] $*"; }
warn() { echo "[bam_split_types] WARN: $*" >&2; }
die()  { echo "[bam_split_types] ERROR: $*" >&2; exit 1; }
have() { command -v "$1" >/dev/null 2>&1; }

usage() {
  cat <<'EOF'
Usage:
  bam_split_types.sh --bam aligned.sorted.bam --outdir DIR [options]
  bam_split_types.sh --bam-glob "*/path/*.bam" --outdir DIR [options]

Required (choose one):
  --bam FILE            Coordinate-sorted BAM
  --bam-glob STR        Glob that expands to BAM(s) (quote it!)

Required:
  --outdir DIR          Output directory

Options:
  --prefix STR          Output prefix (single-BAM mode only; default: basename of --bam minus .bam)
  --threads INT         Threads for samtools sort/view (default: 4)
  --tmpdir DIR          Temp dir root for samtools sort (default: OUTDIR/tmp)
  --resume              Skip outputs that already have .bam + .bai
  --keep-unsorted       Keep intermediate unsorted BAMs (default: delete)
  --gzip-other          Gzip any non-.gz files in OUTDIR at end (off by default)

Emit toggles (defaults all ON):
  --no-unmapped
  --no-supplementary
  --no-primary
  --no-secondary
  --no-splitters
  --no-discordant

MAPQ filters (defaults: none):
  --min-mapq-primary INT
  --min-mapq-supp INT
  --min-mapq-secondary INT
  --min-mapq-splitters INT
  --min-mapq-discordant INT

Splitters heuristic:
  A record is considered a "splitter" if it has SA:Z: tag.
  (This typically captures split-read evidence; includes many supplementary records.)

Supplementary exclusivity:
  supplementary output = flag 0x800 AND NOT splitters (no SA tag)

Discordant heuristic (paired-end):
  paired (0x1), mapped read (not 0x4), mapped mate (not 0x8),
  primary only (exclude 0x100+0x800), NOT properly paired (not 0x2),
  NOT splitters (no SA tag)

Batch output layout:
  In --bam-glob mode, each input BAM writes to its own subdir:
    OUTDIR/<prefix>/
  In --bam mode, writes directly to OUTDIR/ unless you set --prefix.

EOF
}

# -----------------------------
# Defaults
# -----------------------------
THREADS=4
RESUME=0
KEEP_UNSORTED=0
GZIP_OTHER=0

EMIT_UNMAPPED=1
EMIT_SUPP=1
EMIT_PRIMARY=1
EMIT_SECONDARY=1
EMIT_SPLITTERS=1
EMIT_DISCORDANT=1

MIN_MAPQ_PRIMARY=""
MIN_MAPQ_SUPP=""
MIN_MAPQ_SECONDARY=""
MIN_MAPQ_SPLITTERS=""
MIN_MAPQ_DISCORDANT=""

BAM=""
BAM_GLOB=""
OUTDIR=""
PREFIX=""
TMPDIR=""

# -----------------------------
# Arg parse
# -----------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --bam) BAM="$2"; shift 2;;
    --bam-glob) BAM_GLOB="$2"; shift 2;;
    --outdir) OUTDIR="$2"; shift 2;;
    --prefix) PREFIX="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --tmpdir) TMPDIR="$2"; shift 2;;
    --resume) RESUME=1; shift;;
    --keep-unsorted) KEEP_UNSORTED=1; shift;;
    --gzip-other) GZIP_OTHER=1; shift;;

    --no-unmapped) EMIT_UNMAPPED=0; shift;;
    --no-supplementary) EMIT_SUPP=0; shift;;
    --no-primary) EMIT_PRIMARY=0; shift;;
    --no-secondary) EMIT_SECONDARY=0; shift;;
    --no-splitters) EMIT_SPLITTERS=0; shift;;
    --no-discordant) EMIT_DISCORDANT=0; shift;;

    --min-mapq-primary) MIN_MAPQ_PRIMARY="$2"; shift 2;;
    --min-mapq-supp) MIN_MAPQ_SUPP="$2"; shift 2;;
    --min-mapq-secondary) MIN_MAPQ_SECONDARY="$2"; shift 2;;
    --min-mapq-splitters) MIN_MAPQ_SPLITTERS="$2"; shift 2;;
    --min-mapq-discordant) MIN_MAPQ_DISCORDANT="$2"; shift 2;;

    -h|--help) usage; exit 0;;
    *) die "Unknown argument: $1";;
  esac
done

# -----------------------------
# Validate mode
# -----------------------------
if [[ -n "$BAM" && -n "$BAM_GLOB" ]]; then
  die "Use only one of --bam or --bam-glob"
fi
if [[ -z "$BAM" && -z "$BAM_GLOB" ]]; then
  usage
  die "One of --bam or --bam-glob is required"
fi
[[ -n "$OUTDIR" ]] || { usage; die "--outdir is required"; }

have samtools || die "samtools not found in PATH"

mkdir -p "$OUTDIR"

# tmp root
if [[ -z "$TMPDIR" ]]; then
  TMPDIR="$OUTDIR/tmp"
fi
mkdir -p "$TMPDIR"

# -----------------------------
# Expand inputs
# -----------------------------
BAMS=()
if [[ -n "$BAM_GLOB" ]]; then
  shopt -s nullglob
  BAMS=( $BAM_GLOB )
  shopt -u nullglob
  [[ "${#BAMS[@]}" -gt 0 ]] || die "Glob matched zero BAMs: $BAM_GLOB"
else
  [[ -f "$BAM" ]] || die "BAM not found: $BAM"
  BAMS=( "$BAM" )
fi

# -----------------------------
# Helpers
# -----------------------------
bam_complete() {
  local b="$1"
  [[ -s "$b" && ( -s "$b.bai" || -s "${b%.bam}.bai" ) ]]
}

sort_and_index_bam() {
  local in_bam="$1"
  local out_bam="$2"
  local tmp_base="$3"
  say "sort+index: $(basename "$in_bam") -> $(basename "$out_bam")"
  samtools sort -@ "$THREADS" -T "$tmp_base" -o "$out_bam" "$in_bam"
  samtools index -@ "$THREADS" "$out_bam"
  if [[ "$KEEP_UNSORTED" -eq 0 ]]; then
    rm -f "$in_bam" || true
  fi
}

# build optional MAPQ args array (empty if unset)
view_mapq_args() {
  local mq="$1"
  VIEW_MAPQ_ARGS=()
  if [[ -n "${mq}" ]]; then
    VIEW_MAPQ_ARGS=(-q "$mq")
  fi
}

failures=0

# -----------------------------
# Main loop
# -----------------------------
for BAM in "${BAMS[@]}"; do
  say "Processing: $BAM"
  if [[ ! -f "$BAM" ]]; then
    warn "Missing BAM, skip: $BAM"
    failures=$((failures+1))
    continue
  fi

  if [[ ! -f "$BAM.bai" && ! -f "${BAM%.bam}.bai" ]]; then
    warn "No BAM index detected next to input: $BAM (recommended but not required for view)"
  fi

  bn="$(basename "$BAM")"
  bam_prefix="${bn%.bam}"
  if [[ -n "$BAM_GLOB" ]]; then
    PREFIX_THIS="$bam_prefix"
  else
    PREFIX_THIS="${PREFIX:-$bam_prefix}"
  fi

  if [[ -n "$BAM_GLOB" ]]; then
    BAM_OUTDIR="$OUTDIR/$PREFIX_THIS"
  else
    BAM_OUTDIR="$OUTDIR"
  fi
  mkdir -p "$BAM_OUTDIR"

  TMPDIR_THIS="$TMPDIR/$PREFIX_THIS"
  mkdir -p "$TMPDIR_THIS"

  # -----------------------------
  # Output paths
  # -----------------------------
  UNMAPPED_UNSORTED="$BAM_OUTDIR/${PREFIX_THIS}.unmapped.unsorted.bam"
  UNMAPPED_SORTED="$BAM_OUTDIR/${PREFIX_THIS}.unmapped.sorted.bam"

  SUPP_UNSORTED="$BAM_OUTDIR/${PREFIX_THIS}.supplementary.unsorted.bam"
  SUPP_SORTED="$BAM_OUTDIR/${PREFIX_THIS}.supplementary.sorted.bam"

  PRIMARY_UNSORTED="$BAM_OUTDIR/${PREFIX_THIS}.primary.unsorted.bam"
  PRIMARY_SORTED="$BAM_OUTDIR/${PREFIX_THIS}.primary.sorted.bam"

  SECONDARY_UNSORTED="$BAM_OUTDIR/${PREFIX_THIS}.secondary.unsorted.bam"
  SECONDARY_SORTED="$BAM_OUTDIR/${PREFIX_THIS}.secondary.sorted.bam"

  SPLITTERS_UNSORTED="$BAM_OUTDIR/${PREFIX_THIS}.splitters.unsorted.bam"
  SPLITTERS_SORTED="$BAM_OUTDIR/${PREFIX_THIS}.splitters.sorted.bam"

  DISCORDANT_UNSORTED="$BAM_OUTDIR/${PREFIX_THIS}.discordant.unsorted.bam"
  DISCORDANT_SORTED="$BAM_OUTDIR/${PREFIX_THIS}.discordant.sorted.bam"

  # ============================================================
  # Splitters (SA-tag heuristic)
  # ============================================================
  if [[ "$EMIT_SPLITTERS" -eq 1 ]]; then
    if [[ "$RESUME" -eq 1 ]] && bam_complete "$SPLITTERS_SORTED"; then
      say "Exists (skip): $SPLITTERS_SORTED"
    else
      say "Writing splitters (SA tag heuristic): $SPLITTERS_SORTED"
      view_mapq_args "$MIN_MAPQ_SPLITTERS"
      # Keep header (-h). Match SA:Z: in optional fields.
      samtools view -@ "$THREADS" -h "${VIEW_MAPQ_ARGS[@]}" "$BAM" \
        | awk 'BEGIN{OFS="\t"} /^@/ {print; next} $0 ~ /\tSA:Z:/ {print}' \
        | samtools view -@ "$THREADS" -b -o "$SPLITTERS_UNSORTED" -
      sort_and_index_bam "$SPLITTERS_UNSORTED" "$SPLITTERS_SORTED" "$TMPDIR_THIS/$(basename "$SPLITTERS_SORTED").tmp"
    fi
  fi

  # ============================================================
  # Discordant (paired, mapped, primary, NOT proper pair, NOT splitters)
  # ============================================================
  if [[ "$EMIT_DISCORDANT" -eq 1 ]]; then
    if [[ "$RESUME" -eq 1 ]] && bam_complete "$DISCORDANT_SORTED"; then
      say "Exists (skip): $DISCORDANT_SORTED"
    else
      say "Writing discordant (heuristic: paired+mapped+primary+!proper_pair, excludes SA): $DISCORDANT_SORTED"
      view_mapq_args "$MIN_MAPQ_DISCORDANT"
      # Start with primary only and mapped read+mate:
      #   -f 1  (paired)
      #   -F 4  (exclude unmapped read)
      #   -F 8  (exclude unmapped mate)
      #   -F 2304 (exclude secondary+supp)
      # Then exclude properly paired (0x2) by filtering flag bit with awk.
      # Also exclude splitters by removing SA-tagged records.
      samtools view -@ "$THREADS" -h -f 1 -F 2316 "${VIEW_MAPQ_ARGS[@]}" "$BAM" \
        | awk '
            BEGIN{OFS="\t"}
            /^@/ {print; next}
            # fields: QNAME FLAG RNAME POS MAPQ CIGAR RNEXT PNEXT TLEN SEQ QUAL [TAGS...]
            {
              flag=$2
              # exclude properly paired (0x2)
              if (and(flag,2)) next
              # exclude splitters (SA tag)
              if ($0 ~ /\tSA:Z:/) next
              print
            }' \
        | samtools view -@ "$THREADS" -b -o "$DISCORDANT_UNSORTED" -
      sort_and_index_bam "$DISCORDANT_UNSORTED" "$DISCORDANT_SORTED" "$TMPDIR_THIS/$(basename "$DISCORDANT_SORTED").tmp"
    fi
  fi

  # ============================================================
  # Unmapped
  # ============================================================
  if [[ "$EMIT_UNMAPPED" -eq 1 ]]; then
    if [[ "$RESUME" -eq 1 ]] && bam_complete "$UNMAPPED_SORTED"; then
      say "Exists (skip): $UNMAPPED_SORTED"
    else
      say "Writing unmapped (0x4): $UNMAPPED_SORTED"
      samtools view -@ "$THREADS" -b -f 4 "$BAM" > "$UNMAPPED_UNSORTED"
      sort_and_index_bam "$UNMAPPED_UNSORTED" "$UNMAPPED_SORTED" "$TMPDIR_THIS/$(basename "$UNMAPPED_SORTED").tmp"
    fi
  fi

  # ============================================================
  # Supplementary EXCLUSIVE (0x800 but NOT splitters)
  # ============================================================
  if [[ "$EMIT_SUPP" -eq 1 ]]; then
    if [[ "$RESUME" -eq 1 ]] && bam_complete "$SUPP_SORTED"; then
      say "Exists (skip): $SUPP_SORTED"
    else
      say "Writing supplementary-exclusive (0x800 excluding SA): $SUPP_SORTED"
      view_mapq_args "$MIN_MAPQ_SUPP"
      samtools view -@ "$THREADS" -h -f 2048 "${VIEW_MAPQ_ARGS[@]}" "$BAM" \
        | awk 'BEGIN{OFS="\t"} /^@/ {print; next} $0 !~ /\tSA:Z:/ {print}' \
        | samtools view -@ "$THREADS" -b -o "$SUPP_UNSORTED" -
      sort_and_index_bam "$SUPP_UNSORTED" "$SUPP_SORTED" "$TMPDIR_THIS/$(basename "$SUPP_SORTED").tmp"
    fi
  fi

  # ============================================================
  # Primary-only (exclude 0x100+0x800)
  # ============================================================
  if [[ "$EMIT_PRIMARY" -eq 1 ]]; then
    if [[ "$RESUME" -eq 1 ]] && bam_complete "$PRIMARY_SORTED"; then
      say "Exists (skip): $PRIMARY_SORTED"
    else
      say "Writing primary-only (exclude 0x100+0x800): $PRIMARY_SORTED"
      view_mapq_args "$MIN_MAPQ_PRIMARY"
      samtools view -@ "$THREADS" -b -F 2304 "${VIEW_MAPQ_ARGS[@]}" "$BAM" > "$PRIMARY_UNSORTED"
      sort_and_index_bam "$PRIMARY_UNSORTED" "$PRIMARY_SORTED" "$TMPDIR_THIS/$(basename "$PRIMARY_SORTED").tmp"
    fi
  fi

  # ============================================================
  # Secondary-only (0x100)
  # ============================================================
  if [[ "$EMIT_SECONDARY" -eq 1 ]]; then
    if [[ "$RESUME" -eq 1 ]] && bam_complete "$SECONDARY_SORTED"; then
      say "Exists (skip): $SECONDARY_SORTED"
    else
      say "Writing secondary-only (0x100): $SECONDARY_SORTED"
      view_mapq_args "$MIN_MAPQ_SECONDARY"
      samtools view -@ "$THREADS" -b -f 256 "${VIEW_MAPQ_ARGS[@]}" "$BAM" > "$SECONDARY_UNSORTED"
      sort_and_index_bam "$SECONDARY_UNSORTED" "$SECONDARY_SORTED" "$TMPDIR_THIS/$(basename "$SECONDARY_SORTED").tmp"
    fi
  fi

done

# Optional: gzip any non-.gz files in OUTDIR (excluding .bam/.bai)
if [[ "$GZIP_OTHER" -eq 1 ]]; then
  say "Gzipping non-.gz regular files in OUTDIR (excluding .bam/.bai)"
  shopt -s nullglob
  for f in "$OUTDIR"/*; do
    [[ -f "$f" ]] || continue
    [[ "$f" == *.gz ]] && continue
    [[ "$f" == *.bam || "$f" == *.bai ]] && continue
    gzip -f "$f"
  done
  shopt -u nullglob
fi

if [[ "$failures" -gt 0 ]]; then
  warn "Completed with $failures input BAM(s) skipped/missing."
fi

say "Done."
