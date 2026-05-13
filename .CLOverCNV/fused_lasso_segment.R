#!/usr/bin/env Rscript

suppressWarnings(suppressMessages({
  library(data.table)
  library(changepoint)
}))

usage <- function() {
  cat("
Usage:
  fused_lasso_segment.R --probe-tsv <probe_table.tsv.gz> --out-prefix <OUTDIR>
    [--lambda <penalty>] [--min-probes-per-seg 5] [--tmpdir <DIR>]

Notes:
  This v1 uses changepoint::cpt.mean(method='PELT'). --lambda is treated as penalty.
  Outputs:
    <out-prefix>/segments.tsv.gz
    <out-prefix>/segments.bed.gz
    <out-prefix>/boundaries.bed.gz
    <out-prefix>/boundaries.tsv.gz
\n", file=stderr())
}

arg_val <- function(args, key, default=NA) {
  hit <- which(args == key)
  if (length(hit) == 0) return(default)
  if (hit == length(args)) stop(paste("Missing value after", key))
  args[hit + 1]
}
has_flag <- function(args, key) any(args == key)

args <- commandArgs(trailingOnly=TRUE)
if (length(args) == 0 || has_flag(args, "-h") || has_flag(args, "--help")) {
  usage(); quit(status=0)
}

probe_tsv <- arg_val(args, "--probe-tsv", NA)
out_prefix <- arg_val(args, "--out-prefix", NA)
lambda <- as.numeric(arg_val(args, "--lambda", 25))
min_probes <- as.integer(arg_val(args, "--min-probes-per-seg", 5))
tmpdir <- arg_val(args, "--tmpdir", NA)

if (is.na(probe_tsv) || is.na(out_prefix)) {
  usage(); stop("ERROR: --probe-tsv and --out-prefix are required")
}

if (!is.na(tmpdir)) {
  dir.create(tmpdir, showWarnings=FALSE, recursive=TRUE)
  Sys.setenv(TMPDIR=tmpdir)
}

dir.create(out_prefix, showWarnings=FALSE, recursive=TRUE)

message(sprintf("[segment] probe_tsv: %s", probe_tsv))
message(sprintf("[segment] out_prefix: %s", out_prefix))
message(sprintf("[segment] penalty(lambda): %s", lambda))
message(sprintf("[segment] min_probes_per_seg: %s", min_probes))

dt <- fread(probe_tsv)

needed <- c("chrom","start","end","y")
missing <- setdiff(needed, names(dt))
if (length(missing) > 0) stop(paste("ERROR: probe table missing cols:", paste(missing, collapse=",")))

dt[, start := as.integer(start)]
dt[, end := as.integer(end)]
dt[, y := as.numeric(y)]
setorder(dt, chrom, start, end)

merge_small_segments <- function(seg_dt, min_probes=5) {
  if (nrow(seg_dt) == 0) return(seg_dt)
  changed <- TRUE
  while (changed) {
    changed <- FALSE
    small_idx <- which(seg_dt$n_probes < min_probes)
    if (length(small_idx) == 0) break
    i <- small_idx[1]
    if (nrow(seg_dt) == 1) break

    if (i == 1) {
      j <- 2
    } else if (i == nrow(seg_dt)) {
      j <- nrow(seg_dt) - 1
    } else {
      left <- i - 1
      right <- i + 1
      dl <- abs(seg_dt$y[i] - seg_dt$y[left])
      dr <- abs(seg_dt$y[i] - seg_dt$y[right])
      j <- ifelse(dl <= dr, left, right)
    }

    s <- min(seg_dt$start[i], seg_dt$start[j])
    e <- max(seg_dt$end[i], seg_dt$end[j])
    y_new <- (seg_dt$y[i]*seg_dt$n_probes[i] + seg_dt$y[j]*seg_dt$n_probes[j]) /
             (seg_dt$n_probes[i] + seg_dt$n_probes[j])
    n_new <- seg_dt$n_probes[i] + seg_dt$n_probes[j]

    seg_dt$start[j] <- s
    seg_dt$end[j] <- e
    seg_dt$y[j] <- y_new
    seg_dt$n_probes[j] <- n_new

    seg_dt <- seg_dt[-i,]
    seg_dt <- seg_dt[order(seg_dt$start, seg_dt$end),]
    rownames(seg_dt) <- NULL
    changed <- TRUE
  }
  seg_dt
}

segments_all <- list()
boundaries_all <- list()

chrs <- unique(dt$chrom)

for (chr in chrs) {
  sub <- dt[chrom == chr]
  if (nrow(sub) == 0) next

  yvec <- sub$y
  n <- length(yvec)

  if (n < 2) {
    seg <- data.table(
      chrom=chr,
      start=min(sub$start),
      end=max(sub$end),
      y=mean(sub$y),
      n_probes=n
    )
    seg[, ratio := 2^y]
    seg[, seg_idx := sprintf("%06d", 1L)]
    seg[, name := paste0(chrom, "_seg_", seg_idx)]
    segments_all[[chr]] <- seg
    next
  }

  cp <- cpt.mean(
    yvec,
    method="PELT",
    penalty="Manual",
    pen.value=lambda,
    class=TRUE
  )

  cpts_out <- cpts(cp)
  ends <- unique(as.integer(c(cpts_out, n)))
  ends <- ends[ends > 0 & ends <= n]
  ends <- sort(ends)
  starts <- c(1, head(ends, -1) + 1)

  seg <- data.table(
    chrom=chr,
    start=integer(length(starts)),
    end=integer(length(starts)),
    y=numeric(length(starts)),
    n_probes=integer(length(starts))
  )

  for (k in seq_along(starts)) {
    a <- starts[k]
    b <- ends[k]
    seg$start[k] <- min(sub$start[a:b])
    seg$end[k] <- max(sub$end[a:b])
    seg$y[k] <- mean(sub$y[a:b])
    seg$n_probes[k] <- (b - a + 1L)
  }

  seg <- merge_small_segments(seg, min_probes=min_probes)

  seg[, ratio := 2^y]
  seg[, seg_idx := sprintf("%06d", seq_len(.N))]
  seg[, name := paste0(chrom, "_seg_", seg_idx)]

  if (nrow(seg) >= 2) {
    bdt <- data.table(
      chrom=chr,
      pos=integer(nrow(seg)-1),
      left_seg=character(nrow(seg)-1),
      right_seg=character(nrow(seg)-1)
    )
    for (i in 1:(nrow(seg)-1)) {
      p <- as.integer((seg$end[i] + seg$start[i+1]) %/% 2)
      bdt$pos[i] <- p
      bdt$left_seg[i] <- seg$name[i]
      bdt$right_seg[i] <- seg$name[i+1]
    }
    bdt[, name := left_seg]  # boundary keyed by left segment name
    boundaries_all[[chr]] <- bdt
  }

  segments_all[[chr]] <- seg
}

seg_dt <- rbindlist(segments_all, use.names=TRUE, fill=TRUE)
setorder(seg_dt, chrom, start, end)

seg_tsv <- file.path(out_prefix, "segments.tsv.gz")
fwrite(seg_dt, seg_tsv, sep="\t")

seg_bed <- seg_dt[, .(
  chrom,
  start,
  end,
  name,
  bed_score = 0L,
  strand = ".",
  y,
  ratio,
  depth_mean = NA_real_,
  n_probes
)]
seg_bed_path <- file.path(out_prefix, "segments.bed.gz")
fwrite(seg_bed, seg_bed_path, sep="\t", col.names=FALSE)

b_all <- if (length(boundaries_all) > 0) rbindlist(boundaries_all, use.names=TRUE, fill=TRUE) else data.table()
setorder(b_all, chrom, pos)

b_tsv <- file.path(out_prefix, "boundaries.tsv.gz")
fwrite(b_all, b_tsv, sep="\t")

b_bed <- b_all[, .(
  chrom,
  start = pmax(0L, pos - 1L),
  end = pos,
  name,
  score = 0L,
  strand = "."
)]
b_bed_path <- file.path(out_prefix, "boundaries.bed.gz")
fwrite(b_bed, b_bed_path, sep="\t", col.names=FALSE)

message(sprintf("[segment] wrote: %s", seg_tsv))
message(sprintf("[segment] wrote: %s", seg_bed_path))
message(sprintf("[segment] wrote: %s", b_tsv))
message(sprintf("[segment] wrote: %s", b_bed_path))

