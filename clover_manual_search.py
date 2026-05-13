#!/usr/bin/env python3
import argparse
import gzip
import math
import os
import re
import sys
from collections import defaultdict

try:
    import pysam
except ImportError:
    pysam = None


def die(msg):
    sys.stderr.write(f"ERROR: {msg}\n")
    sys.exit(1)


def say(msg):
    sys.stderr.write(f"[clover_manual_search] {msg}\n")


def open_text(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def sample_from_primary_bedgraph(path):
    b = os.path.basename(path)
    b = re.sub(r"\.gz$", "", b)
    b = re.sub(r"\.(bedGraph|bedgraph|bdg|tsv)$", "", b)

    m = re.match(r"(.+?)\.primary(?:\.|_|$)", b)
    if m:
        return m.group(1)

    return b.split(".")[0]


def load_bed(path, require_name=False):
    out = []
    with open_text(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            chrom = f[0]
            start = int(f[1])
            end = int(f[2])
            name = f[3] if len(f) >= 4 else f"{chrom}:{start}-{end}"
            out.append((chrom, start, end, name))

    if require_name:
        seen = defaultdict(int)
        fixed = []
        for chrom, start, end, name in out:
            seen[name] += 1
            if seen[name] > 1:
                name = f"{name}__{seen[name]}"
            fixed.append((chrom, start, end, name))
        out = fixed

    return out


def load_bedgraph(path):
    data = defaultdict(list)

    with open_text(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#") or line.startswith("track"):
                continue

            f = line.rstrip("\n").split()
            if len(f) < 4:
                continue

            chrom = f[0]
            start = int(f[1])
            end = int(f[2])
            val = float(f[3])

            if end > start and math.isfinite(val):
                data[chrom].append((start, end, val))

    for chrom in data:
        data[chrom].sort()

    return data


def merge_intervals(intervals):
    if not intervals:
        return []

    intervals = sorted(intervals)
    merged = [list(intervals[0])]

    for s, e in intervals[1:]:
        if s <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])

    return [(s, e) for s, e in merged]


def subtract_exclusions_from_interval(start, end, exclusions):
    pieces = [(start, end)]

    for ex_s, ex_e in exclusions:
        new_pieces = []

        for s, e in pieces:
            if ex_e <= s or ex_s >= e:
                new_pieces.append((s, e))
            else:
                if ex_s > s:
                    new_pieces.append((s, ex_s))
                if ex_e < e:
                    new_pieces.append((ex_e, e))

        pieces = new_pieces

        if not pieces:
            break

    return pieces


def weighted_mean_bedgraph(intervals, query_start=None, query_end=None, exclusions=None):
    total = 0.0
    bases = 0
    exclusions = exclusions or []

    for s, e, val in intervals:
        qs = s if query_start is None else max(s, query_start)
        qe = e if query_end is None else min(e, query_end)

        if qe <= qs:
            continue

        pieces = subtract_exclusions_from_interval(qs, qe, exclusions)

        for ps, pe in pieces:
            n = pe - ps
            total += n * val
            bases += n

    if bases == 0:
        return None, 0

    return total / bases, bases


def load_exclusions(path):
    ex = defaultdict(list)

    if not path:
        return ex

    for chrom, start, end, _ in load_bed(path):
        ex[chrom].append((start, end))

    for chrom in ex:
        ex[chrom] = merge_intervals(ex[chrom])

    return ex


def find_files(directory, pattern):
    rx = re.compile(pattern)
    hits = []

    for root, _, files in os.walk(directory):
        for f in files:
            p = os.path.join(root, f)
            if rx.search(f):
                hits.append(p)

    return sorted(hits)


def find_matching_bam(bam_dir, sample, kind):
    bams = find_files(
        bam_dir,
        rf"(?i).*{re.escape(sample)}.*{kind}.*\.bam$"
    )

    if bams:
        return bams[0]

    if kind == "splitters":
        alt = r"(?i)(split|splitter|supp)"
    else:
        alt = r"(?i)(disc|discord|discordant)"

    bams = [
        p for p in find_files(bam_dir, rf"(?i).*{re.escape(sample)}.*\.bam$")
        if re.search(alt, os.path.basename(p))
    ]

    return sorted(bams)[0] if bams else None


def find_matching_primary_bam(bam_dir, sample):
    if not bam_dir:
        return None

    bams = [
        p for p in find_files(bam_dir, rf"(?i).*{re.escape(sample)}.*\.bam$")
        if re.search(r"(?i)(primary|aligned|sorted)", os.path.basename(p))
        and not re.search(
            r"(?i)(split|splitter|supp|supplementary|disc|discord|discordant)",
            os.path.basename(p)
        )
    ]

    return sorted(bams)[0] if bams else None


def count_bam_reads(
    bam_path,
    chrom,
    start,
    end,
    count_all_alignments=False,
    drop_dup=True
):
    if pysam is None:
        die("pysam is required for BAM counting.")

    if not bam_path or not os.path.exists(bam_path):
        return "NA"

    start = max(0, int(start))
    end = max(start, int(end))

    count = 0

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        try:
            iterator = bam.fetch(chrom, start, end)
        except ValueError:
            return "NA"

        for read in iterator:
            if read.is_unmapped:
                continue
            if drop_dup and read.is_duplicate:
                continue
            if not count_all_alignments and read.is_paired and not read.is_read1:
                continue
            if read.reference_end is None:
                continue
            if read.reference_end <= start or read.reference_start >= end:
                continue

            count += 1

    return count


def total_or_na(a, b):
    if a == "NA" and b == "NA":
        return "NA"
    return (0 if a == "NA" else a) + (0 if b == "NA" else b)


def main():
    ap = argparse.ArgumentParser(
        description=(
            "Re-score user BED regions using primary bedGraphs and "
            "split/discordant BAM boundary support."
        )
    )

    ap.add_argument("--regions", required=True, help="BED file of target regions: chrom start end [name]")
    ap.add_argument("--bedgraph-dir", required=True, help="Directory containing *.primary* bedGraph files")
    ap.add_argument("--bam-dir", required=True, help="Directory containing splitters and discordant BAMs")
    ap.add_argument("--primary-bam-dir", default=None, help="Directory containing primary/aligned BAMs for raw region read counts")
    ap.add_argument("--exclude-bed", default=None, help="BED regions to exclude from chromosome baseline")
    ap.add_argument("--out", required=True, help="Output TSV")

    ap.add_argument(
        "--primary-glob-regex",
        default=r"(?i).*primary.*\.(bedGraph|bedgraph|bdg|tsv)(\.gz)?$"
    )
    ap.add_argument("--window", type=int, default=1000, help="Total boundary window size. Default: 1000")
    ap.add_argument("--cn-base", type=float, default=1.0, help="Neutral CN baseline. Default: 1.0")

    ap.add_argument("--count-all-alignments", action="store_true",
                    help="Count all alignments instead of read1-only for paired reads")
    ap.add_argument("--keep-duplicates", action="store_true",
                    help="Do not drop duplicate reads during boundary counting")

    args = ap.parse_args([x for x in sys.argv[1:] if x.strip()])

    if args.window <= 0:
        die("--window must be > 0")

    half = args.window // 2

    regions = load_bed(args.regions, require_name=True)
    exclusions = load_exclusions(args.exclude_bed)

    primary_files = find_files(args.bedgraph_dir, args.primary_glob_regex)

    if not primary_files:
        die(f"No primary bedGraph files found in {args.bedgraph_dir}")

    say(f"Loaded {len(regions)} regions")
    say(f"Found {len(primary_files)} primary bedGraph files")

    header = [
        "sample",
        "chrom",
        "start",
        "end",
        "region_name",
        "region_bp",
        "region_covered_bp",
        "chrom_baseline_covered_bp",
        "region_primary_mean",
        "chrom_primary_baseline_mean",
        "ratio",
        "log2_ratio",
        "cn",
        "primary_bam",
        "region_primary_reads",
        "left_boundary",
        "right_boundary",
        "boundary_window",
        "splitters_bam",
        "discordant_bam",
        "left_splitter_reads",
        "right_splitter_reads",
        "left_discordant_reads",
        "right_discordant_reads",
        "left_total_support_reads",
        "right_total_support_reads",
    ]

    primary_bam_search_dir = args.primary_bam_dir or args.bam_dir

    with open(args.out, "w") as out:
        out.write("\t".join(header) + "\n")

        for bg in primary_files:
            sample = sample_from_primary_bedgraph(bg)
            say(f"Processing sample {sample}")

            split_bam = find_matching_bam(args.bam_dir, sample, "splitters")
            disco_bam = find_matching_bam(args.bam_dir, sample, "discordant")
            primary_bam = find_matching_primary_bam(primary_bam_search_dir, sample)

            if not split_bam:
                say(f"WARNING: no splitters BAM found for {sample}")
            if not disco_bam:
                say(f"WARNING: no discordant BAM found for {sample}")
            if not primary_bam:
                say(f"WARNING: no primary BAM found for {sample}")

            bg_data = load_bedgraph(bg)
            baseline_cache = {}

            for chrom, start, end, name in regions:
                region_bp = end - start

                if chrom not in bg_data:
                    vals = [
                        sample,
                        chrom,
                        start,
                        end,
                        name,
                        region_bp,
                        0,
                        0,
                        "NA",
                        "NA",
                        "NA",
                        "NA",
                        "NA",
                        primary_bam or "NA",
                        "NA",
                        start,
                        end,
                        args.window,
                        split_bam or "NA",
                        disco_bam or "NA",
                        "NA",
                        "NA",
                        "NA",
                        "NA",
                        "NA",
                        "NA",
                    ]
                    out.write("\t".join(map(str, vals)) + "\n")
                    continue

                if chrom not in baseline_cache:
                    chrom_mean, chrom_cov_bp = weighted_mean_bedgraph(
                        bg_data[chrom],
                        exclusions=exclusions.get(chrom, [])
                    )
                    baseline_cache[chrom] = (chrom_mean, chrom_cov_bp)

                chrom_mean, chrom_cov_bp = baseline_cache[chrom]

                region_mean, region_cov_bp = weighted_mean_bedgraph(
                    bg_data[chrom],
                    query_start=start,
                    query_end=end,
                    exclusions=None
                )

                if region_mean is None or chrom_mean is None or chrom_mean == 0:
                    ratio = "NA"
                    log2_ratio = "NA"
                    cn = "NA"
                else:
                    ratio_val = region_mean / chrom_mean
                    ratio = ratio_val
                    log2_ratio = math.log(ratio_val, 2) if ratio_val > 0 else "NA"
                    cn = args.cn_base * ratio_val

                region_primary_reads = count_bam_reads(
                    primary_bam,
                    chrom,
                    start,
                    end,
                    count_all_alignments=args.count_all_alignments,
                    drop_dup=True
                )

                left_s = max(0, start - half)
                left_e = start + half
                right_s = max(0, end - half)
                right_e = end + half

                left_split = count_bam_reads(
                    split_bam,
                    chrom,
                    left_s,
                    left_e,
                    count_all_alignments=args.count_all_alignments,
                    drop_dup=not args.keep_duplicates
                )
                right_split = count_bam_reads(
                    split_bam,
                    chrom,
                    right_s,
                    right_e,
                    count_all_alignments=args.count_all_alignments,
                    drop_dup=not args.keep_duplicates
                )
                left_disco = count_bam_reads(
                    disco_bam,
                    chrom,
                    left_s,
                    left_e,
                    count_all_alignments=args.count_all_alignments,
                    drop_dup=not args.keep_duplicates
                )
                right_disco = count_bam_reads(
                    disco_bam,
                    chrom,
                    right_s,
                    right_e,
                    count_all_alignments=args.count_all_alignments,
                    drop_dup=not args.keep_duplicates
                )

                vals = [
                    sample,
                    chrom,
                    start,
                    end,
                    name,
                    region_bp,
                    region_cov_bp,
                    chrom_cov_bp,
                    region_mean if region_mean is not None else "NA",
                    chrom_mean if chrom_mean is not None else "NA",
                    ratio,
                    log2_ratio,
                    cn,
                    primary_bam or "NA",
                    region_primary_reads,
                    start,
                    end,
                    args.window,
                    split_bam or "NA",
                    disco_bam or "NA",
                    left_split,
                    right_split,
                    left_disco,
                    right_disco,
                    total_or_na(left_split, left_disco),
                    total_or_na(right_split, right_disco),
                ]

                out.write("\t".join(map(str, vals)) + "\n")

    say(f"Wrote {args.out}")


if __name__ == "__main__":
    main()
