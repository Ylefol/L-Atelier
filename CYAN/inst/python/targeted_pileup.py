#!/usr/bin/env python
"""
Targeted per-base pileup at a fixed, sparse list of genomic positions.

This is a lightweight alternative to REDItools2 for the germline-variant
confirmation step (CYAN_flag_germline_variants()): that step only needs
coverage and per-base (A/C/G/T) read counts at a short list of already-known
candidate positions from the RNA side, not REDItools2's general whole-genome
editing-detection scan. REDItools2's -B/bed_file option only filters which
positions get WRITTEN to output -- it still walks every read, one base at a
time, across the entire BAM (confirmed by profiling: a 150-site test took
hours because the underlying scan was genome-wide, not site-restricted).

Each BED site is answered via pysam's indexed pileup (an O(sites) indexed
seek instead of an O(genome length) base-by-base scan), but at real
candidate-list scale (millions of known RNA-editing sites) a naive one
pileup()-call-per-site loop is itself the bottleneck: profiling a live
2.19M-site run showed ~33 sites/sec (~18h/sample), because each call pays
htslib's index-seek + BGZF-decompression + pysam object-construction
overhead regardless of how little work the single base underneath it is.
Candidate sites are not uniformly spread, though -- a real candidate list's
positions cluster tightly (median gap ~70bp, 85% of consecutive sites on
the same chromosome within 500bp of each other), so two optimizations are
applied before falling back to the original per-site semantics:

  1. Windowing: nearby sites (within --merge-gap bp of each other) are
     merged into a single region and answered with ONE pileup() call that
     covers all of them, instead of one call per site. This cut pileup()
     call volume by >10x on a real candidate list (benchmarked on this
     project's own S327 candidate BED, chr1).
  2. Multiprocessing: windows are split into --processes roughly-equal
     (by total site count) chunks and processed by independent worker
     processes, each with its own BAM/reference file handle -- embarrassingly
     parallel since windows never share reads across worker boundaries.

Per-site/per-read filtering logic is UNCHANGED from the original
implementation (see passes_read_filters()) -- windowing only changes how
many pileup() calls are made to answer the same set of requested positions,
not which reads or bases count towards each position's result. Zero-coverage
requested positions are still emitted (initialized to 0 before the window's
pileup columns are applied), same as before.

Read/base-level filtering mirrors REDItools2's own filter_read()/filter_base()
logic (see reditools.py) so germline-status calls stay consistent with how
coverage is computed elsewhere in this package:
  - read must be mapped, not QC-fail, not secondary/supplementary, not a
    PCR duplicate (SAM FLAG 0x400 -- relies on upstream duplicate-marking,
    same caveat as REDItools2)
  - mapping quality >= min_read_quality
  - read length >= min_read_length
  - if paired, must be in one of the canonical proper-pair flag
    combinations (99, 147, 83, 163)
  - must not carry an 'SA' tag (chimeric/split-read alignment)
  - base quality at the aligned position >= min_base_quality

Output is a tab-delimited file with the same core column names REDItools2
uses (Region, Position, Reference, Coverage-q<N>, BaseCount[A,C,G,T]) so it
can be parsed with the same BaseCount-splitting approach, but every
candidate site is written exactly once (including zero-coverage sites),
unlike REDItools2 which only emits positions it happens to cover. Row order
is NOT guaranteed to match the input BED's order (it is grouped by merged
window, windows sorted by position) -- downstream parsing
(.cyan_parse_targeted_pileup_output()) does not depend on row order.

Usage: python targeted_pileup.py -f <bam> -r <reference> -B <bed> -o <output>
       [-q min_read_quality] [-bq min_base_quality] [-mrl min_read_length]
       [-p processes] [--merge-gap bp] [-V]
"""

import argparse
import multiprocessing
import os
import sys
import tempfile

import pysam


def passes_read_filters(aln, min_mapq, min_read_length):
    f = aln.flag

    # Unmapped (FLAG 77 or 141)
    if f == 77 or f == 141:
        return False
    if aln.is_qcfail:
        return False
    if aln.mapping_quality < min_mapq:
        return False
    if aln.query_length < min_read_length:
        return False
    if aln.is_secondary or aln.is_supplementary:
        return False
    if aln.is_duplicate:
        return False
    # If paired, must be one of the canonical proper-pair orientations
    # (99=+-/read1, 147=+-/read2, 83=-+/read1, 163=-+/read2)
    if aln.is_paired and f not in (99, 147, 83, 163):
        return False
    if aln.has_tag("SA"):
        return False

    return True


def load_sites(bed_path):
    sites = []
    with open(bed_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2]) if len(fields) > 2 else start + 1
            sites.append((chrom, start, end))
    return sites


def merge_into_windows(sites, merge_gap):
    """Group sorted (chrom, start, end) sites into (chrom, win_start,
    win_end, [positions]) windows, merging consecutive same-chromosome
    sites that are within merge_gap bp of each other. Expects `sites`
    already sorted by (chrom, start)."""
    windows = []
    cur_chrom = None
    cur_start = None
    cur_end = None
    cur_positions = []

    for chrom, start, end in sites:
        if cur_chrom == chrom and start - cur_end <= merge_gap:
            cur_end = max(cur_end, end)
            cur_positions.append(start)
        else:
            if cur_chrom is not None:
                windows.append((cur_chrom, cur_start, cur_end, cur_positions))
            cur_chrom, cur_start, cur_end = chrom, start, end
            cur_positions = [start]

    if cur_chrom is not None:
        windows.append((cur_chrom, cur_start, cur_end, cur_positions))

    return windows


def split_windows(windows, n_chunks):
    """Split windows into n_chunks contiguous groups balanced by total
    site count (not window count) -- windows vary widely in how many
    requested positions they carry."""
    total_sites = sum(len(w[3]) for w in windows)
    target = max(1, total_sites / n_chunks)

    chunks = []
    current = []
    current_sites = 0
    for w in windows:
        current.append(w)
        current_sites += len(w[3])
        if current_sites >= target and len(chunks) < n_chunks - 1:
            chunks.append(current)
            current = []
            current_sites = 0
    if current:
        chunks.append(current)
    return chunks


_MITO_SYNONYMS = {"M": "MT", "MT": "M"}


def chrom_candidates(chrom):
    """Names to try, in order, for a chromosome that may not match a BAM/
    FASTA's own naming convention: the name as given, its chr-prefix swapped
    form, and (since this project's liftOver output uses 'M' for the
    mitochondrial contig while this reference's BAM/FASTA use 'MT' --
    confirmed the cause of a worker crash on 2026-11-04, see CYAN_STATUS.txt)
    the M/MT synonym of both of those, with or without a chr prefix."""
    candidates = [chrom]
    candidates.append(chrom[3:] if chrom.startswith("chr") else "chr" + chrom)

    for c in list(candidates):
        base = c[3:] if c.startswith("chr") else c
        if base in _MITO_SYNONYMS:
            swapped = _MITO_SYNONYMS[base]
            candidates.append(swapped)
            candidates.append("chr" + swapped)

    seen = set()
    ordered = []
    for c in candidates:
        if c not in seen:
            seen.add(c)
            ordered.append(c)
    return ordered


def open_pileup(samfile, chrom, start, end):
    last_err = None
    for candidate in chrom_candidates(chrom):
        try:
            return samfile.pileup(
                candidate, start, end, truncate=True,
                min_base_quality=0, min_mapping_quality=0,
                ignore_overlaps=False, stepper="all",
            )
        except ValueError as e:
            last_err = e
            continue
    raise ValueError(
        "Could not find chromosome '{}' in the BAM under any of these names: {}. "
        "Last error: {}".format(chrom, chrom_candidates(chrom), last_err)
    )


def ref_base_at(reference_reader, chrom, pos0):
    if reference_reader is None:
        return "N"
    for name in chrom_candidates(chrom):
        try:
            return reference_reader.fetch(name, pos0, pos0 + 1).upper() or "N"
        except (KeyError, ValueError):
            continue
    return "N"


def process_chunk(args):
    (chunk, bam_path, reference_path, out_path, min_read_quality,
     min_base_quality, min_read_length, coverage_label, worker_id,
     verbose, progress_every) = args

    samfile = pysam.AlignmentFile(bam_path, "rb")
    reference_reader = pysam.FastaFile(reference_path) if reference_path else None

    n_done = 0
    total = sum(len(w[3]) for w in chunk)

    with open(out_path, "w") as out:
        for chrom, win_start, win_end, positions in chunk:
            counts_by_pos = {p: {"A": 0, "C": 0, "G": 0, "T": 0} for p in positions}

            for col in open_pileup(samfile, chrom, win_start, win_end):
                pos0 = col.reference_pos
                if pos0 not in counts_by_pos:
                    continue
                counts = counts_by_pos[pos0]
                for pread in col.pileups:
                    if pread.is_del or pread.is_refskip:
                        continue
                    aln = pread.alignment
                    if not passes_read_filters(aln, min_read_quality, min_read_length):
                        continue
                    qpos = pread.query_position
                    if qpos is None:
                        continue
                    if aln.query_qualities[qpos] < min_base_quality:
                        continue
                    base = aln.query_sequence[qpos].upper()
                    if base in counts:
                        counts[base] += 1

            for pos0 in positions:
                counts = counts_by_pos[pos0]
                coverage = sum(counts.values())
                ref = ref_base_at(reference_reader, chrom, pos0)
                out.write("{}\t{}\t{}\t{}\t[{}, {}, {}, {}]\n".format(
                    chrom, pos0 + 1, ref, coverage,
                    counts["A"], counts["C"], counts["G"], counts["T"],
                ))

            n_done += len(positions)
            if verbose and progress_every and n_done % progress_every < len(positions):
                sys.stderr.write("[TARGETED_PILEUP] worker {}: {}/{} sites done\n".format(
                    worker_id, n_done, total))

    samfile.close()
    if reference_reader is not None:
        reference_reader.close()

    return out_path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--bam", required=True)
    parser.add_argument("-r", "--reference", default=None)
    parser.add_argument("-B", "--bed", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("-q", "--min-read-quality", type=int, default=20)
    parser.add_argument("-bq", "--min-base-quality", type=int, default=30)
    parser.add_argument("-mrl", "--min-read-length", type=int, default=30)
    parser.add_argument("-p", "--processes", type=int, default=1,
                         help="Number of worker processes (default 1, sequential).")
    parser.add_argument("--merge-gap", type=int, default=1000,
                         help="Merge candidate sites within this many bp of each "
                              "other into one pileup() call (default 1000).")
    parser.add_argument("-V", "--verbose", action="store_true")
    args = parser.parse_args()

    sites = sorted(load_sites(args.bed))
    if args.verbose:
        sys.stderr.write("[TARGETED_PILEUP] {} candidate sites\n".format(len(sites)))

    windows = merge_into_windows(sites, args.merge_gap)
    if args.verbose:
        sys.stderr.write("[TARGETED_PILEUP] merged into {} windows (gap<={}bp)\n".format(
            len(windows), args.merge_gap))

    coverage_label = "Coverage-q{}".format(args.min_base_quality)
    n_proc = max(1, args.processes)
    chunks = split_windows(windows, n_proc)
    n_proc = len(chunks)  # don't spin up more workers than there are chunks

    tmp_dir = tempfile.mkdtemp(prefix="cyan_targeted_pileup_")
    chunk_paths = [os.path.join(tmp_dir, "chunk_{}.txt".format(i)) for i in range(n_proc)]

    total_sites = len(sites)
    progress_every = max(1, total_sites // (n_proc * 20))

    job_args = [
        (chunks[i], args.bam, args.reference, chunk_paths[i],
         args.min_read_quality, args.min_base_quality, args.min_read_length,
         coverage_label, i, args.verbose, progress_every)
        for i in range(n_proc)
    ]

    if n_proc == 1:
        process_chunk(job_args[0])
    else:
        if args.verbose:
            sys.stderr.write("[TARGETED_PILEUP] running {} worker processes\n".format(n_proc))
        with multiprocessing.Pool(processes=n_proc) as pool:
            pool.map(process_chunk, job_args)

    with open(args.output, "w") as out:
        out.write("Region\tPosition\tReference\t{}\tBaseCount[A,C,G,T]\n".format(coverage_label))
        for p in chunk_paths:
            with open(p) as chunk_f:
                for line in chunk_f:
                    out.write(line)

    for p in chunk_paths:
        os.remove(p)
    os.rmdir(tmp_dir)

    sys.stderr.write("[TARGETED_PILEUP] Done. {} sites written to {}\n".format(
        total_sites, args.output))


if __name__ == "__main__":
    main()
