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

Here, each BED site is looked up directly via pysam's indexed pileup
(samfile.pileup(contig, start, end, truncate=True)), which seeks straight to
that interval instead of scanning everything in between -- the same
has-coverage-at-this-exact-site question, answered in O(sites) indexed seeks
instead of O(genome length) base-by-base iteration.

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
unlike REDItools2 which only emits positions it happens to cover.

Usage: python targeted_pileup.py -f <bam> -r <reference> -B <bed> -o <output>
       [-q min_read_quality] [-bq min_base_quality] [-mrl min_read_length] [-V]
"""

import argparse
import sys

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


def open_pileup(samfile, chrom, start, end):
    try:
        return samfile.pileup(
            chrom, start, end, truncate=True,
            min_base_quality=0, min_mapping_quality=0,
            ignore_overlaps=False, stepper="all",
        )
    except ValueError:
        alt_chrom = chrom[3:] if chrom.startswith("chr") else "chr" + chrom
        return samfile.pileup(
            alt_chrom, start, end, truncate=True,
            min_base_quality=0, min_mapping_quality=0,
            ignore_overlaps=False, stepper="all",
        )


def ref_base_at(reference_reader, chrom, pos0):
    if reference_reader is None:
        return "N"
    for name in (chrom, chrom[3:] if chrom.startswith("chr") else "chr" + chrom):
        try:
            return reference_reader.fetch(name, pos0, pos0 + 1).upper() or "N"
        except (KeyError, ValueError):
            continue
    return "N"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--bam", required=True)
    parser.add_argument("-r", "--reference", default=None)
    parser.add_argument("-B", "--bed", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("-q", "--min-read-quality", type=int, default=20)
    parser.add_argument("-bq", "--min-base-quality", type=int, default=30)
    parser.add_argument("-mrl", "--min-read-length", type=int, default=30)
    parser.add_argument("-V", "--verbose", action="store_true")
    args = parser.parse_args()

    sites = load_sites(args.bed)
    if args.verbose:
        sys.stderr.write("[TARGETED_PILEUP] {} candidate sites\n".format(len(sites)))

    samfile = pysam.AlignmentFile(args.bam, "rb")
    reference_reader = pysam.FastaFile(args.reference) if args.reference else None

    coverage_label = "Coverage-q{}".format(args.min_base_quality)

    with open(args.output, "w") as out:
        out.write("Region\tPosition\tReference\t{}\tBaseCount[A,C,G,T]\n".format(coverage_label))

        n_done = 0
        for chrom, start, end in sites:
            for pos0 in range(start, end):
                counts = {"A": 0, "C": 0, "G": 0, "T": 0}

                for col in open_pileup(samfile, chrom, pos0, pos0 + 1):
                    if col.reference_pos != pos0:
                        continue
                    for pread in col.pileups:
                        if pread.is_del or pread.is_refskip:
                            continue
                        aln = pread.alignment
                        if not passes_read_filters(aln, args.min_read_quality, args.min_read_length):
                            continue
                        qpos = pread.query_position
                        if qpos is None:
                            continue
                        if aln.query_qualities[qpos] < args.min_base_quality:
                            continue
                        base = aln.query_sequence[qpos].upper()
                        if base in counts:
                            counts[base] += 1

                coverage = sum(counts.values())
                ref = ref_base_at(reference_reader, chrom, pos0)
                out.write("{}\t{}\t{}\t{}\t[{}, {}, {}, {}]\n".format(
                    chrom, pos0 + 1, ref, coverage,
                    counts["A"], counts["C"], counts["G"], counts["T"],
                ))

                n_done += 1
                if args.verbose and n_done % 5000 == 0:
                    sys.stderr.write("[TARGETED_PILEUP] {}/{} sites done\n".format(n_done, len(sites)))

    samfile.close()
    if reference_reader is not None:
        reference_reader.close()

    sys.stderr.write("[TARGETED_PILEUP] Done. {} sites written to {}\n".format(len(sites), args.output))


if __name__ == "__main__":
    main()
