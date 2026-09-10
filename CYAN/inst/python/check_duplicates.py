#!/usr/bin/env python
"""
Scan a BAM file for the SAM duplicate flag (0x400 / is_duplicate).

Exits as soon as one flagged read is found (near-instant for any BAM that
was actually duplicate-marked). Scans at most `max_reads` reads before
concluding none are flagged, so a genuinely duplicate-free BAM has a bounded
worst-case cost.

Usage: python check_duplicates.py <bam_path> <max_reads>
Prints "TRUE" (a duplicate-flagged read was found) or "FALSE" (none found in
the scanned reads) as the last line of stdout.
"""

import sys
import pysam


def main():
    bam_path = sys.argv[1]
    max_reads = int(sys.argv[2])

    bamfile = pysam.AlignmentFile(bam_path, "rb")
    found = False
    n = 0
    for read in bamfile:
        n += 1
        if read.is_duplicate:
            found = True
            break
        if n >= max_reads:
            break
    bamfile.close()

    print("TRUE" if found else "FALSE")


if __name__ == "__main__":
    main()
