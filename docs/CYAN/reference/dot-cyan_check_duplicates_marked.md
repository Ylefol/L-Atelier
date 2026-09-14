# Check a BAM for duplicate-flagged reads (SAM FLAG 0x400)

REDItools2's own duplicate filter (`read.is_duplicate`, hardcoded in
`reditools.py`) only excludes reads that already carry this flag – it
does not mark duplicates itself. On a BAM that was never run through a
duplicate-marking tool (e.g. `samtools markdup`, Picard
`MarkDuplicates`), the flag is never set and that filter is a silent
no-op, letting PCR-duplicate reads inflate editing-frequency estimates
uncaught. This scans for at least one flagged read as a cheap upstream
sanity check, exiting immediately once one is found.

## Usage

``` r
.cyan_check_duplicates_marked(bam, max_reads = 5e+06, verbose = TRUE)
```

## Arguments

- bam:

  Character. Path to a sorted, indexed BAM file.

- max_reads:

  Integer. Maximum number of reads to scan before concluding none are
  flagged (default 5,000,000). Only reached for a genuinely
  duplicate-free BAM – any real duplicate-marked BAM exits on the first
  flagged read.

- verbose:

  Logical. Print `[CYAN]` progress messages.

## Value

Logical. TRUE if at least one duplicate-flagged read was found.
