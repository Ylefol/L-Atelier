# Convert a processed BAM to a fragment-level BED file

Produces a 6-column fragment BED from a paired-end BAM. Optionally
applies the standard Tn5 insertion-site correction (+4 bp on the +
strand, −5 bp on the − strand) to shift each fragment toward the actual
cut site.

## Usage

``` r
HORIZON_bam_to_bed(
  sample_sheet,
  sample_id,
  bam_file = NULL,
  tag = NULL,
  shift_reads = TRUE,
  threads = 4L,
  remove_tmp = TRUE,
  force = FALSE
)
```

## Arguments

- sample_sheet:

  Validated sample sheet data.frame. Required when `bam_file` is `NULL`.

- sample_id:

  Character. Sample ID to process.

- bam_file:

  Character or `NULL`. Path to a coordinate-sorted, indexed BAM file.
  When supplied, `.latest_bam()` is bypassed and this BAM is used
  directly. Output is written to `dirname(bam_file)`. Default `NULL`.

- tag:

  Character or `NULL`. Optional suffix inserted before `_fragments.bed`
  in the output filename, producing `<sample_id>_<tag>_fragments.bed`.
  Use this to avoid overwriting when calling the function twice for the
  same sample at different pipeline stages (e.g.
  `tag = "quantification"` for the pre-downsampled BED used by DESeq2,
  leaving the untagged file for peak calling). Default `NULL` (no tag —
  output is `<sample_id>_fragments.bed`).

- shift_reads:

  Logical. Apply Tn5 insertion-site correction (+4/−5 bp). Set `TRUE`
  for ATAC-seq, `FALSE` for all other assays. Default `TRUE`.

- threads:

  Integer. Threads for `samtools sort`. Default 4.

- remove_tmp:

  Logical. Remove intermediate name-sorted BAM and BEDPE after
  conversion. Default `TRUE`.

- force:

  Logical. If `FALSE` (default), skip conversion when the output BED
  already exists. Set `TRUE` to regenerate.

## Value

Character. Path to the fragment BED file, invisibly.

## Details

The output BED is directly compatible with `ELEUTHIA_load_bed()` and
`ELEUTHIA_merge_fragments()` in GAIA.

Internally:

1.  Name-sorts the BAM (`samtools sort -n`).

2.  Converts to BEDPE (`bedtools bamtobed -bedpe`).

3.  Streams the BEDPE through `awk` to build fragment coordinates and
    optionally apply the Tn5 shift, then sorts with `sort`.

**When to use `shift_reads`:**

- ATAC-seq — `TRUE`: the Tn5 transposase cuts and ligates adapters, so
  the read start is offset from the actual insertion site by 4/5 bp.

- ChIP-seq / CUT&RUN / CUT&TAG — `FALSE`: no insertion-site correction
  is needed; fragment coordinates are used as-is.
