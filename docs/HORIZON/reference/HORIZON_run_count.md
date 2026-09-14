# Count reads per feature using featureCounts

Wraps
[`featureCounts`](https://rdrr.io/pkg/Rsubread/man/featureCounts.html)
to quantify reads mapping to genomic features for a single sample.
Strandedness is mapped to featureCounts integer codes (0 = unstranded, 1
= forward, 2 = reverse).

## Usage

``` r
HORIZON_run_count(
  sample_sheet = NULL,
  sample_id = NULL,
  bam_file = NULL,
  output_dir = NULL,
  paired_end = NULL,
  strandedness = NULL,
  annotation,
  feature_type = "exon",
  attribute_type = "gene_id",
  threads = 4,
  min_mapping_quality = 10,
  ...
)
```

## Arguments

- sample_sheet:

  Validated sample sheet data.frame from
  [`HORIZON_validate_sample_sheet`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_validate_sample_sheet.md).
  Required in sample-sheet mode; must be `NULL` when `bam_file` is
  supplied.

- sample_id:

  Character. Sample ID to process (sample-sheet mode). Optional in
  direct BAM mode — inferred from the filename if omitted.

- bam_file:

  Character or `NULL`. Path to a sorted BAM file. When supplied,
  `sample_sheet` must be `NULL`. Default `NULL`.

- output_dir:

  Character or `NULL`. Output directory for direct BAM mode. Ignored in
  sample-sheet mode. Defaults to the directory containing `bam_file`.

- paired_end:

  Logical or `NULL`. Whether the BAM is paired-end. Required in direct
  BAM mode; read from the sample sheet in sample-sheet mode.

- strandedness:

  Character or `NULL`. One of `"unstranded"`, `"forward"`, `"reverse"`.
  Required in direct BAM mode; read from the sample sheet in
  sample-sheet mode.

- annotation:

  Character. Path to a GTF/GFF annotation file. Passed to `annot.ext` in
  [`featureCounts`](https://rdrr.io/pkg/Rsubread/man/featureCounts.html).

- feature_type:

  Character. Feature type in the GTF to count reads against. Default
  `"exon"`.

- attribute_type:

  Character. GTF attribute used to group features into meta-features
  (genes). Default `"gene_id"`.

- threads:

  Integer. Number of threads. Default 4.

- min_mapping_quality:

  Integer. Minimum mapping quality score for a read to be counted.
  Default 10.

- ...:

  Additional arguments passed to
  [`featureCounts`](https://rdrr.io/pkg/Rsubread/man/featureCounts.html)
  (e.g., `allowMultiOverlap`, `countMultiMappingReads`,
  `requireBothEndsMapped`, `minFragLength`, `maxFragLength`).

## Value

Character. Path to the counts text file, invisibly.

## Details

**Two usage modes:**

1.  **Sample-sheet mode** (default): provide `sample_sheet` and
    `sample_id`. `paired_end` and `strandedness` are read from the
    sample sheet row, and the sorted BAM is located at the standard
    HORIZON path
    (`<output_dir>/<sample_id>/aligned/<sample_id>_sorted.bam`). Output
    is written to `<output_dir>/<sample_id>/counts/`.

2.  **Direct BAM mode**: provide `bam_file` (path to a sorted BAM), plus
    `paired_end` and `strandedness` explicitly, since neither can be
    reliably inferred without a sample sheet row. `sample_id` is
    inferred from the filename if omitted. Output goes to `output_dir`
    (defaults to the directory containing the BAM).

Two files are written to the resolved output directory:

- `<sample_id>_counts.txt` — tab-delimited table of gene ID and raw
  count (used by
  [`HORIZON_aggregate_counts`](https://ylefol.github.io/L-Atelier/HORIZON/reference/HORIZON_aggregate_counts.md))

- `<sample_id>_featurecounts.rds` — full featureCounts result object
  including alignment statistics and annotation
