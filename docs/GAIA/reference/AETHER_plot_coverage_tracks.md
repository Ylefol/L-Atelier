# Plot Genome Coverage Tracks

Produces a genome browser-style stacked coverage track plot over a
specified genomic region. Each named track shows the read pileup signal
as a filled area curve. Multiple files per track are averaged into a
single signal; an optional variance ribbon can be displayed when
replicates are provided.

Two input modes are supported — provide exactly one:

- **BigWig mode** (`bigwig_files`): reads signal from pre-computed
  BigWig files. Use when BigWigs already exist (e.g. from
  [`ELEUTHIA_bed_to_bigwig()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ELEUTHIA_bed_to_bigwig.md)).

- **BED mode** (`bed_files`): computes coverage on-the-fly for the
  requested region directly from fragment-level BED files via
  [`ELEUTHIA_bed_region_coverage()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ELEUTHIA_bed_region_coverage.md).
  More efficient when only a small region is needed and BigWig files
  have not been pre-generated.

An optional peak annotation track can be added at the bottom. Additional
annotation tracks (e.g. CpG islands, methylation sites) can be added via
`annotation_tracks`.

## Usage

``` r
AETHER_plot_coverage_tracks(
  bigwig_files = NULL,
  bed_files = NULL,
  region,
  chrom_sizes = NULL,
  normalize = "CPM",
  bin_size = 10L,
  peak_bed = NULL,
  annotation_tracks = NULL,
  group_colors = NULL,
  scale_y = "free",
  show_variance = FALSE,
  variance_type = "sd",
  track_height = 1,
  annotation_height = 0.2,
  peak_color = "steelblue4",
  vline = NULL,
  vline_color = "black",
  vline_type = "dotted",
  title = NULL,
  verbose = FALSE
)
```

## Arguments

- bigwig_files:

  Track input. Two formats are accepted:

  - **Flat** (current default): a named list or named character vector
    where each element is a character vector of BigWig file paths for
    one track. All tracks share a common y-axis when
    `scale_y = "fixed"`.

  - **Grouped**: a named list of named lists/character vectors. Each
    outer element defines a group of tracks (e.g. one omics type). When
    `scale_y = "fixed"`, the y-axis is fixed *within* each group but
    allowed to differ across groups, so tracks of different dynamic
    range (e.g. ATAC vs RNA) are not forced onto the same scale.
    Example:

        list(
          ATAC = list(ATAC_WT = c("wt1.bw","wt2.bw"), ATAC_KO = "ko.bw"),
          RNA  = list(RNA_WT  = "rna_wt.bw",           RNA_KO  = "rna_ko.bw")
        )

  Mutually exclusive with `bed_files`.

- bed_files:

  Named list of character vectors, one element per track. Each element
  is a character vector of fragment-level BED file paths. Coverage is
  computed on-the-fly for the requested region only. A named character
  vector is accepted as shorthand. Mutually exclusive with
  `bigwig_files`. Requires `chrom_sizes`.

- region:

  Genomic region to display. One of:

  - Character string: `"chr1:1000000-2000000"`

  - Named character vector:
    `c(chr="chr1", start="1000000", end="2000000")`

  - A `GRanges` object (single range)

- chrom_sizes:

  Required when using `bed_files`. Exact chromosome lengths as a named
  numeric vector or
  [`APOLLO_get_chromosome_sizes()`](https://ylefol.github.io/L-Atelier/GAIA/reference/APOLLO_get_chromosome_sizes.md)
  data.frame. Ignored in BigWig mode.

- normalize:

  Character. Normalization for BED mode only: `"CPM"` (default) or
  `"raw"`. Ignored in BigWig mode (signal is already normalized in the
  file).

- bin_size:

  Integer. Bin width in bp for BED mode only. Default `10`. Ignored in
  BigWig mode.

- peak_bed:

  Optional. Peak annotation track displayed below coverage tracks.
  Either a path to a BED file or a data.frame with at least three
  columns (chr, start, end). Equivalent to adding a `"region"` type
  entry to `annotation_tracks` named `"Peaks"`.

- annotation_tracks:

  Optional named list of additional annotation tracks to display below
  coverage tracks (and below `peak_bed` if provided). Each element is a
  named list describing one track. Two types are supported:

  `type = "region"`

  :   Genomic intervals shown as filled rectangles (e.g. CpG islands,
      repeats). Required field: `data` (BED file path or data.frame with
      chr/start/end columns). Optional fields: `color` (fill colour,
      default `"steelblue4"`), `height` (relative patchwork height,
      default = `annotation_height`).

  `type = "score"`

  :   Per-site scores shown as a lollipop plot — a vertical segment from
      0 to the score, capped with a point, both coloured by score value
      (e.g. methylation fractions). Required field: `data` (BigWig file
      path). Optional fields: `color_low` (colour at score=0, default
      `"steelblue"`), `color_high` (colour at score=1, default
      `"firebrick"`), `score_limits` (numeric(2), y-axis and colour
      scale limits; auto-detected from data if `NULL`), `height`
      (relative patchwork height, default = `annotation_height * 2`).

  `type = "genes"`

  :   Collapsed IGV-style gene model track. Shows the intron backbone,
      exon blocks (UTR height), CDS blocks (taller), strand direction
      arrows, and italic gene name labels. Multiple genes are
      automatically stacked into rows to avoid overlap. Required field:
      either `txdb` (TxDb object from
      [`APOLLO_make_txdb()`](https://ylefol.github.io/L-Atelier/GAIA/reference/APOLLO_make_txdb.md)
      — fast, preferred) or `gtf_file` (path to GTF — slower, reads
      genome-wide). Optional fields: `chr_mapping` (e.g. `"T2T"`,
      applied when using `gtf_file`), `show_cds` (logical, default
      `TRUE`), `label_genes` (logical, default `TRUE`), `label_size`
      (default `2.5`), `color` (default `"black"`), `height` (default =
      `annotation_height * 3`).

  Example:

      annotation_tracks = list(
        "CpG Islands" = list(type="region", data="cpg.bed", color="darkgreen"),
        "GpC Meth"    = list(type="score",  data="meth.bw",
                             color_low="steelblue", color_high="firebrick")
      )

- group_colors:

  Optional named character vector mapping track names to colors. If
  `NULL`, colors are assigned automatically — known group names (WT, KO,
  etc.) use the AETHER default palette; others are generated.

- scale_y:

  Character. Y-axis scaling: `"free"` (default, each track autoscales
  independently) or `"fixed"` (tracks share the same y-axis maximum —
  globally when `bigwig_files` is flat, or within each group when
  `bigwig_files` is grouped).

- show_variance:

  Logical. When a track contains multiple files, display a shaded ribbon
  showing the spread across replicates. Default `FALSE`. Has no effect
  on single-file tracks.

- variance_type:

  Character. Type of variance ribbon when `show_variance = TRUE`: `"sd"`
  (default, mean ± 1 SD) or `"range"` (min to max across replicates).

- track_height:

  Numeric. Relative height unit for each coverage track. Default `1`.

- annotation_height:

  Numeric. Default relative height unit for annotation tracks
  (`peak_bed` and `annotation_tracks` entries that do not specify their
  own `height`). Default `0.2`.

- peak_color:

  Character. Fill color for peak rectangles in the `peak_bed` annotation
  track. Default `"steelblue4"`.

- vline:

  Optional numeric vector of genomic positions (bp) at which to draw a
  vertical reference line through every track and annotation panel. Each
  position is drawn as an independent line within its panel, so panels
  are not visually connected. Accepts either plain numerics or
  `"chr:position"` strings (e.g. `"chr10:74073596"`) — the chromosome is
  ignored since all panels share the same region. Default `NULL` (no
  lines).

- vline_color:

  Character. Colour of the vertical reference line(s). Default
  `"black"`.

- vline_type:

  Character. Line type of the vertical reference line(s) (any value
  accepted by `linetype`). Default `"dotted"`.

- title:

  Optional character string. Plot title placed above the first track.

- verbose:

  Logical. Print progress messages. Default `FALSE`.

## Value

A `patchwork` plot object (combination of `ggplot` panels).

## Details

When multiple files are provided for a track, per-bin mean is computed
across all files. Bins absent in a file are treated as zero before
averaging.

BigWig mode requires: `GenomicRanges`, `IRanges`, `rtracklayer`. BED
mode requires: `data.table`, `GenomicRanges`, `IRanges`, `GenomeInfoDb`.
Both modes require: `patchwork`.
