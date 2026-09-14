# Plot Signal Profile Heatmap (scale-regions layout)

Produces a deepTools `computeMatrix scale-regions`-style figure. Each
panel shows a mean signal profile curve above a row-sorted heatmap. The
x-axis spans from `window` bp upstream of the TSS through the scaled
gene body to `window` bp downstream of the TES. The gene body is scaled
to `body_bins` equal bins regardless of gene length, so genes of
different sizes are directly comparable. Rows are sorted by descending
total signal (highest at top). Multiple samples are displayed
side-by-side in one composite figure.

## Usage

``` r
AETHER_plot_profile_heatmap(
  bigwig_files,
  txdb,
  genes = NULL,
  peaks = NULL,
  window = 3000L,
  bin_size = 50L,
  body_bins = 100L,
  colors = c("#A50026", "#D73027", "#F46D43", "#FDAE61", "#FEE090", "#FFFFBF", "#E0F3F8",
    "#ABD9E9", "#74ADD1", "#4575B4", "#313695"),
  profile_color = "#2166AC",
  raster_interpolate = TRUE,
  cap_quantile = 0.99,
  cap_val = NULL,
  profile_ylim = NULL,
  max_genes = 5000L,
  profile_height = 1,
  heatmap_height = 5,
  title = NULL,
  verbose = TRUE
)
```

## Arguments

- bigwig_files:

  Named list of BigWig file paths. Each element is a character vector.
  If multiple files are supplied they are averaged into one track.
  Element names become the panel labels. Example:
  `list(WT = c("wt1.bw","wt2.bw"), KO = "ko.bw")`

- txdb:

  TxDb object (from
  [`APOLLO_make_txdb()`](https://ylefol.github.io/L-Atelier/GAIA/reference/APOLLO_make_txdb.md)).
  Used to extract TSS and TES positions for all (or selected) genes.

- genes:

  Optional character vector of gene identifiers to include. If `NULL`
  (default) all genes in `txdb` are used. Matched against raw TxDb gene
  IDs and NCBI "gene-" stripped forms.

- peaks:

  Optional `GRanges` object or path to a BED/narrowPeak file of called
  peaks (e.g. from ATAC or ChIP). When supplied, only genes whose
  extended region (body ± `window`) overlaps at least one peak are
  retained. This is more principled than the `max_genes` heuristic:
  rather than guessing the top-N genes by TSS signal, you provide the
  exact loci of interest. `max_genes` still acts as a hard cap on the
  result.

- window:

  Integer. Base pairs upstream of TSS and downstream of TES to include
  as flanking regions. Default `3000`.

- bin_size:

  Integer. Bin width in bp for the flanking regions. The number of flank
  bins on each side = `ceiling(window / bin_size)`. Default `50`.

- body_bins:

  Integer. Number of bins the gene body is scaled to, regardless of gene
  length. Default `100`.

- colors:

  Character vector of \\\ge\\ 2 colours, ordered from zero/low signal to
  maximum signal, passed to
  [`ggplot2::scale_fill_gradientn()`](https://ggplot2.tidyverse.org/reference/scale_gradient.html).
  Default is an 11-stop red-yellow-blue diverging ramp (matching the
  ColorBrewer/deepTools `"RdYlBu"` palette) that gives low signal a
  distinct red, mid-range signal a pale yellow, and high signal a
  distinct dark blue — reproducing the classic deepTools `plotHeatmap`
  look. Supply exactly 2 colours for a plain low/high linear gradient
  (the previous behaviour).

- profile_color:

  Character. Line and fill colour for the profile curve. Default
  `"#2166AC"` (blue).

- raster_interpolate:

  Logical. Passed to
  [`ggplot2::geom_raster()`](https://ggplot2.tidyverse.org/reference/geom_tile.html)'s
  `interpolate` argument. `TRUE` (default) bilinearly blends adjacent
  cells, which avoids an aliased/speckled look when many gene rows are
  downsampled into a shorter rendered image (typical for `max_genes` in
  the thousands). Set `FALSE` for crisp, unblended per-bin/per-gene
  pixels.

- cap_quantile:

  Numeric in (0,1). Signal values above this quantile are capped before
  display (shared cap across all panels). Default `0.99`.

- cap_val:

  Optional numeric. Fixed colour-scale cap value shared across calls
  (e.g. a mark-level global cap). When `NULL` (default) the cap is
  computed per-call from `cap_quantile`.

- profile_ylim:

  Optional numeric. Fixed y-axis upper limit for the mean signal
  profile, shared across calls. When `NULL` (default) the limit is
  computed per-call from the data.

- max_genes:

  Integer or `NULL`. Hard cap on the number of genes displayed. When
  `peaks` is supplied this limits the peak-filtered set; when `peaks` is
  `NULL` the TSS-signal pre-filter keeps the top `max_genes` genes. Set
  `NULL` to disable. Default `5000`.

- profile_height:

  Numeric. Relative height of the profile panel. Default `1`.

- heatmap_height:

  Numeric. Relative height of the heatmap panel. Default `5`.

- title:

  Optional character string. Overall figure title.

- verbose:

  Logical. Print progress messages. Default `TRUE`.

## Value

A `patchwork` composite plot.

## Details

The gene body is scaled independently per gene: a 500 bp gene and a 50
kb gene both contribute `body_bins` columns to the matrix. The flanking
regions use a fixed `bin_size`. Minus-strand genes are flipped so that
upstream is always on the left. TSS and TES positions are marked with
dashed vertical lines on both the profile and heatmap panels.

Requires: `GenomicRanges`, `IRanges`, `S4Vectors`, `rtracklayer`,
`GenomicFeatures`, `EnrichedHeatmap`, `patchwork`.
