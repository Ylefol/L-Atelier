# Export ORA results to files

Writes CSV tables, an RDS object, and optional fold-enrichment dotplots
for a `keraunos_ora` or `keraunos_ora_multi` object.

## Usage

``` r
TALARIA_export_ora(
  result,
  output_dir,
  prefix = "ora",
  plots = TRUE,
  top_n = 20L,
  plot_width = 8,
  plot_height = 7,
  verbose = TRUE
)
```

## Arguments

- result:

  A `keraunos_ora` or `keraunos_ora_multi` object.

- output_dir:

  Character. Directory to write into (created if absent).

- prefix:

  Character. Filename prefix. Default `"ora"`.

- plots:

  Logical. Generate dotplots. Default `TRUE`.

- top_n:

  Integer. Total terms shown per dotplot, distributed across sources.
  Default `20L`.

- plot_width:

  Numeric. PNG width in inches. Default `8`.

- plot_height:

  Numeric. PNG height in inches. Default `7`.

- verbose:

  Logical. Print a summary. Default `TRUE`.

## Value

Invisibly, a character vector of file paths written.

## Details

For `keraunos_ora_multi` each cluster gets its own subdirectory. Up- and
down-regulated directions are combined into a single results CSV
(distinguished by a `direction` column) and a single dotplot image with
the "up" and "down" panels shown side by side, each keeping its own
fold-enrichment axis and significance colour gradient (red for up, blue
for down) rather than a shared scale.
