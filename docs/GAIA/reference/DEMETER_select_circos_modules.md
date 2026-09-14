# Select WGCNA modules for circos display

Identifies modules worth displaying in the circos plot based on multiple
criteria: trait correlation significance and strength, module size,
overlap with PART clusters, and (optionally) enrichment results. Returns
module names suitable for passing to
[`DEMETER_filter_circos_modules()`](https://ylefol.github.io/L-Atelier/GAIA/reference/DEMETER_filter_circos_modules.md).

## Usage

``` r
DEMETER_select_circos_modules(
  circos_data,
  enrichment = NULL,
  min_cor = 0.3,
  require_sig_cor = TRUE,
  min_size = 30,
  min_overlap = 1,
  exclude_grey = TRUE,
  verbose = TRUE
)
```

## Arguments

- circos_data:

  A `part_wgcna_circos` object.

- enrichment:

  Optional. Enrichment results from
  [`APOLLO_enrich_gost()`](https://ylefol.github.io/L-Atelier/GAIA/reference/APOLLO_enrich_gost.md)
  or similar. Must have a `module` column. Modules with at least one
  significant term pass the filter.

- min_cor:

  Numeric. Minimum absolute trait correlation with at least one trait.
  Default: 0.3. Set to 0 to disable.

- require_sig_cor:

  Logical. Require at least one significant trait correlation (based on
  the \_sig columns in module_df). Default: TRUE. Ignored if no trait
  data is available.

- min_size:

  Integer. Minimum module size (gene count). Default: 30.

- min_overlap:

  Integer. Minimum number of genes overlapping with any PART cluster
  (from the association matrix). Default: 1.

- exclude_grey:

  Logical. Exclude the grey (unassigned) module. Default: TRUE.

- verbose:

  Logical. Print filtering summary. Default: TRUE.

## Value

Character vector of module names passing all criteria.

## Details

All enabled criteria are combined with AND logic: a module must pass
every applicable filter. Criteria are skipped when their data is
unavailable (e.g., enrichment filtering is skipped if
`enrichment = NULL`).

## Examples

``` r
if (FALSE) { # \dontrun{
selected <- DEMETER_select_circos_modules(circos_data, min_cor = 0.4)
circos_data <- DEMETER_filter_circos_modules(circos_data, selected)

} # }
```
