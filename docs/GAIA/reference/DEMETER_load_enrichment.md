# Load saved enrichment results

Loads a `gost_enrichment` object (from
[`APOLLO_enrich_gost()`](https://ylefol.github.io/L-Atelier/GAIA/reference/APOLLO_enrich_gost.md))
that was saved by the user via
[`saveRDS()`](https://rdrr.io/r/base/readRDS.html), or reconstructs a
partial object from a flat CSV export. RDS is the recommended format: it
preserves the complete object including per-module gost objects and
metadata. CSV reconstruction is provided as a fallback — only the
`$combined` data.frame is restored; per-module gost objects and
`$metadata` are unavailable.

## Usage

``` r
DEMETER_load_enrichment(file_path, verbose = TRUE)
```

## Arguments

- file_path:

  Character. Path to an `.rds` or `.csv` file. For CSV, the file must
  contain at minimum the columns: `term_id`, `term_name`, `source`,
  `p_value`, `module` (i.e., the `$combined` output from
  [`APOLLO_enrich_gost()`](https://ylefol.github.io/L-Atelier/GAIA/reference/APOLLO_enrich_gost.md)).

- verbose:

  Logical. Print loading summary. Default: `TRUE`.

## Value

A `gost_enrichment` S3 object (list) with elements:

- results:

  Named list of raw gprofiler2 gost objects (RDS only; empty list when
  loading from CSV)

- combined:

  Combined data.frame of all results across all modules

- summary:

  Summary data.frame of term counts (RDS only; NULL from CSV)

- metadata:

  List of analysis parameters including organism and sources (RDS only;
  partially inferred from CSV)

## See also

[`APOLLO_enrich_gost`](https://ylefol.github.io/L-Atelier/GAIA/reference/APOLLO_enrich_gost.md),
[`AETHER_plot_go_treemap`](https://ylefol.github.io/L-Atelier/GAIA/reference/AETHER_plot_go_treemap.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Save enrichment result after running the analysis
gost_res <- APOLLO_enrich_gost(gene_lists, sources = c("GO:BP", "KEGG"))
saveRDS(gost_res, "results/enrichment/gost_results.rds")

# Reload in a fresh session
gost_res <- DEMETER_load_enrichment("results/enrichment/gost_results.rds")

# Plot GO treemap from reloaded object
AETHER_plot_go_treemap(gost_res, ont = "BP")
} # }
```
