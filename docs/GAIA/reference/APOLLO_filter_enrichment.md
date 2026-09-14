# Filter enrichment results by source or module

Extract specific subsets from a gost_enrichment result object.

## Usage

``` r
APOLLO_filter_enrichment(
  enrich_result,
  sources = NULL,
  modules = NULL,
  top_n = NULL
)
```

## Arguments

- enrich_result:

  A gost_enrichment object from APOLLO_enrich_gost().

- sources:

  Character vector of sources to keep (e.g., c("GO:BP", "KEGG")).
  Default: NULL (all sources).

- modules:

  Character vector of modules to keep. Default: NULL (all modules).

- top_n:

  Integer. Keep only top N terms per module per source. Default: NULL
  (all).

## Value

Filtered data.frame of enrichment results.

## Examples

``` r
if (FALSE) { # \dontrun{
# Get only KEGG results from blue module
kegg_blue <- APOLLO_filter_enrichment(results, sources = "KEGG", modules = "blue")

# Get top 10 GO:BP terms per module
top_bp <- APOLLO_filter_enrichment(results, sources = "GO:BP", top_n = 10)

} # }
```
