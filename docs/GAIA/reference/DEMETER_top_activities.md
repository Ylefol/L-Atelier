# Select top N sources from activity results

Filters activity results to keep only the top N sources (TFs/pathways)
based on variance, mean absolute activity, or other criteria.

## Usage

``` r
DEMETER_top_activities(result, n = 20, by = "variance", verbose = TRUE)
```

## Arguments

- result:

  A decoupler_result object, or a matrix of activities (sources as rows,
  samples as columns).

- n:

  Integer. Number of top sources to keep. Default: 20.

- by:

  Character. Criterion for selecting top sources:

  - "variance": highest variance across samples (most variable)

  - "activity": highest mean absolute activity (strongest signals)

  - "max": highest maximum absolute activity Default: "variance".

- verbose:

  Logical. Print selection info. Default: TRUE.

## Value

If input was decoupler_result, returns filtered decoupler_result. If
input was matrix, returns filtered matrix.

## Details

Use this to reduce the number of TFs/pathways before:

- Plotting heatmaps (too many rows are unreadable

- Using as traits in WGCNA module-trait correlation

- Any analysis where you want to focus on the most informative sources

## Examples

``` r
if (FALSE) { # \dontrun{
# Filter to top 20 most variable TFs
tf_top <- DEMETER_top_activities(tf_result, n = 20, by = "variance")

# Use with WGCNA
tf_trait_cor <- ARTEMIS_wgcna_module_traits(
  modules = modules,
  traits = t(tf_top$activities)
)

# Filter to top 30 by activity strength
tf_top <- DEMETER_top_activities(tf_result, n = 30, by = "activity")

} # }
```
