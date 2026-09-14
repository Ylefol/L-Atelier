# Load decoupleR activity results

Loads previously saved decoupleR results from RDS files. Use this to
reload activity inference results without re-running the analysis.

## Usage

``` r
DEMETER_load_activity(file_path, verbose = TRUE)
```

## Arguments

- file_path:

  Path to the RDS file containing decoupler results.

- verbose:

  Logical. Print progress. Default: TRUE.

## Value

The loaded decoupler_result or decoupler_comparison object.

## Examples

``` r
if (FALSE) { # \dontrun{
# Save results (in your analysis script)
saveRDS(tf_result, "results/tf_activity.rds")

# Load later
tf_result <- DEMETER_load_activity("results/tf_activity.rds")

} # }
```
