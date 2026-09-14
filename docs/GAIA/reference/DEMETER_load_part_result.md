# Load a PART clustering result

Loads a PART result saved by
[`ELEUTHIA_export_timeseries_results()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ELEUTHIA_export_timeseries_results.md).
Tries the full RDS object first; falls back to reconstructing from the
cluster assignment CSV if the RDS is not present.

## Usage

``` r
DEMETER_load_part_result(exp_dir, prefix, verbose = TRUE)
```

## Arguments

- exp_dir:

  Character. Directory containing the saved PART files.

- prefix:

  Character. Filename prefix used when the result was exported (e.g.,
  the experiment name).

- verbose:

  Logical. Print loading messages. Default: TRUE.

## Value

An `artemis_part` object, or NULL if no PART files are found. Objects
reconstructed from CSV will have `dendrogram = NULL` and `data = NULL` —
sufficient for downstream plotting and gene extraction but not for
re-running tree-based operations.

## Examples

``` r
if (FALSE) { # \dontrun{
part <- DEMETER_load_part_result("results/timeseries/experiment1/",
                                  prefix = "experiment1")
} # }
```
