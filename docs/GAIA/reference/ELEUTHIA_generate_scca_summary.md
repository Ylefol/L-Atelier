# Generate sCCA Analysis Summary (Plain Text)

Internal function to generate a plain text summary of sCCA
cross-validation analysis.

## Usage

``` r
ELEUTHIA_generate_scca_summary(
  coarse_cv_results,
  fine_cv_results,
  pilot_results,
  method,
  X,
  metadata,
  saved_files
)
```

## Arguments

- coarse_cv_results:

  List from MINERVA_cv_scca

- fine_cv_results:

  List from MINERVA_cv_scca

- pilot_results:

  List from MINERVA_pilot_compare_methods (optional)

- method:

  Character string

- X:

  List of datasets

- metadata:

  List of additional metadata

- saved_files:

  List of saved file paths

## Value

Character vector with plain text content (one element per line)
