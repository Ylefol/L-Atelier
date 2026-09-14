# Convert Columns to Factors Based on QC Report

Helper function to convert potential categorical columns (identified by
HADES_qc_mixed_data) to factors. Can also be used to convert any
specified columns.

## Usage

``` r
HADES_convert_to_factor(
  data,
  qc_report = NULL,
  columns = NULL,
  ordered = FALSE,
  exclude = NULL,
  verbose = TRUE
)
```

## Arguments

- data:

  Data frame to modify.

- qc_report:

  Optional. A hades_qc object from HADES_qc_mixed_data(). If provided,
  will use the potential_categorical list by default.

- columns:

  Character vector. Column names to convert. If NULL and qc_report is
  provided, uses all potential_categorical columns. If both qc_report
  and columns are provided, columns takes precedence.

- ordered:

  Logical or named list. If TRUE, creates ordered factors (for ordinal
  data). If a named list, specifies which columns should be ordered
  (e.g., list(severity = TRUE, group = FALSE)). Default = FALSE.

- exclude:

  Character vector. Column names to exclude from conversion even if
  flagged. Useful when you've reviewed the QC and determined some
  flagged columns are actually continuous.

- verbose:

  Logical. Print information about conversions. Default = TRUE.

## Value

The data frame with specified columns converted to factors.

## Details

This function is intended to be used after reviewing the QC report from
HADES_qc_mixed_data(). The typical workflow is:

1.  Run HADES_qc_mixed_data() and review potential_categorical columns

2.  Decide which ones are truly categorical

3.  Use this function to convert them, excluding any that are actually
    continuous

For ordinal data (e.g., severity scores 0-3), set ordered = TRUE to
create ordered factors. The levels will be sorted numerically for
numeric columns.

## Examples

``` r
if (FALSE) { # \dontrun{
# Convert all flagged potential categorical columns
data_clean <- HADES_convert_to_factor(my_data, qc_report)

# Convert specific columns only
data_clean <- HADES_convert_to_factor(my_data, columns = c("severity", "stage"))

# Exclude some flagged columns (they're actually continuous)
data_clean <- HADES_convert_to_factor(my_data, qc_report, exclude = c("age_group"))

# Create ordered factors for ordinal data
data_clean <- HADES_convert_to_factor(my_data, qc_report, ordered = TRUE)

} # }
```
