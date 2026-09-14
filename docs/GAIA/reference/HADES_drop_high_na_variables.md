# Remove Variables with High Missing Data

Simple helper to remove columns (variables) that exceed a missing data
threshold. This is a straightforward filter that doesn't require a QC
report.

## Usage

``` r
HADES_drop_high_na_variables(
  data,
  na_threshold = 0.2,
  protect = NULL,
  na_strings = c("NA", "N/A", "na", "n/a", "None", "none", "NONE", "NULL", "null", "NaN",
    "nan", "", " ", ".", "-", "--", "missing", "Missing", "MISSING", "#N/A", "#NA",
    "#VALUE!"),
  verbose = TRUE
)
```

## Arguments

- data:

  Data frame to filter.

- na_threshold:

  Numeric (0-1). Columns with NA proportion exceeding this threshold
  will be removed. Default = 0.2 (20%).

- protect:

  Character vector. Column names that should never be removed,
  regardless of their NA proportion (e.g., sample ID column).

- na_strings:

  Character vector. Strings to treat as NA values (e.g., "N/A", "NA",
  "None"). These will be counted as missing AND converted to actual NA
  in the returned data. Set to NULL to disable. Default includes common
  patterns.

- verbose:

  Logical. Print information about removed columns. Default = TRUE.

## Value

The data frame with high-NA columns removed (and NA-like strings
converted to NA if na_strings is not NULL).

## Details

This is a simple filtering function for quick data cleaning. For more
comprehensive QC assessment, use HADES_qc_mixed_data() first.

The function calculates the proportion of NA values in each column and
removes those exceeding the threshold. Protected columns are never
removed.

When na_strings is provided, matching strings in character/factor
columns are treated as missing data for threshold calculation and
converted to actual NA in the output.

## Examples

``` r
if (FALSE) { # \dontrun{
# Remove columns with >20% missing
data_clean <- HADES_drop_high_na_variables(my_data)

# Stricter threshold (>10% missing)
data_clean <- HADES_drop_high_na_variables(my_data, na_threshold = 0.1)

# Protect the ID column from removal
data_clean <- HADES_drop_high_na_variables(my_data, protect = "sample_id")

# Disable NA-string conversion
data_clean <- HADES_drop_high_na_variables(my_data, na_strings = NULL)

} # }
```
