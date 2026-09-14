# Quality Control Assessment for Mixed Data

Assess a mixed (quantitative + qualitative) dataset and report QC
metrics to inform filtering decisions. This function reports and flags
potential issues but does NOT automatically remove anything - the user
reviews the report and decides what to filter.

## Usage

``` r
HADES_qc_mixed_data(
  data,
  na_threshold = 0.2,
  variance_threshold = 1e-10,
  dominance_threshold = 0.99,
  correlation_threshold = 0.9,
  check_correlations = TRUE,
  check_samples = TRUE,
  check_type_issues = TRUE,
  na_strings = c("NA", "N/A", "na", "n/a", "Na", "N/a", "None", "none", "NONE", "NULL",
    "null", "Null", "NaN", "nan", "NAN", "", " ", "  ", ".", "-", "--", "---", "missing",
    "Missing", "MISSING", "not available", "Not Available", "NOT AVAILABLE", "unknown",
    "Unknown", "UNKNOWN", "unk", "UNK", "#N/A", "#NA", "#VALUE!", "#REF!", "#DIV/0!",
    "#NAME?", "n.a.", "N.A.", "n.a", "N.A", "nil", "NIL", "Nil", "undefined",
    "UNDEFINED"),
  numeric_threshold = 0.5,
  categorical_threshold = 10,
  sample_id_col = NULL,
  verbose = TRUE
)
```

## Arguments

- data:

  Data frame with mixed variable types.

- na_threshold:

  Numeric (0-1). Variables/samples with NA proportion exceeding this
  threshold will be flagged. Default = 0.2 (20%).

- variance_threshold:

  Numeric. Quantitative variables with variance below this threshold
  will be flagged as near-zero variance. Default = 1e-10.

- dominance_threshold:

  Numeric (0-1). Qualitative variables where one level represents more
  than this proportion will be flagged as near-constant. Default = 0.99
  (99%).

- correlation_threshold:

  Numeric (0-1). Pairs of quantitative variables with absolute
  correlation above this threshold will be reported. Default = 0.9.

- check_correlations:

  Logical. Whether to compute correlation matrix for quantitative
  variables. Default = TRUE.

- check_samples:

  Logical. Whether to compute sample-level QC metrics. Default = TRUE.

- check_type_issues:

  Logical. Whether to check for type inconsistencies and NA-like
  strings. Default = TRUE.

- na_strings:

  Character vector. Strings to detect as potential NA values. Default
  includes common patterns: "N/A", "NA", "None", "NULL", "", etc. Set to
  NULL to disable NA-string detection.

- numeric_threshold:

  Numeric (0-1). For character columns, if this proportion of non-NA
  values can be parsed as numeric, flag as potential type issue. Default
  = 0.5 (50%).

- categorical_threshold:

  Integer. Numeric columns with this many or fewer unique values will be
  flagged as potentially categorical (might need to be converted to
  factor for FAMD). Default = 10.

- sample_id_col:

  Character or integer. Column name or index containing sample IDs.
  Recommended to specify this for meaningful sample identification. When
  specified, the column is validated for uniqueness, missing values, and
  empty strings. If NULL (default), row names or row numbers are used.

- verbose:

  Logical. Print progress messages. Default = TRUE.

## Value

A list with class "hades_qc" containing:

- variables:

  Data frame with per-variable QC metrics and flags

- samples:

  Data frame with per-sample QC metrics and flags (if check_samples =
  TRUE)

- correlations:

  Data frame of highly correlated variable pairs (if check_correlations
  = TRUE)

- type_issues:

  List with type inconsistencies, NA-like strings, and special
  characters found (if check_type_issues = TRUE). Includes:
  potential_numeric, potential_categorical, na_like_strings,
  mixed_na_representations, special_chars_in_values

- id_issues:

  List with sample ID validation results (if sample_id_col specified)

- summary:

  List with counts of flagged variables/samples

- thresholds:

  List of thresholds used for flagging

## Details

This function is designed for use before FAMD or other mixed-data
analyses. It helps identify:

- Variables with excessive missing data

- Variables with near-zero variance (uninformative)

- Highly correlated variable pairs (redundancy)

- Samples with excessive missing data

- Type inconsistencies (character columns that look numeric)

- NA-like strings that weren't parsed as NA (e.g., "N/A", "None")

- Numeric columns that may be categorical (few unique values)

- Special characters in categorical values (spaces, underscores) that
  may cause issues with FAMD category naming

The user should review the output and make informed decisions about what
to filter before proceeding with analysis.

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic QC assessment
qc_report <- HADES_qc_mixed_data(my_data)
print(qc_report)

# View flagged variables
qc_report$variables[qc_report$variables$flagged_any, ]

# View highly correlated pairs
qc_report$correlations

# View type issues (NA-like strings, potential numeric columns)
qc_report$type_issues

} # }
```
