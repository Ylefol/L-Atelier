# Factor Analysis of Mixed Data (FAMD)

Wrapper around FactoMineR::FAMD for dimensionality reduction of datasets
containing both quantitative (continuous) and qualitative (categorical)
variables.

## Usage

``` r
ARTEMIS_famd(
  data,
  na_action = "fail",
  na_threshold = 0.2,
  ncp = 5,
  sup_quanti = NULL,
  sup_quali = NULL,
  ncp_impute = 5,
  prefix_levels = TRUE,
  verbose = TRUE
)
```

## Arguments

- data:

  Data frame with mixed variable types. Must contain at least one
  quantitative and one qualitative variable.

- na_action:

  Character. How to handle missing values:

  - "fail" (default): Error if any NA values present after column
    filtering

  - "omit": Remove rows with any NA values

  - "impute": Use missMDA::imputeFAMD() for iterative imputation

- na_threshold:

  Numeric (0-1). Columns with proportion of NA values exceeding this
  threshold are dropped before other NA handling. Default = 0.2 (20%).

- ncp:

  Integer. Number of dimensions to keep in results. Default = 5.

- sup_quanti:

  Integer vector. Indices of quantitative supplementary variables (not
  used in computation, but projected onto results).

- sup_quali:

  Integer vector. Indices of qualitative supplementary variables (not
  used in computation, but projected onto results).

- ncp_impute:

  Integer. Number of components used for imputation when na_action =
  "impute". Default = 5.

- prefix_levels:

  Logical. If TRUE (default), prefix factor levels with variable names
  (e.g., "0" becomes "gender_0"). This prevents ambiguity when multiple
  factors have overlapping numeric levels and makes FAMD output more
  readable.

- verbose:

  Logical. Print progress messages. Default = TRUE.

## Value

A list with class "artemis_famd" containing:

- famd:

  The FactoMineR FAMD result object

- data_used:

  The data frame used for analysis (after NA handling)

- eigenvalues:

  Data frame with eigenvalue, variance percent, and cumulative percent

- ind:

  List with individual (sample) coordinates, cos2, and contributions

- var:

  List with variable information (quanti and quali separately)

- na_report:

  Summary of NA handling actions taken

- level_mapping:

  List mapping original factor levels to prefixed versions (NULL if
  prefix_levels = FALSE)

- call:

  The function call

## Details

FAMD is particularly useful for exploratory analysis of sample metadata
that contains both continuous variables (age, signal intensity) and
categorical variables (treatment, batch, sex).

The function performs the following steps:

1.  Validate input (check for mixed types)

2.  Drop columns exceeding NA threshold

3.  Handle remaining NAs based on na_action

4.  Run FAMD analysis

5.  Structure and return results

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic usage with sample metadata
result <- ARTEMIS_famd(sample_metadata)

# With imputation for missing values
result <- ARTEMIS_famd(sample_metadata, na_action = "impute")

# With supplementary variables
result <- ARTEMIS_famd(sample_metadata, sup_quali = c(1))  # First column as supplementary

} # }
```
