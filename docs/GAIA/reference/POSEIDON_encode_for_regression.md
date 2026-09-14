# Encode Mixed Data for Regression Analysis

Prepares a mixed (categorical + quantitative) dataset for regression
analysis by separating the target variable, encoding categorical
predictors, and creating group mappings for Group LASSO.

## Usage

``` r
POSEIDON_encode_for_regression(
  data,
  target_var,
  encoding = "dummy",
  reference = "first",
  ordinal_vars = NULL,
  ordinal_treatment = "categorical",
  drop_target_na = TRUE,
  verbose = TRUE
)
```

## Arguments

- data:

  A data.frame with mixed data types.

- target_var:

  Character string. Name of the target variable column.

- encoding:

  Character string. Encoding method for categorical variables:

  - "dummy" (default): Treatment/dummy coding. Creates k-1 binary
    columns for k levels, with one level as reference.

  - "effect": Effect coding. Similar to dummy but reference level coded
    as -1 instead of 0 (compares to grand mean).

- reference:

  Character string or named list. How to choose reference level:

  - "first" (default): First level alphabetically

  - "last": Last level alphabetically

  - Named list: Specific reference per variable, e.g., list(sex =
    "male", treatment = "placebo")

- ordinal_vars:

  Character vector. Names of variables that are ordinal (ordered
  categorical). These can be treated differently from nominal
  categoricals.

- ordinal_treatment:

  Character string. How to treat ordinal variables:

  - "categorical" (default): Treat as nominal (dummy coding)

  - "continuous": Treat as numeric (1, 2, 3, ...)

  - "polynomial": Use polynomial contrasts (linear, quadratic, etc.)

- drop_target_na:

  Logical. Remove rows where target is NA (default = TRUE).

- verbose:

  Logical. Print encoding summary (default = TRUE).

## Value

A list containing:

- X:

  Numeric matrix of encoded predictors (samples x features)

- y:

  Vector of target variable values

- var_mapping:

  List with encoding information:

  - original_vars: Original variable names

  - var_types: Type of each original variable
    (quantitative/categorical/ordinal)

  - groups: Integer vector mapping each column in X to its variable
    group (for Group LASSO - all columns from same variable share group
    number)

  - group_names: Variable name for each group number

  - col_to_var: Maps each X column name to original variable

  - reference_levels: Reference level for each categorical variable

  - encoding: Encoding method used

- sample_ids:

  Row identifiers from original data

- target_name:

  Name of target variable

- n_removed_na:

  Number of rows removed due to NA in target

## Details

This function is the first step in a supervised variable selection
workflow. The output is designed to work with Group LASSO (via grpreg
package), where the `groups` vector ensures that all dummy columns from
a single categorical variable are treated as a group (either all in or
all out).

For ordinal variables (like severity scores 1-9), the choice of
treatment affects interpretation:

- "categorical": Most flexible, allows non-linear effects, each level
  can have independent effect (recommended for selection stage)

- "continuous": Assumes linear effect across ordered levels

- "polynomial": Captures non-linear ordered effects with fewer
  parameters

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic encoding
encoded <- POSEIDON_encode_for_regression(
  data = clinical_data,
  target_var = "outcome"
)

# With specific reference levels
encoded <- POSEIDON_encode_for_regression(
  data = clinical_data,
  target_var = "response",
  reference = list(treatment = "placebo", sex = "female")
)

# With ordinal variables treated as continuous
encoded <- POSEIDON_encode_for_regression(
  data = clinical_data,
  target_var = "survival",
  ordinal_vars = c("severity_score", "stage"),
  ordinal_treatment = "continuous"
)

} # }
```
