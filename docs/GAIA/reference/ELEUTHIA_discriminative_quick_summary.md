# Quick Summary of Discriminative Analysis

Print a quick one-page summary without file export.

## Usage

``` r
ELEUTHIA_discriminative_quick_summary(glasso_fit, cv_result, selected)
```

## Arguments

- glasso_fit:

  An artemis_group_lasso object.

- cv_result:

  Output from ARTEMIS_select_lambda().

- selected:

  Output from ARTEMIS_extract_selected_variables().

## Value

Invisibly returns NULL.
