# Fit Final Interpretable Model on Selected Variables

After variable selection with Group LASSO, fits an unpenalized
regression model on only the selected variables. This provides unbiased
coefficient estimates with proper p-values and confidence intervals for
inference and interpretation.

## Usage

``` r
ARTEMIS_fit_final_model(
  data,
  target_var,
  selected_vars,
  family = "binomial",
  verbose = TRUE
)
```

## Arguments

- data:

  Original data.frame (before encoding). Must contain the target
  variable and all selected predictor variables.

- target_var:

  Character string. Name of the target variable column.

- selected_vars:

  Character vector. Names of selected variables from
  `ARTEMIS_extract_selected_variables()$selected_names`.

- family:

  Character string. Response type:

  - "binomial": Binary outcome (logistic regression)

  - "multinomial": Multi-class outcome (3+ categories)

  - "gaussian": Continuous outcome (linear regression)

- verbose:

  Logical. Print model summary (default = TRUE).

## Value

A list of class "artemis_final_model" containing:

- model:

  The fitted glm object

- formula:

  The model formula used

- family:

  Response family

- selected_vars:

  Variables included in model

- n_samples:

  Number of samples used

- converged:

  Whether the model converged

## Details

Why refit without penalty?

- Group LASSO coefficients are shrunk toward zero (biased)

- P-values and CIs from penalized models are not valid for inference

- Refitting on selected variables gives unbiased estimates

- Standard errors, p-values, and CIs are now interpretable

This is the "post-selection inference" step. The selected variables were
chosen by Group LASSO; now we estimate their effects properly.

Note: P-values should be interpreted cautiously since variables were
pre-selected. The selection step already established these variables are
predictive; the final model quantifies their effects.

## Examples

``` r
if (FALSE) { # \dontrun{
# After selection
selected <- ARTEMIS_extract_selected_variables(fit, cv_result$lambda.1se)

# Fit final model on original data
final <- ARTEMIS_fit_final_model(
  data = clinical_data,
  target_var = "outcome",
  selected_vars = selected$selected_names,
  family = "binomial"
)

# View standard model summary
summary(final$model)

} # }
```
