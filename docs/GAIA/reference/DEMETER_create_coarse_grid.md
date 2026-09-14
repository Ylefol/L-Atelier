# Create Coarse Parameter Grid

Generates a coarse grid of parameter combinations for initial
exploration in cross-validation or grid search. Creates all combinations
of the specified values across n parameters.

## Usage

``` r
DEMETER_create_coarse_grid(n_params, values = c(0.3, 0.5, 0.7))
```

## Arguments

- n_params:

  Integer. Number of parameters to create grid for.

- values:

  Numeric vector. Values to use for each parameter (default = c(0.3,
  0.5, 0.7)).

## Value

Data frame with n_params columns (param1, param2, ..., paramN)
containing all combinations of the specified values. Each row represents
one parameter combination. Total rows = length(values)^n_params.

## Details

This function is useful for initial broad exploration of parameter space
before refining the search with a fine grid around promising regions.

For example, with 3 parameters and 3 values each, this creates 27 (3^3)
combinations. This coarse grid helps identify the general region of
optimal parameters without excessive computation.

## Examples

``` r
if (FALSE) { # \dontrun{
# Create coarse grid for 3 tau parameters
coarse_grid <- DEMETER_create_coarse_grid(n_params = 3)
# Returns 27 combinations of (0.3, 0.5, 0.7)

# Custom coarse values
coarse_grid <- DEMETER_create_coarse_grid(
  n_params = 2,
  values = c(0.1, 0.5, 0.9)
)

} # }
```
