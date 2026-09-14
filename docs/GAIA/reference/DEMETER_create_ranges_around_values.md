# Create Ranges Around Best Parameter Values

Helper function to automatically generate ranges for fine grid search
centered around best parameter values from coarse grid search.

## Usage

``` r
DEMETER_create_ranges_around_values(best_values, width = 0.2, bounds = c(0, 1))
```

## Arguments

- best_values:

  Numeric vector. Best parameter values from coarse search.

- width:

  Numeric. Total width of range around each value (default = 0.2). Range
  will be `value - width/2` to `value + width/2`.

- bounds:

  Numeric vector of length 2. Lower and upper bounds to clip ranges
  (default = c(0, 1)). Useful for parameters with valid ranges (e.g.,
  tau in 0 to 1).

## Value

List of numeric vectors. Each element is c(min, max) for one parameter,
suitable for input to DEMETER_create_fine_grid().

## Details

This convenience function automates the creation of fine grid ranges by:

- Centering each range around the corresponding best value

- Ensuring symmetric exploration (±width/2)

- Clipping to valid parameter bounds

This is particularly useful in two-stage grid search:

1.  Run coarse grid CV to find best parameters

2.  Use this function to create ranges around best values

3.  Run fine grid CV for refined tuning

## Examples

``` r
if (FALSE) { # \dontrun{
# After coarse CV, best tau values were (0.5, 0.7, 0.5)
ranges <- DEMETER_create_ranges_around_values(
  best_values = c(0.5, 0.7, 0.5),
  width = 0.2,
  bounds = c(0, 1)
)
# Returns: list(c(0.4, 0.6), c(0.6, 0.8), c(0.4, 0.6))

# Then create fine grid
fine_grid <- DEMETER_create_fine_grid(ranges, n_points = 5)

# Larger search width
ranges <- DEMETER_create_ranges_around_values(
  best_values = c(0.3, 0.8),
  width = 0.4,
  bounds = c(0, 1)
)
# Returns: list(c(0.1, 0.5), c(0.6, 1.0))

} # }
```
