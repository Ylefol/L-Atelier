# Compute Sparsity Quality Score

Scores sparsity based on how well it falls within ideal range. Penalizes
both extreme sparsity (too few features) and lack of sparsity (too many
features).

## Usage

``` r
MINERVA_sparsity_quality(
  sparsity,
  ideal_min = 0.05,
  ideal_max = 0.3,
  penalty_steepness = 2
)
```

## Arguments

- sparsity:

  Numeric. Proportion of non-zero features (0 to 1)

- ideal_min:

  Numeric. Minimum ideal sparsity (default = 0.05, or 5%)

- ideal_max:

  Numeric. Maximum ideal sparsity (default = 0.30, or 30%)

- penalty_steepness:

  Numeric. How quickly to penalize outside ideal range (default = 2)

## Value

Numeric quality score between 0 and 1, where 1 = ideal sparsity

## Details

The quality score is:

- 1.0 if sparsity is within `ideal_min` to `ideal_max`

- Decreases as sparsity moves outside this range

- Uses exponential decay based on penalty_steepness
