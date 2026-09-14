# Create K-Fold Cross-Validation Indices

Creates k-fold CV splits, ensuring same samples are held out across all
datasets (maintains sample pairing)

## Usage

``` r
MINERVA_create_folds(n_samples, k = 5, seed = NULL)
```

## Arguments

- n_samples:

  Number of samples

- k:

  Number of folds (default = 5)

- seed:

  Random seed for reproducibility

## Value

List of k vectors, each containing test indices for that fold
