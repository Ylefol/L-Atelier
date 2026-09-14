# Update Canonical Vector for Multi-Dataset ConvCCA

Updates one canonical vector in a multi-dataset setting by maximizing
correlation with all other datasets

## Usage

``` r
HEPHAESTUS_updateWmulti_convCCA(K, w, whichDSet, tau, penalty = "LASSO")
```

## Arguments

- K:

  List of correlation matrices between datasets

- w:

  List of current canonical vectors for all datasets

- whichDSet:

  Index of the dataset being updated

- tau:

  Sparsity tuning parameter for this dataset

- penalty:

  Penalty function: "LASSO" (default) or "SCAD"

## Value

Updated canonical vector for dataset whichDSet
