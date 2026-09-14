# Cluster Stability Analysis for Consensus Clustering

Assesses the stability of a chosen k by running ConsensusClusterPlus
`n_runs` times on the same imputed dataset with different random seeds.
For each patient, counts how often they are assigned to the same (modal)
cluster across all runs.

## Usage

``` r
ARTEMIS_consensus_stability(
  imputed_obj,
  k,
  n_runs = 10,
  vars = NULL,
  log_vars = NULL,
  scale_data = TRUE,
  n_resample = 100,
  resample_prop = 0.8,
  clusterAlg = "km",
  distance = "euclidean",
  imputation_idx = 1,
  seed = 42,
  verbose = TRUE
)
```

## Arguments

- imputed_obj:

  Output of
  [`POSEIDON_mice_impute()`](https://ylefol.github.io/L-Atelier/GAIA/reference/POSEIDON_mice_impute.md)
  or a list of data frames.

- k:

  Integer. Number of clusters to assess.

- n_runs:

  Integer. Number of independent CC runs. Default = 10.

- vars:

  Character vector. Variables to use. If NULL, all columns used. Default
  = NULL.

- log_vars:

  Character vector. Variables to log1p-transform before scaling. Default
  = NULL.

- scale_data:

  Logical. Scale data before clustering. Default = TRUE.

- n_resample:

  Integer. Resamples per CC run. Default = 100.

- resample_prop:

  Numeric. Subsample fraction. Default = 0.8.

- clusterAlg:

  Character. "km" or "hc". Default = "km".

- distance:

  Character. Distance metric. Default = "euclidean".

- imputation_idx:

  Integer. Which imputed dataset to use. Default = 1.

- seed:

  Integer. Base seed; run i uses seed+i. Default = 42.

- verbose:

  Logical. Print progress. Default = TRUE.

## Value

A list containing:

- assignment_matrix:

  Integer matrix (n_samples × n_runs) of cluster assignments per run

- stability_counts:

  Integer vector — for each patient, number of runs in their modal
  cluster

- pct_stable:

  Numeric — percentage of patients in the same cluster across ALL runs

- k:

  Number of clusters assessed

- n_runs:

  Number of runs performed

- sample_names:

  Sample identifiers

## Details

Requires package `ConsensusClusterPlus`.
