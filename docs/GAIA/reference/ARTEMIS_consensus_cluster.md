# Consensus Clustering Across Multiple Imputed Datasets

Runs ConsensusClusterPlus on each imputed dataset and averages the
pairwise consensus matrices across all imputations. This extends the
standard single-dataset approach to account for imputation uncertainty,
following the method of Wick et al. (Oslo sepsis clustering, 2024) who
averaged matrices from 100 imputed datasets.

## Usage

``` r
ARTEMIS_consensus_cluster(
  imputed_obj,
  vars = NULL,
  log_vars = NULL,
  scale_data = TRUE,
  k_max = 10,
  n_resample = 100,
  resample_prop = 0.8,
  clusterAlg = "km",
  distance = "euclidean",
  seed = 42,
  verbose = TRUE
)
```

## Arguments

- imputed_obj:

  Either the output of
  [`POSEIDON_mice_impute()`](https://ylefol.github.io/L-Atelier/GAIA/reference/POSEIDON_mice_impute.md)
  or a plain named list of data frames (one per imputation).

- vars:

  Character vector. Variables to use for clustering. If NULL, all
  columns are used. Default = NULL.

- log_vars:

  Character vector. Variables to log1p-transform before scaling. Should
  match the variables used in
  [`POSEIDON_mice_impute()`](https://ylefol.github.io/L-Atelier/GAIA/reference/POSEIDON_mice_impute.md).
  Default = NULL.

- scale_data:

  Logical. Scale variables to mean=0, SD=1 before clustering. Default =
  TRUE.

- k_max:

  Integer. Maximum number of clusters to test. Default = 10.

- n_resample:

  Integer. Number of subsampling resamples per dataset in
  ConsensusClusterPlus. Default = 100.

- resample_prop:

  Numeric. Fraction of samples drawn per resample. Default = 0.8.

- clusterAlg:

  Character. Internal clustering algorithm: "km" (k-means) or "hc"
  (hierarchical). Default = "km".

- distance:

  Character. Distance metric. Default = "euclidean".

- seed:

  Integer. Base random seed; each imputation uses seed+i. Default = 42.

- verbose:

  Logical. Print progress. Default = TRUE.

## Value

A list with class "artemis_cc" containing:

- consensus_matrices:

  List (indexed by k) of averaged n×n consensus matrices

- cluster_assignments:

  List (indexed by k) of named integer cluster assignment vectors

- wcss:

  Named vector of within-cluster sum-of-squares per k

- cdf_area:

  Named vector of mean consensus value (area proxy) per k

- delta_area:

  Named vector of relative CDF area change between consecutive k

- cdf_values:

  List (indexed by k) of sorted upper-triangle values for CDF plotting

- k_max:

  Maximum k tested

- k_range:

  Integer vector 2:k_max

- n_imputations:

  Number of imputed datasets used

- n_samples:

  Number of samples

- sample_names:

  Sample identifiers

- vars_used:

  Variables used for clustering

- scaled_mean_data:

  Mean preprocessed data matrix (for centroid calculation)

## Details

Requires package `ConsensusClusterPlus`
(BiocManager::install("ConsensusClusterPlus")).

Cluster assignments are derived from the averaged consensus matrix via
hierarchical clustering (UPGMA) on (1 - consensus_matrix), then cutree
at k. WCSS is computed on the mean log1p-transformed, scaled dataset.
