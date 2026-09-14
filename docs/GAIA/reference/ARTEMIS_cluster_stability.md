# Assess Cluster Stability via Bootstrap

Evaluates how stable cluster assignments are under bootstrap resampling.
High stability suggests the clustering captures real structure rather
than noise.

## Usage

``` r
ARTEMIS_cluster_stability(cluster_result, n_boot = 100, verbose = TRUE)
```

## Arguments

- cluster_result:

  An artemis_cluster object from ARTEMIS_cluster_mixed().

- n_boot:

  Integer. Number of bootstrap iterations. Default = 100.

- verbose:

  Logical. Print progress. Default = TRUE.

## Value

A list with class "artemis_stability" containing:

- jaccard_mean:

  Mean Jaccard similarity across bootstrap samples

- jaccard_per_cluster:

  Jaccard similarity for each cluster

- sample_stability:

  Per-sample stability score (proportion of times assigned to same
  cluster as in original)

- n_boot:

  Number of bootstrap iterations performed

- interpretation:

  Text interpretation of stability

## Details

For each bootstrap iteration:

1.  Resample data with replacement

2.  Re-cluster using same parameters

3.  Compare cluster assignments for samples present in both

Jaccard similarity measures overlap between original and bootstrap
clusters. Values \> 0.75 indicate stable clusters; \< 0.5 indicates
instability.
