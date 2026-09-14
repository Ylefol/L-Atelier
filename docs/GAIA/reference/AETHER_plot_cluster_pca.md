# PCA Plot Coloured by Consensus Cluster Assignment

Projects patients onto the first two principal components of the mean
preprocessed data (log1p-transformed, z-scored, averaged across all
imputations) and colours points by their consensus cluster assignment at
the chosen k.

## Usage

``` r
AETHER_plot_cluster_pca(
  cc_result,
  k,
  point_size = 1.5,
  alpha = 0.7,
  show_ellipse = TRUE,
  ellipse_level = 0.95,
  show_centroids = TRUE,
  centroid_size = 4,
  title = NULL
)
```

## Arguments

- cc_result:

  Output of
  [`ARTEMIS_consensus_cluster()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_consensus_cluster.md).

- k:

  Integer. Cluster solution to display.

- point_size:

  Numeric. Point size. Default = 1.5.

- alpha:

  Numeric. Point transparency (0–1). Default = 0.7.

- show_ellipse:

  Logical. Draw a 95% confidence ellipse per cluster. Default = TRUE.

- ellipse_level:

  Numeric. Confidence level for ellipses (0–1). Default = 0.95.

- show_centroids:

  Logical. Mark the per-cluster centroid with a filled symbol. Default =
  TRUE.

- centroid_size:

  Numeric. Size of the centroid symbol. Default = 4.

- title:

  Character. Plot title. If NULL, auto-generated. Default = NULL.

## Value

A ggplot object.

## Details

PCA is computed on `cc_result$scaled_mean_data`, which is already
centred and scaled, so `prcomp` is called with `center=FALSE` and
`scale.=FALSE`. Cluster separation in PC space is an indirect view of
the consensus structure: clusters defined by the consensus matrix (co-
clustering frequency) need not be linearly separable in PC space, so
partial overlap is expected for clinical data.
