# Cluster samples and detect outliers

Performs hierarchical clustering of samples to visualize relationships
and identify potential outliers. Can optionally remove outlier samples
based on a height cutoff.

## Usage

``` r
ARTEMIS_wgcna_cluster_samples(
  wgcna_data,
  method = "average",
  cut_height = NULL,
  min_cluster_size = 10,
  plot = TRUE,
  verbose = TRUE
)
```

## Arguments

- wgcna_data:

  A wgcna_data object or expression matrix (samples as rows).

- method:

  Clustering method for hclust(). Default: "average". Options:
  "average", "complete", "single", "ward.D", "ward.D2", etc.

- cut_height:

  Optional height cutoff for outlier removal. Samples that cluster above
  this height may be removed. Default: NULL (no removal).

- min_cluster_size:

  Minimum cluster size to keep when cutting. Default: 10.

- plot:

  Logical. Generate dendrogram plot. Default: TRUE.

- verbose:

  Logical. Print progress messages. Default: TRUE.

## Value

A list with class "wgcna_cluster" containing:

- sample_tree:

  hclust object

- datExpr:

  Expression matrix (possibly with outliers removed)

- datTraits:

  Traits data (if wgcna_data provided)

- outliers_removed:

  Names of removed outlier samples

- method:

  Clustering method used

- cut_height:

  Height cutoff used (if any)

## Examples

``` r
if (FALSE) { # \dontrun{
wgcna_data <- ARTEMIS_wgcna_prepare(counts, traits)
clust <- ARTEMIS_wgcna_cluster_samples(wgcna_data)

# With outlier removal
clust <- ARTEMIS_wgcna_cluster_samples(wgcna_data, cut_height = 100)

} # }
```
