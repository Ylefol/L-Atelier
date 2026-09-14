# Match cluster labels between two clusterings

Internal helper to handle label switching problem. Finds the best
mapping from cluster labels in clust2 to clust1.

## Usage

``` r
.match_cluster_labels(clust1, clust2)
```

## Arguments

- clust1:

  Original cluster assignments

- clust2:

  New cluster assignments to map

## Value

Named integer vector: mapping from clust2 labels to clust1 labels
