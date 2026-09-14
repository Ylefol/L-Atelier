# Evaluate Clustering Across a Range of k Values

Systematically evaluates clustering quality across multiple k values and
tracks how samples move between clusters. This helps identify the
optimal number of clusters and visualize cluster stability.

## Usage

``` r
ARTEMIS_evaluate_k_range(
  data,
  k_range = 2:8,
  method = "pam",
  hclust_method = "ward.D2",
  stand = FALSE,
  verbose = TRUE
)
```

## Arguments

- data:

  Data frame with mixed variable types, OR a pre-computed distance
  matrix (class "dist").

- k_range:

  Integer vector. Range of k values to evaluate. Default = 2:8.

- method:

  Character. Clustering method: "pam" or "hierarchical". Default =
  "pam".

- hclust_method:

  Character. Linkage method for hierarchical clustering. Default =
  "ward.D2".

- stand:

  Logical. Standardize variables before computing Gower distance?
  Default = FALSE.

- verbose:

  Logical. Print progress messages. Default = TRUE.

## Value

A list with class "artemis_k_evaluation" containing:

- metrics:

  Data frame with k, silhouette_avg, and cluster size info

- assignments:

  Data frame with sample assignments at each k (wide format)

- transitions:

  Data frame tracking sample flows between consecutive k values

- optimal_k:

  Suggested optimal k based on silhouette

- distance:

  The distance matrix used

- k_range:

  The k values evaluated

## Details

This function is useful for:

- Determining optimal number of clusters (via silhouette scores)

- Understanding cluster stability (do samples stay together as k
  increases?)

- Visualizing sample flows with AETHER_plot_cluster_sankey()

The transitions output tracks how clusters at k split into clusters at
k+1, which can reveal natural substructure or artificial splits.

## Examples

``` r
if (FALSE) { # \dontrun{
# Evaluate k from 2 to 6
k_eval <- ARTEMIS_evaluate_k_range(my_data, k_range = 2:6)

# View metrics
print(k_eval)

# Visualize transitions
AETHER_plot_cluster_sankey(k_eval)

} # }
```
