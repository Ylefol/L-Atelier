# Cluster Mixed Data Using Gower Distance

Clustering for datasets containing both quantitative and qualitative
variables. Uses Gower distance which handles mixed types natively, then
applies PAM (partitioning around medoids) or hierarchical clustering.

## Usage

``` r
ARTEMIS_cluster_mixed(
  data,
  k = "auto",
  k_range = 2:10,
  method = "pam",
  hclust_method = "ward.D2",
  stand = FALSE,
  verbose = TRUE
)
```

## Arguments

- data:

  Data frame with mixed variable types. Should be pre-cleaned (no
  high-NA columns, appropriate factor conversions done).

- k:

  Integer or "auto". Number of clusters. If "auto", determines optimal k
  using silhouette scores over k_range. Default = "auto".

- k_range:

  Integer vector. Range of k values to test when k = "auto". Default =
  2:10.

- method:

  Character. Clustering method: "pam" (partitioning around medoids) or
  "hierarchical". Default = "pam".

- hclust_method:

  Character. Linkage method for hierarchical clustering. One of
  "ward.D2", "complete", "average", "single". Default = "ward.D2".

- stand:

  Logical. Standardize variables before computing Gower distance?
  Default = FALSE (Gower handles scale differences internally).

- verbose:

  Logical. Print progress messages. Default = TRUE.

## Value

A list with class "artemis_cluster" containing:

- clusters:

  Integer vector of cluster assignments

- k:

  Number of clusters used

- k_selection:

  If k="auto", data frame with silhouette scores per k

- k_evaluation:

  If k="auto", full evaluation data (class "artemis_k_evaluation")
  containing assignments at each k and transition data for Sankey
  visualization. Pass to AETHER_plot_cluster_sankey() or
  AETHER_plot_k_selection(). NULL if k was specified directly.

- silhouette:

  Silhouette information for final clustering

- silhouette_avg:

  Average silhouette width (cluster quality measure)

- medoids:

  For PAM: indices of medoid samples

- distance:

  The Gower distance matrix

- method:

  Clustering method used

- data_used:

  The data frame used for clustering

- call:

  The function call

## Details

Gower distance handles mixed data by computing appropriate distances for
each variable type:

- Quantitative: Manhattan distance, normalized by range

- Qualitative: Simple matching (0 if same, 1 if different)

- Ordinal: Ranked, then treated as quantitative

PAM is generally preferred over k-means for mixed data because:

- Works with any distance metric (not just Euclidean)

- More robust to outliers (uses medoids, not centroids)

- Medoids are actual data points, aiding interpretation

## Examples

``` r
if (FALSE) { # \dontrun{
# Auto-select k
result <- ARTEMIS_cluster_mixed(my_data)

# Visualize k selection and cluster transitions (when k="auto")
AETHER_plot_k_selection(result$k_evaluation)
AETHER_plot_cluster_sankey(result$k_evaluation)

# Specify k directly (no k_evaluation generated)
result <- ARTEMIS_cluster_mixed(my_data, k = 3)

# Use hierarchical clustering
result <- ARTEMIS_cluster_mixed(my_data, method = "hierarchical")

} # }
```
