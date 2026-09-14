# Plot Sankey Diagram for Cluster Transitions Across k Values

Visualize how samples flow between clusters as k increases. This helps
understand cluster stability and identify natural substructure.

## Usage

``` r
AETHER_plot_cluster_sankey(
  k_eval,
  colors = NULL,
  node_width = 30,
  font_size = 12,
  title = "Cluster Transitions Across k",
  save_html = NULL
)
```

## Arguments

- k_eval:

  An artemis_k_evaluation object from ARTEMIS_evaluate_k_range().

- colors:

  Named character vector or NULL. Custom colors for clusters. If NULL
  (default), uses a generated palette.

- node_width:

  Numeric. Width of cluster nodes. Default = 30.

- font_size:

  Numeric. Font size for labels. Default = 12.

- title:

  Character. Plot title. Default = "Cluster Transitions Across k".

- save_html:

  Character. If provided, saves plot as standalone HTML file. Default =
  NULL (no save).

## Value

A plotly/networkD3 Sankey diagram object

## Details

The Sankey diagram shows:

- Columns: different k values (left to right)

- Nodes: clusters at each k value

- Flows: how samples move from clusters at k to clusters at k+1

- Width: number of samples in each flow

Interpretation:

- Clean splits (one cluster splits into two) suggest natural
  substructure

- Messy redistributions suggest artificial divisions

- Stable cores (samples staying together) indicate robust clusters

## Examples

``` r
if (FALSE) { # \dontrun{
k_eval <- ARTEMIS_evaluate_k_range(my_data, k_range = 2:5)
AETHER_plot_cluster_sankey(k_eval)

# Save as HTML for sharing
AETHER_plot_cluster_sankey(k_eval, save_html = "cluster_transitions.html")

} # }
```
