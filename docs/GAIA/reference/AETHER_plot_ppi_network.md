# Static PPI Network Plot

Creates a publication-quality network plot using ggraph. Nodes can be
colored by cluster membership, fold change, module, or any named vector.
Node size reflects degree or custom metric.

## Usage

``` r
AETHER_plot_ppi_network(
  graph,
  color_by = "degree",
  size_by = "degree",
  layout = "fr",
  show_labels = "hubs",
  hub_n = 15,
  colors = NULL,
  edge_alpha = 0.3,
  edge_color = "gray70",
  size_range = c(2, 10),
  group_labels = NULL,
  group_label_size = 3.5,
  group_label_color = "black",
  title = NULL,
  seed = 42
)
```

## Arguments

- graph:

  An igraph object (e.g., from APOLLO_get_ppi()).

- color_by:

  Node coloring. One of:

  - Named vector: gene -\> value (numeric for gradient, character/factor
    for discrete). E.g., cluster assignments or log2FC values.

  - "degree": Color by node degree (connectivity).

  - "betweenness": Color by betweenness centrality.

  - Single color string: uniform color (e.g., "steelblue").

  - NULL: default gray.

- size_by:

  Node sizing. One of:

  - "degree" (default): Size by connectivity.

  - "betweenness": Size by betweenness centrality.

  - Named numeric vector: gene -\> value.

  - Single number: uniform size.

- layout:

  Character. igraph layout algorithm. Default: "fr"
  (Fruchterman-Reingold). Other options: "kk" (Kamada-Kawai), "circle",
  "tree", "grid", "stress", "dh", "gem", "graphopt".

- show_labels:

  Logical or character. TRUE shows all labels, FALSE hides all, "hubs"
  shows labels only for top hub nodes. Default: "hubs".

- hub_n:

  Integer. Number of top hubs to label when show_labels="hubs". Default:
  15.

- colors:

  Color specification for nodes. For discrete: named vector of colors.
  For continuous: vector of 2-3 colors for gradient. Default: NULL
  (auto-selected based on color_by type).

- edge_alpha:

  Numeric. Edge transparency (0-1). Default: 0.3.

- edge_color:

  Character. Edge color. Default: "gray70".

- size_range:

  Numeric vector of length 2. Min and max node sizes. Default: c(2, 10).

- group_labels:

  Named character vector. Maps a discrete `color_by` group (e.g. PART
  cluster "C1", "C2", ...) to a label string (e.g. its top enriched GO
  term from
  [`APOLLO_enrich_gost()`](https://ylefol.github.io/L-Atelier/GAIA/reference/APOLLO_enrich_gost.md)).
  One label is placed at the centroid of each group's nodes. Only
  applies when `color_by` produces a discrete grouping (e.g. a named
  character vector of cluster assignments) — ignored with a warning for
  continuous coloring ("degree", "betweenness", or a numeric named
  vector). Default: NULL (no group labels).

- group_label_size:

  Numeric. Font size for group labels. Default: 3.5.

- group_label_color:

  Character. Text color for group labels. Default: "black".

- title:

  Character. Plot title. Default: NULL (auto-generated).

- seed:

  Integer. Random seed for reproducible layouts. Default: 42.

## Value

A ggplot object.

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic network colored by degree
p <- AETHER_plot_ppi_network(graph, color_by = "degree")

# Color by cluster membership
clusters <- c(TP53 = "C1", BRCA1 = "C1", EGFR = "C2", MYC = "C2")
p <- AETHER_plot_ppi_network(graph, color_by = clusters)

# Color by fold change
fc <- c(TP53 = 2.1, BRCA1 = -1.5, EGFR = 0.3, MYC = 3.2)
p <- AETHER_plot_ppi_network(graph, color_by = fc, colors = c("blue", "white", "red"))

# Label each PART cluster with its top enriched GO term
top_terms <- c(C1 = "immune response", C2 = "coagulation cascade")
p <- AETHER_plot_ppi_network(graph, color_by = clusters,
                              group_labels = top_terms)

} # }
```
