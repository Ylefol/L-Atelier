# Interactive PPI Network Plot

Creates an interactive network visualization using visNetwork. Nodes can
be dragged, zoomed, and hovered for information. Useful for exploring
large networks.

## Usage

``` r
AETHER_plot_ppi_network_interactive(
  graph,
  color_by = "degree",
  size_by = "degree",
  colors = NULL,
  size_range = c(10, 40),
  show_labels = TRUE,
  physics = TRUE,
  layout = "layout_with_fr",
  group_labels = NULL,
  title = NULL,
  save_html = NULL,
  seed = 42
)
```

## Arguments

- graph:

  An igraph object (e.g., from APOLLO_get_ppi()).

- color_by:

  Node coloring. Same options as AETHER_plot_ppi_network(): named
  vector, "degree", "betweenness", single color, or NULL.

- size_by:

  Node sizing. Same options as AETHER_plot_ppi_network().

- colors:

  Color specification. Named vector for discrete, or gradient vector for
  continuous. Default: NULL (auto-selected).

- size_range:

  Numeric vector of length 2. Min and max node sizes. Default: c(10,
  40).

- show_labels:

  Logical. Show node labels. Default: TRUE.

- physics:

  Logical. Enable physics simulation for layout. Default: TRUE. Set
  FALSE for static layout (faster for large networks).

- layout:

  Character. Layout algorithm when physics=FALSE. Options:
  "layout_with_fr", "layout_nicely", "layout_in_circle". Default:
  "layout_with_fr".

- group_labels:

  Named character vector. Maps a discrete `color_by` group (e.g. PART
  cluster "C1", "C2", ...) to a label string (e.g. its top enriched GO
  term from
  [`APOLLO_enrich_gost()`](https://ylefol.github.io/L-Atelier/GAIA/reference/APOLLO_enrich_gost.md)).
  Added to each matching node's hover tooltip as "Cluster" and "Top
  term" lines (no separate legend is drawn — see
  [`AETHER_plot_ppi_network`](https://ylefol.github.io/L-Atelier/GAIA/reference/AETHER_plot_ppi_network.md)'s
  static `group_labels` for that). Only applies when `color_by` produces
  a discrete grouping — ignored with a warning otherwise. Default: NULL.

- title:

  Character. Plot title. Default: NULL.

- save_html:

  Character. Path to save as standalone HTML file. Default: NULL.

- seed:

  Integer. Random seed for layout. Default: 42.

## Value

A visNetwork object (rendered in viewer or saved to HTML).

## Examples

``` r
if (FALSE) { # \dontrun{
# Interactive exploration
AETHER_plot_ppi_network_interactive(graph, color_by = "degree")

# Save to HTML
AETHER_plot_ppi_network_interactive(graph, save_html = "network.html")

# Show cluster + top enriched term in each node's tooltip
top_terms <- c(C1 = "immune response", C2 = "coagulation cascade")
AETHER_plot_ppi_network_interactive(graph, color_by = clusters,
                                     group_labels = top_terms)

} # }
```
