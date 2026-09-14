# GO DAG Plot coloured by module

Builds the GO Directed Acyclic Graph (DAG) induced by a set of enriched
GO terms, traverses upward toward the root to include ancestor context
nodes, and renders the result as a hierarchical network plot using
ggraph. Enriched nodes are coloured by module of origin; ancestor
(non-enriched) context nodes are shown in grey.

## Usage

``` r
AETHER_plot_go_dag(
  enrichment_result,
  ont = "BP",
  top_n = 20L,
  min_ancestor_freq = 2L,
  max_depth = 4L,
  show_shared = FALSE,
  edge_types = c("is_a", "part_of", "regulates", "positively_regulates",
    "negatively_regulates"),
  node_size_range = c(3, 10),
  palette = NULL,
  title = NULL,
  verbose = TRUE
)
```

## Arguments

- enrichment_result:

  A `gost_enrichment` object from
  [`APOLLO_enrich_gost`](https://ylefol.github.io/L-Atelier/GAIA/reference/APOLLO_enrich_gost.md)
  /
  [`DEMETER_load_enrichment`](https://ylefol.github.io/L-Atelier/GAIA/reference/DEMETER_load_enrichment.md),
  or a plain data.frame with columns `term_id`, `source`, `p_value`,
  `module`.

- ont:

  Character. GO ontology: `"BP"` (default), `"MF"`, or `"CC"`. Only one
  ontology per call.

- top_n:

  Integer. Top N enriched GO terms per module (by p-value) to include
  before building the DAG. Default: `20L`.

- min_ancestor_freq:

  Integer. Minimum number of enriched terms an ancestor node must be an
  ancestor of in order to be retained. Enriched terms themselves are
  always kept. Lowering this value shows more context; raising it
  focuses on shared hubs only. Default: `2L`.

- max_depth:

  Integer or `NULL`. Maximum number of levels to traverse upward from
  enriched terms. `NULL` traverses all the way to the ontology root.
  Default: `4L`.

- show_shared:

  Logical. If `FALSE` (default), each enriched term is coloured by its
  single best-p module. If `TRUE`, terms found in multiple modules get a
  combined label (e.g., `"C1-C2"`) and a distinct auto-assigned colour.

- edge_types:

  Character vector. GO relationship types to include as edges. Any
  combination of `"is_a"`, `"part_of"`, `"regulates"`,
  `"positively_regulates"`, `"negatively_regulates"`. Default: all five.

- node_size_range:

  Numeric vector of length 2. Min and max point sizes for enriched nodes
  (scaled by `-log10(p)`). Default: `c(3, 10)`.

- palette:

  Named character vector. Module name → hex color mapping. `NULL`
  (default) auto-assigns from internal palette.

- title:

  Character. Plot title. `NULL` (default) auto-generates.

- verbose:

  Logical. Print progress messages. Default: `TRUE`.

## Value

A `ggplot` / `ggraph` object. Save with
[`ggplot2::ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html).

## Details

**Algorithm:**

1.  Filters `enrichment_result$combined` to the requested ontology and
    takes the top `top_n` terms per module by p-value.

2.  Traverses the GO DAG upward (child → parent) from enriched terms
    using
    [`GO.db::GOBPPARENTS`](https://rdrr.io/pkg/GO.db/man/GOBPPARENTS.html)
    (or MF/CC), up to `max_depth` levels.

3.  For each ancestor node, counts how many enriched terms it is an
    ancestor of. Ancestors below `min_ancestor_freq` are removed.

4.  Renders the filtered DAG with ggraph using the Sugiyama
    (layered/hierarchical) layout. Enriched nodes are coloured by module
    group; ancestor nodes are grey. Edge line type encodes GO
    relationship.

**Required packages:** GO.db (Bioconductor), ggraph, igraph (both
already in Suggests).

## See also

[`APOLLO_enrich_gost`](https://ylefol.github.io/L-Atelier/GAIA/reference/APOLLO_enrich_gost.md),
[`AETHER_plot_go_treemap`](https://ylefol.github.io/L-Atelier/GAIA/reference/AETHER_plot_go_treemap.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic GO:BP DAG
p <- AETHER_plot_go_dag(gost_res, ont = "BP")
ggplot2::ggsave("go_dag_bp.png", p, width = 14, height = 10)

# Tighter focus: only high-frequency hubs, shallow traversal
p <- AETHER_plot_go_dag(gost_res, ont = "BP",
                         min_ancestor_freq = 3, max_depth = 3)

# is_a edges only (cleaner)
p <- AETHER_plot_go_dag(gost_res, ont = "BP", edge_types = "is_a")
} # }
```
