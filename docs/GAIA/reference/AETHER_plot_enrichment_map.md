# Enrichment Map (Term-Term Network)

Visualises enriched pathway/GO terms as a network where nodes are terms
and edges connect terms that share a substantial fraction of their gene
sets (Jaccard similarity \\\ge\\ `jaccard_threshold`). Node size encodes
the best \\-\log\_{10}(p)\\ across all modules; node fill encodes the
primary module (lowest p-value). Edge width is proportional to Jaccard
similarity. A force-directed layout naturally clusters related terms
into "super-groups", making it easy to spot shared biology across
modules.

Terms with no edges above the threshold still appear as isolated nodes —
they represent module-specific biology not shared with other enriched
terms.

## Usage

``` r
AETHER_plot_enrichment_map(
  enrichment_result,
  source = "REAC",
  top_n = 20L,
  jaccard_threshold = 0.2,
  layout = "fr",
  node_size_range = c(3, 12),
  edge_width_range = c(0.3, 2),
  max_label_width = 40L,
  label_size = 2.5,
  palette = NULL,
  title = NULL,
  verbose = TRUE
)
```

## Arguments

- enrichment_result:

  A `gost_enrichment` object from
  [`APOLLO_enrich_gost()`](https://ylefol.github.io/L-Atelier/GAIA/reference/APOLLO_enrich_gost.md)
  or a plain data.frame with columns `term_id`, `term_name`, `source`,
  `module`, `p_value`, and optionally `intersection`.

- source:

  Character(1). Database to display. Default `"REAC"`.

- top_n:

  Integer. Top N terms per module (by p-value) pooled before
  deduplication. Default `20L`. Reduce to 10–15 for cleaner graphs.

- jaccard_threshold:

  Numeric in `[0, 1]`. Minimum Jaccard similarity to draw an edge.
  Default `0.2`. Increase to reduce edge density; decrease if too few
  edges appear.

- layout:

  Character. igraph/ggraph layout algorithm. `"fr"`
  (Fruchterman-Reingold, default) and `"kk"` (Kamada-Kawai) both work
  well; `"fr"` tends to produce rounder, more separated clusters.

- node_size_range:

  Numeric(2). Min and max point size for nodes. Default `c(3, 12)`.

- edge_width_range:

  Numeric(2). Min and max edge linewidth. Default `c(0.3, 2)`.

- max_label_width:

  Integer. Maximum label characters before truncation. Default `40L`.

- label_size:

  Numeric. `ggrepel` text size. Default `2.5`.

- palette:

  Named character vector of colours, one per module. `NULL` (default)
  auto-assigns from the internal 20-colour palette.

- title:

  Character. Plot title. `NULL` auto-generates.

- verbose:

  Logical. Default `TRUE`.

## Value

A `ggplot` object (ggsave-compatible).

## Details

Jaccard similarity is computed from the `intersection` gene lists in the
`gost_enrichment` object (requires `evcodes = TRUE` in
`APOLLO_enrich_gost`). The union of intersection gene sets across all
modules is used per term, so similarity reflects shared biology rather
than query-specific overlap.

Requires ggraph and igraph (both already in Suggests).

## Examples

``` r
if (FALSE) { # \dontrun{
p <- AETHER_plot_enrichment_map(enrich, source = "REAC", top_n = 20)
ggplot2::ggsave("emap.png", p, width = 12, height = 10)

# Fewer terms, tighter edges
p <- AETHER_plot_enrichment_map(enrich, top_n = 10, jaccard_threshold = 0.3)
} # }
```
