# Get protein-protein interaction network directly from STRING

Queries the STRING REST API directly (rather than via OmniPath) for a
set of genes, returning both an igraph network object (with the same
shape as
[`APOLLO_get_ppi`](https://ylefol.github.io/L-Atelier/GAIA/reference/APOLLO_get_ppi.md),
so it can be used as a drop-in replacement) and, optionally, STRING's
own rendered network image. STRING's image renderer uses its own layout
engine, which is often better at avoiding node/label overlap than
re-laying out the same edges with
[`AETHER_plot_ppi_network`](https://ylefol.github.io/L-Atelier/GAIA/reference/AETHER_plot_ppi_network.md).

## Usage

``` r
APOLLO_get_string_ppi(
  genes,
  organism = 9606,
  score_threshold = 400,
  network_type = c("functional", "physical"),
  add_nodes = 0,
  image_format = c("image", "highres_image", "svg", "none"),
  image_path = NULL,
  white_background = TRUE,
  verbose = TRUE
)
```

## Arguments

- genes:

  Character vector of gene symbols to include in the network.

- organism:

  Integer. NCBI taxonomy ID. Default: 9606 (human).

- score_threshold:

  Integer (0-1000). STRING's `required_score` (combined confidence score
  cutoff). Default: 400 (STRING's own "medium confidence" default).

- network_type:

  Character. "functional" (default; includes predicted/ indirect
  functional associations, STRING's website default) or "physical"
  (direct physical binding interactions only).

- add_nodes:

  Integer. Number of extra first-shell interactors STRING should add
  beyond the query gene list. Default: 0 (within-list network only,
  matching
  [`APOLLO_get_ppi()`](https://ylefol.github.io/L-Atelier/GAIA/reference/APOLLO_get_ppi.md)'s
  default behavior).

- image_format:

  Character. Which STRING-rendered image to download: "image" (default,
  PNG), "highres_image" (high-resolution PNG), "svg" (vector), or "none"
  (skip image download entirely).

- image_path:

  Character. File path to save the downloaded image. Required when
  `image_format != "none"`.

- verbose:

  Logical. Print progress and network summary. Default: TRUE.

## Value

A list of class "apollo_string_ppi" with:

- graph:

  igraph object — same vertex/edge/graph attribute shape as
  [`APOLLO_get_ppi()`](https://ylefol.github.io/L-Atelier/GAIA/reference/APOLLO_get_ppi.md)'s
  return value.

- string_ids:

  data.frame. Raw STRING ID mapping for the query genes.

- interaction_df:

  data.frame. Edge list with STRING combined scores.

- image_path:

  Character or NULL. Path to the downloaded image, if any.

- params:

  List of the call parameters used.

## Details

This calls the STRING REST API (string-db.org/api) directly over HTTP
using POST requests (via the `curl` package), with identifiers sent in
the request body rather than a GET query string — this avoids the
URL-length ceiling that GET hits once gene lists run into the hundreds
(STRING's own documented recommendation for large identifier lists).
Gene symbols are first resolved to STRING IDs via the `get_string_ids`
endpoint; genes that fail to resolve are reported and kept as isolated
nodes in the returned graph (for parity with
[`APOLLO_get_ppi()`](https://ylefol.github.io/L-Atelier/GAIA/reference/APOLLO_get_ppi.md)'s
`n_query_genes`/`n_mapped_genes` graph attributes).

On a non-200 STRING API response, the actual HTTP status code and
response body are included in the error (rather than the generic
connection-failure message a GET-based
[`url()`](https://rdrr.io/r/base/connections.html) connection would
give), to make transient vs. persistent failures easier to tell apart.

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic STRING network, with STRING's own rendered image
res <- APOLLO_get_string_ppi(c("TP53", "BRCA1", "EGFR", "MYC", "CDK2"),
                              image_path = "string_network.png")

# Data only, no image, stricter confidence threshold
res <- APOLLO_get_string_ppi(genes, score_threshold = 700,
                              image_format = "none")

# Feed the graph into the existing ggraph-based plot for comparison
p <- AETHER_plot_ppi_network(res$graph, layout = "stress")
} # }
```
