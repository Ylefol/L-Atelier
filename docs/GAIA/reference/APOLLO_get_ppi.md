# Get protein-protein interaction network from OmniPath

Retrieves PPI data from OmniPath (aggregating STRING, BioGRID, IntAct,
and other sources) for a set of genes and returns an igraph network
object.

## Usage

``` r
APOLLO_get_ppi(
  genes,
  organism = 9606,
  resources = NULL,
  min_resources = 1,
  drop_isolates = FALSE,
  directed = FALSE,
  verbose = TRUE
)
```

## Arguments

- genes:

  Character vector of gene symbols to include in the network.

- organism:

  Integer. NCBI taxonomy ID. Default: 9606 (human). Common values: 9606
  (human), 10090 (mouse), 10116 (rat).

- resources:

  Character vector. Specific interaction databases to query. Default:
  NULL (all available). Examples: "STRING", "BioGRID", "IntAct",
  "SIGNOR", "PhosphoSite".

- min_resources:

  Integer. Minimum number of databases supporting an interaction for it
  to be included. Default: 1. Higher values = more stringent filtering.

- drop_isolates:

  Logical. Remove query genes with zero interactions (degree 0) from the
  returned graph, instead of keeping them as disconnected nodes.
  Default: FALSE.

- directed:

  Logical. Whether to return directed interactions. Default: FALSE
  (undirected PPI network).

- verbose:

  Logical. Print network summary. Default: TRUE.

## Value

An igraph graph object with:

- Vertex attributes:

  name (gene symbol), degree, betweenness

- Edge attributes:

  source, target, n_resources, resources (database names)

- Graph attributes:

  organism, n_query_genes, n_mapped_genes

Also has attribute "interaction_df" with the raw interaction data.frame.

## Details

OmniPath aggregates protein-protein interactions from dozens of
databases including STRING, BioGRID, IntAct, SIGNOR, PhosphoSite, and
more. This provides broader coverage than any single database.

Only interactions where BOTH endpoints are in the provided gene list are
returned (within-list network). To include first neighbors, add them to
the gene list before calling.

Network metrics (degree, betweenness) are computed automatically and
stored as vertex attributes for use in downstream visualization.

## Examples

``` r
if (FALSE) { # \dontrun{
# PPI network for a set of DEGs
graph <- APOLLO_get_ppi(c("TP53", "BRCA1", "EGFR", "MYC", "CDK2"))

# Mouse network from specific databases
graph <- APOLLO_get_ppi(genes, organism = 10090, resources = c("STRING", "BioGRID"))

# Stricter filtering: require at least 2 supporting databases
graph <- APOLLO_get_ppi(genes, min_resources = 2)

} # }
```
