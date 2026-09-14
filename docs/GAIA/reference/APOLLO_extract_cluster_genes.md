# Extract gene lists from clustering results

Extracts gene lists from clustering results, returning a named list
suitable for enrichment analysis. Accepts `wgcna_modules`, `wgcna_hubs`
(from ARTEMIS_wgcna\_\*), or `artemis_part` objects. When passed a
wgcna_hubs object (from ARTEMIS_wgcna_hub_genes() with n_top = NULL and
a kME threshold), only hub genes passing the threshold are returned —
suitable for enrichment analysis on large modules.

## Usage

``` r
APOLLO_extract_cluster_genes(
  modules,
  modules_of_interest = NULL,
  exclude_grey = TRUE,
  strip_version = FALSE,
  org_db = NULL,
  from_type = "ENSEMBL",
  to_type = "ENTREZID",
  verbose = TRUE
)
```

## Arguments

- modules:

  A `wgcna_modules` object from
  [`ARTEMIS_wgcna_detect_modules()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_wgcna_detect_modules.md),
  a `wgcna_hubs` object from
  [`ARTEMIS_wgcna_hub_genes()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_wgcna_hub_genes.md),
  or an `artemis_part` object from
  [`ARTEMIS_part()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_part.md).

- modules_of_interest:

  Character vector of module/cluster names to extract. Default: NULL
  (all clusters, excluding unassigned).

- exclude_grey:

  Logical. Exclude unassigned genes: the grey module for WGCNA
  (`module_0`), or the outlier cluster for PART (`C0`). Default: TRUE.

- strip_version:

  Logical. Remove version numbers from Ensembl-style IDs (e.g.,
  ENSG00000141510.16 -\> ENSG00000141510). Default: FALSE.

- org_db:

  Optional OrgDb for ID conversion. If provided, converts to ENTREZID.

- from_type:

  Type of input gene IDs (for conversion). Default: "ENSEMBL".

- to_type:

  Type of output gene IDs (for conversion). Default: "ENTREZID".

- verbose:

  Logical. Print progress. Default: TRUE.

## Value

Named list of gene vectors, one per cluster/module.

## Examples

``` r
if (FALSE) { # \dontrun{
# From WGCNA modules
cluster_genes <- APOLLO_extract_cluster_genes(modules)

# From PART clustering (C0 outliers excluded by default)
part_result <- ARTEMIS_part(mat, seed = 42)
cluster_genes <- APOLLO_extract_cluster_genes(part_result)

# Strip Ensembl version numbers for gprofiler
cluster_genes <- APOLLO_extract_cluster_genes(modules, strip_version = TRUE)

# Get ENTREZID from Ensembl IDs
library(org.Hs.eg.db)
cluster_entrez <- APOLLO_extract_cluster_genes(modules, strip_version = TRUE,
                                                org_db = org.Hs.eg.db,
                                                from_type = "ENSEMBL")

} # }
```
