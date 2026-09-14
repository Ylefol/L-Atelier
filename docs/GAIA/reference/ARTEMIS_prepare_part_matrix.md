# Prepare Expression Matrix for PART Clustering

Extracts normalized counts from DEA results and subsets to selected
genes, optionally z-score scaling each gene. This bridges DEA results to
PART input.

## Usage

``` r
ARTEMIS_prepare_part_matrix(
  de_result,
  genes,
  samples = NULL,
  scale = TRUE,
  log_transform = FALSE,
  verbose = TRUE
)
```

## Arguments

- de_result:

  DEA result in one of the following formats:

  - DEA result from
    [`ARTEMIS_differential_counts()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_differential_counts.md)
    (extracts normalized counts from the DESeqDataSet)

  - `artemis_ts_de` object (extracts norm_counts)

  - A numeric matrix of normalized expression (genes x samples)

- genes:

  Character vector of gene IDs to include (e.g., from
  [`ARTEMIS_select_de_genes()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_select_de_genes.md)).

- samples:

  Character vector of sample IDs (column names) to include, in the
  desired order. NULL = all samples. Default: NULL.

- scale:

  Logical. Z-score scale each gene (row) across samples. Recommended for
  PART when using euclidean distance. Default: TRUE.

- log_transform:

  Logical. Apply log2(x + 1) transformation before scaling. Default:
  FALSE (matches original TiSA behavior). Set TRUE if you want log-scale
  distances for normalized counts.

- verbose:

  Logical. Print progress. Default: TRUE.

## Value

Numeric matrix ready for
[`ARTEMIS_part()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_part.md).

## Examples

``` r
if (FALSE) { # \dontrun{
# From standard DEA result (default: no log transform, matches TiSA)
dea <- ARTEMIS_differential_counts(quant, "WT", "KO")
genes <- ARTEMIS_select_de_genes(dea, l2fc_thresh = 1)
part_mat <- ARTEMIS_prepare_part_matrix(dea, genes)
part_result <- ARTEMIS_part(part_mat, seed = 42)

# From time series DEA
genes <- ARTEMIS_select_de_genes(cond_de, l2fc_thresh = 2)
part_mat <- ARTEMIS_prepare_part_matrix(cond_de, genes)

# With log transformation (optional, if you prefer log-scale distances)
part_mat <- ARTEMIS_prepare_part_matrix(dea, genes, log_transform = TRUE)

} # }
```
