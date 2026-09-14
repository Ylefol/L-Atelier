# Load and combine DEA gene lists from multiple experiments

Scans a parent directory for DEA result subdirectories and collects
significant gene IDs across experiments. Useful for building a combined
gene selection for downstream analyses (e.g., PART clustering).

## Usage

``` r
DEMETER_load_dea_genes(
  dea_dir,
  source = "significant",
  prefix = "dea",
  experiments = NULL,
  l2fc_thresh = 1,
  p_thresh = 0.05,
  gene_col = NULL,
  verbose = TRUE
)
```

## Arguments

- dea_dir:

  Path to the parent directory containing DEA result subdirectories (one
  per experiment/comparison).

- source:

  Character. Which file to read genes from:

  - "significant" (default): reads the pre-filtered significant results
    CSV

  - "selected": reads the selected genes CSV (gene_id column)

  - "all": reads the full results CSV and applies threshold filtering

- prefix:

  Character. Filename prefix used during export. Default: "dea".

- experiments:

  Character vector. Specific subdirectory names to include. If NULL
  (default), all subdirectories containing the expected files are used.

- l2fc_thresh:

  Numeric. Absolute log2 fold-change threshold (only used when source =
  "all"). Default: 1.0.

- p_thresh:

  Numeric. Adjusted p-value threshold (only used when source = "all").
  Default: 0.05.

- gene_col:

  Character. Column name containing gene IDs. If NULL (default),
  auto-detected from: feature_id, gene_id, gene, ensembl_id.

- verbose:

  Logical. Print progress messages. Default: TRUE.

## Value

Character vector of unique gene IDs (union across all experiments).
Attributes:

- "per_experiment": named list of gene vectors per experiment

- "n_per_experiment": named integer vector of gene counts

- "source": which file source was used

- "parameters": list of thresholds used (if source = "all")

## Examples

``` r
if (FALSE) { # \dontrun{
# Default: collect from pre-filtered significant files
genes <- DEMETER_load_dea_genes("results/DEA/")

# Only specific experiments
genes <- DEMETER_load_dea_genes("results/DEA/", experiments = c("C5RO", "CS1AN"))

# Custom thresholds on all results
genes <- DEMETER_load_dea_genes("results/DEA/", source = "all",
                                 l2fc_thresh = 1.5, p_thresh = 0.01)

# Use as gene selection for PART
part_mat <- ARTEMIS_prepare_part_matrix(norm_data, genes = genes)

} # }
```
