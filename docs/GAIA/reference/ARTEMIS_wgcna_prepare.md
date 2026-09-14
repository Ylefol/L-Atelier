# Prepare data for WGCNA analysis

Formats expression data and traits for WGCNA. Accepts either pre-loaded
data objects (from ELEUTHIA functions) or file paths.

## Usage

``` r
ARTEMIS_wgcna_prepare(
  counts,
  traits,
  counts_format = "auto",
  sample_col = NULL,
  gene_col = 1,
  transpose = TRUE,
  remove_zero_variance = TRUE,
  verbose = TRUE
)
```

## Arguments

- counts:

  Expression matrix/data.frame OR path to count file/directory. If
  matrix/data.frame: genes as rows, samples as columns (will be
  transposed). If path: will be loaded based on counts_format parameter.

- traits:

  Data.frame OR path to traits file. Must have sample identifiers that
  match count data.

- counts_format:

  How to interpret counts if it's a path:

  - "file": single CSV/TSV file with genes as rows

  - "directory": folder of individual count files (one per sample)

  - "auto": attempt to detect (default when counts is a path) Ignored if
    counts is already a matrix/data.frame.

- sample_col:

  Column name in traits containing sample identifiers. If NULL, uses
  rownames of traits. Default: NULL.

- gene_col:

  Column name/index for gene identifiers in count data. Default: 1
  (first column) or rownames.

- transpose:

  Logical. If TRUE, transpose counts so samples are rows. Default: TRUE
  (standard gene x sample input becomes sample x gene).

- remove_zero_variance:

  Logical. Remove genes with zero variance. Default: TRUE.

- verbose:

  Logical. Print progress messages. Default: TRUE.

## Value

A list with class "wgcna_data" containing:

- datExpr:

  Expression matrix with samples as rows, genes as columns

- datTraits:

  Traits data.frame with samples as rows, matching datExpr

- gene_names:

  Character vector of gene names

- sample_names:

  Character vector of sample names

- n_genes:

  Number of genes

- n_samples:

  Number of samples

- removed_genes:

  Genes removed due to zero variance (if any)

- removed_samples:

  Samples removed due to missing data (if any)

## Examples

``` r
if (FALSE) { # \dontrun{
# From pre-loaded ELEUTHIA data
counts <- ELEUTHIA_load_rnaseq_from_sheet(sample_sheet)
traits <- sample_sheet[, c("sample_id", "age", "condition")]
wgcna_data <- ARTEMIS_wgcna_prepare(counts, traits, sample_col = "sample_id")

# From file paths
wgcna_data <- ARTEMIS_wgcna_prepare(
  counts = "path/to/counts.csv",
  traits = "path/to/metadata.csv",
  sample_col = "sample_id"
)

} # }
```
