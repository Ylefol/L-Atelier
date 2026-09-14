# Differential Analysis for Log2-Scale Data (limma)

Performs differential analysis using limma linear models on any
log2-scale numeric matrix. Suitable for Olink NPX data, log2-transformed
proteomics, or microarray expression. For count-based assays (RNA-seq,
ATAC-seq) use
[`ARTEMIS_differential_counts()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_differential_counts.md)
instead.

## Usage

``` r
ARTEMIS_limma_de(
  matrix,
  sample_meta,
  group_col,
  reference,
  experiment,
  sample_col = NULL,
  covariates = NULL,
  block_col = NULL,
  feature_meta = NULL,
  feature_id_col = NULL,
  fdr = 0.05,
  lfc = 0,
  verbose = TRUE
)
```

## Arguments

- matrix:

  Numeric matrix, features × samples (log2-scale). Rownames are used as
  feature identifiers and must be present.

- sample_meta:

  Data.frame with one row per sample. Rows must correspond to columns of
  `matrix`, matched by `sample_col` (or rownames if
  `sample_col = NULL`).

- group_col:

  Character. Column in `sample_meta` defining the comparison groups.

- reference:

  Character. Reference/baseline group label (denominator in fold
  change). Positive logFC = higher in `experiment` vs `reference`.

- experiment:

  Character. Experimental group label (numerator in fold change).

- sample_col:

  Character or NULL. Column in `sample_meta` holding sample IDs matching
  `colnames(matrix)`. NULL = use `rownames(sample_meta)`. Default: NULL.

- covariates:

  Character vector or NULL. Column names in `sample_meta` to include as
  additive covariates in the design matrix. Default: NULL.

- block_col:

  Character or NULL. Column in `sample_meta` identifying the blocking
  variable for repeated measures (e.g. SubjectID for longitudinal data).
  When provided,
  [`limma::duplicateCorrelation()`](https://rdrr.io/pkg/limma/man/dupcor.html)
  estimates the within-block correlation and incorporates it into the
  model. Default: NULL.

- feature_meta:

  Data.frame or NULL. Optional per-feature annotations to join onto the
  results table. Matched first by rownames (if meaningful), then by any
  column whose values match the matrix rownames. Default: NULL.

- feature_id_col:

  Character or NULL. Column in `feature_meta` to use as `feature_id` in
  the results instead of the matrix rownames. The original matrix
  rownames are retained as `original_id`. Useful for replacing internal
  IDs (e.g. OlinkID) with gene symbols (e.g. `"Assay"` from
  `ol$assay_meta`). A warning is issued if the column contains
  duplicates. Default: NULL.

- fdr:

  Numeric. FDR threshold used for the summary and the `sig` flag in
  results. Default: 0.05.

- lfc:

  Numeric. Minimum absolute log2 fold change for the `sig` flag (applied
  in addition to FDR). Default: 0.

- verbose:

  Logical. Print progress and summary. Default: TRUE.

## Value

An S3 object of class `"artemis_limma"` containing:

- results:

  Data.frame sorted by pvalue: feature_id, (original_id if
  `feature_id_col` used), (feature_meta columns if provided),
  log2FoldChange, AveExpr, t, pvalue, padj, B, sig (TRUE/FALSE).

- fit:

  The limma MArrayLM object after `eBayes()` (for advanced use).

- summary:

  List: n_features, n_tested, n_sig_up, n_sig_down, fdr, lfc.

- comparison:

  Character string "experiment vs reference".

- norm_matrix:

  The input matrix subset to the two groups analysed.

- params:

  List of parameters used (for provenance).

## Details

The design matrix uses a no-intercept parameterisation (`~ 0 + group`)
with an explicit contrast between the two groups. Covariates are added
as additive terms: `~ 0 + group + covariate1 + ...`.

For Olink NPX data, NPX is already log2-scale — do NOT transform the
matrix before passing it here. For mass spectrometry data, ensure values
are log2-transformed and normalised before calling this function.

The `sig` column combines both FDR and LFC thresholds:
`padj < fdr & abs(log2FoldChange) >= lfc`. When `lfc = 0` (default),
only FDR filtering applies.

Background for ORA/GSEA: all tested features are present in `$results`
(including non-significant). Use `results$feature_id` as the ORA
background, consistent with how
[`ARTEMIS_differential_counts()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_differential_counts.md)
is used.

## Examples

``` r
if (FALSE) { # \dontrun{
# Olink: two-group comparison
ol <- ELEUTHIA_load_olink("data.parquet", metadata_file = "layout.xlsx")
de <- ARTEMIS_limma_de(
  matrix      = ol$wide,
  sample_meta = ol$sample_meta,
  group_col   = "Group",
  reference   = "control",
  experiment  = "patient",
  feature_meta = ol$assay_meta
)
print(de)
sig_proteins <- subset(de$results, sig)

# Longitudinal design (block by subject)
de_long <- ARTEMIS_limma_de(
  matrix      = patient_matrix,
  sample_meta = patient_meta,
  group_col   = "Timepoint",
  reference   = "1",
  experiment  = "3",
  block_col   = "SubjectID",
  feature_meta = ol$assay_meta
)

# With covariate adjustment
de_adj <- ARTEMIS_limma_de(
  matrix      = ol$wide,
  sample_meta = ol$sample_meta,
  group_col   = "Group",
  reference   = "control",
  experiment  = "patient",
  covariates  = c("Age", "Sex")
)

# Use gene symbols as feature_id instead of OlinkIDs
de <- ARTEMIS_limma_de(
  matrix        = ol$wide,
  sample_meta   = ol$sample_meta,
  group_col     = "Group",
  reference     = "control",
  experiment    = "patient",
  feature_meta  = ol$assay_meta,
  feature_id_col = "Assay"
)
} # }
```
