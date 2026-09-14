# Time Series Conditional Differential Expression (limma)

Compares two groups at each timepoint independently using limma linear
models. For example, Treatment vs Control at TP1, TP2, TP3. Operates on
the log2-scale matrix stored in an `artemis_ts_norm_limma` object — no
transformation is applied during the comparison.

## Usage

``` r
ARTEMIS_timeseries_conditional_limma(
  ts_norm,
  reference,
  experiment,
  block_col = NULL,
  covariates = NULL,
  alpha = 0.05,
  verbose = TRUE
)
```

## Arguments

- ts_norm:

  An `artemis_ts_norm_limma` object from
  [`ARTEMIS_normalize_timeseries_limma()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_normalize_timeseries_limma.md).

- reference:

  Character. Reference/baseline group (denominator in FC).

- experiment:

  Character. Experimental group (numerator in FC).

- block_col:

  Character or NULL. Column in targets for the blocking variable
  (repeated measures within each timepoint comparison, e.g. a crossover
  design). When provided,
  [`limma::duplicateCorrelation()`](https://rdrr.io/pkg/limma/man/dupcor.html)
  estimates the within-block correlation per timepoint. Default: NULL.

- covariates:

  Character vector or NULL. Additional covariate columns from targets to
  include in the design. If `batch_col` was set in
  [`ARTEMIS_normalize_timeseries_limma()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_normalize_timeseries_limma.md),
  it is automatically prepended. Default: NULL.

- alpha:

  Numeric. FDR threshold for counting significant features. Default:
  0.05.

- verbose:

  Logical. Print progress. Default: TRUE.

## Value

An S3 object of class `"artemis_ts_de"` (same class as the DESeq2
workflow) containing:

- results:

  Named list of per-timepoint result lists, each with `$results`
  (data.frame: feature_id, log2FoldChange, AveExpr, t, pvalue, padj, B,
  sig), `$fit`, `$summary`

- summary:

  Data.frame: timepoint, experiment_name, n_tested, n_sig_up, n_sig_down

- type:

  "conditional"

- comparison:

  Character "experiment vs reference"

- norm_counts:

  The log2 matrix (named for compatibility with
  [`ARTEMIS_prepare_part_matrix()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_prepare_part_matrix.md)
  and
  [`ARTEMIS_select_de_genes()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_select_de_genes.md))

- parameters:

  List of parameters used

## Examples

``` r
if (FALSE) { # \dontrun{
ts_norm <- ARTEMIS_normalize_timeseries_limma(ol$wide, ol$sample_meta,
  group_col = "Group", time_col = "Timepoint")

cond_de <- ARTEMIS_timeseries_conditional_limma(
  ts_norm, reference = "Control", experiment = "Sepsis"
)
} # }
```
