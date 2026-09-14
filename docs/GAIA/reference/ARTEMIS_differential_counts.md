# Differential Analysis for Count Data

Performs differential analysis on count data using DESeq2. Requires an
`artemis_norm` object from
[`ARTEMIS_normalize_counts()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_normalize_counts.md)
as input, ensuring size factors and batch correction are set
consistently across all comparisons before any DEA is run.

## Usage

``` r
ARTEMIS_differential_counts(
  norm_data,
  reference,
  experiment,
  alpha = 0.05,
  verbose = TRUE
)
```

## Arguments

- norm_data:

  An `artemis_norm` object from
  [`ARTEMIS_normalize_counts()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_normalize_counts.md).
  Size factors and batch correction settings are inherited from this
  object.

- reference:

  Character. The reference/baseline group label (e.g., "WT"). This is
  the denominator in fold change calculations.

- experiment:

  Character. The experimental/comparison group label (e.g., "KO"). This
  is the numerator in fold change calculations.

- alpha:

  Numeric. FDR threshold for summary statistics (default = 0.05).

- verbose:

  Logical. Print progress and summary (default = TRUE).

## Value

An `artemis_dea` object (list) containing:

- results:

  Data.frame with per-feature DESeq2 results: feature_id, baseMean,
  log2FoldChange, lfcSE, stat, pvalue, padj

- dds:

  The DESeqDataSet object for further analysis

- summary:

  List with summary statistics: n_features, n_sig_up, n_sig_down,
  n_tested, alpha

- comparison:

  Character describing the comparison (experiment vs reference)

- norm_counts:

  Normalized count matrix (useful for downstream analysis)

Calling [`print()`](https://rdrr.io/r/base/print.html) on the returned
object reprints the comparison and significance summary (and top 5 hits)
without rerunning DESeq2 – see
[`print.artemis_dea`](https://ylefol.github.io/L-Atelier/GAIA/reference/print.artemis_dea.md).

## Details

Always call
[`ARTEMIS_normalize_counts()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_normalize_counts.md)
on the full dataset first, then pass the result here for each pairwise
comparison. This ensures:

- Size factors are estimated on all samples, not just the two being
  compared

- Batch correction is configured in one place and applied consistently

- Normalized counts are available for visualization before running DEA

The function subsets to the two groups of interest, preserves the
pre-computed size factors, estimates dispersions, and runs the Wald
test. Positive log2FC means higher in experiment relative to reference.

## Examples

``` r
if (FALSE) { # \dontrun{
# Step 1: Normalize full dataset (with optional batch correction)
norm_data <- ARTEMIS_normalize_counts(counts, sample_sheet, batch_col = "batch")

# Step 2: Run each pairwise comparison
dea_KO_vs_WT  <- ARTEMIS_differential_counts(norm_data, "WT", "KO")
dea_OE_vs_WT  <- ARTEMIS_differential_counts(norm_data, "WT", "OE")

# View significant features
sig_features <- subset(dea_KO_vs_WT$results, padj < 0.05)
} # }
```
