# Artemis - Peak Annotation Distribution Comparison

Statistical comparison of genomic feature distributions between two
groups of samples, using per-sample proportions. Compare Peak Annotation
Distributions Between Two Groups

For each sample, computes the percentage of peaks falling in each
genomic feature category (Promoter, Intron, etc.), then tests whether
those percentages differ between two groups using a Wilcoxon rank-sum or
two-sample t-test. Multiple-testing correction is applied across
features.

This per-sample approach correctly treats each biological replicate as
independent, unlike pooling all peaks per group (which ignores
replication).

## Usage

``` r
ARTEMIS_compare_annotation_distribution(
  annotated_list,
  metadata,
  group_col = "group",
  feature_col = "annotation_simple",
  method = "wilcoxon",
  p_adj_method = "BH",
  verbose = TRUE
)
```

## Arguments

- annotated_list:

  Named list of data.frames from
  [`APOLLO_annotate_peaks()`](https://ylefol.github.io/L-Atelier/GAIA/reference/APOLLO_annotate_peaks.md).
  Each element corresponds to one sample; the name is used as the sample
  identifier and must match a value in `metadata$sample_id` (or
  `rownames(metadata)`).

- metadata:

  Data.frame with at least columns `sample_id` and `group_col`. If
  `sample_id` is absent, `rownames(metadata)` is used.

- group_col:

  Character. Column in `metadata` containing group labels. Exactly two
  unique values are required. Default: `"group"`.

- feature_col:

  Character. Column in each annotated data.frame to use as the feature
  category. Default: `"annotation_simple"`. If absent and `"annotation"`
  exists, a simplified version is derived automatically.

- method:

  Character. Statistical test to compare group proportions. One of
  `"wilcoxon"` (default; more robust with small n) or `"t.test"`.

- p_adj_method:

  Character. Multiple-testing correction method passed to
  [`p.adjust()`](https://rdrr.io/r/stats/p.adjust.html). Default: `"BH"`
  (Benjamini-Hochberg).

- verbose:

  Logical. Print progress. Default: `TRUE`.

## Value

A named list with class `"artemis_annotation_test"` containing:

- results:

  Data.frame of per-feature test results: feature, group mean counts,
  mean difference, test statistic, p.value, p.adj.

- counts:

  Long-format data.frame of per-sample peak counts with group labels
  attached (for plotting).

- groups:

  Character vector of the two group labels.

- group_col:

  The group column name used.

- feature_col:

  The feature column name used.

- method:

  The test method used.

## Details

**Why per-sample counts?** Each sample's per-feature peak count is
treated as the unit of inference. This preserves absolute peak numbers,
which may be biologically meaningful when total peak counts differ
between groups.

**Confounding note:** If samples differ in total peak count, raw counts
can be confounded. This is not a concern when peak libraries have been
balanced (e.g., via spike-in normalization) within each batch.

**Power note:** With few replicates (n = 2-3 per group), statistical
power is limited regardless of method. Results should be interpreted
alongside effect sizes (mean difference) rather than p-values alone.

## See also

[`AETHER_plot_annotation_comparison`](https://ylefol.github.io/L-Atelier/GAIA/reference/AETHER_plot_annotation_comparison.md)
for visualization.

## Examples

``` r
if (FALSE) { # \dontrun{
# Annotate peaks per sample
peak_list <- ELEUTHIA_load_peaks_from_sheet(sample_sheet, "ATACseq")
annotated_list <- lapply(peak_list, APOLLO_annotate_peaks, txdb = txdb)

# Build minimal metadata
meta <- sample_sheet[, c("sample_id", "group")]

# Run comparison
test_result <- ARTEMIS_compare_annotation_distribution(
  annotated_list = annotated_list,
  metadata       = meta,
  group_col      = "group"
)

# Visualise
p <- AETHER_plot_annotation_comparison(test_result)
} # }
```
