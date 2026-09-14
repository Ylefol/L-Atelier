# Normalize Count Data with DESeq2

Creates a DESeq2 object with size factors estimated on the FULL dataset.
Use this when performing multiple pairwise comparisons from the same
experiment - normalizing once ensures consistent size factors across all
subsequent comparisons.

## Usage

``` r
ARTEMIS_normalize_counts(
  counts,
  targets,
  group_col = "group",
  batch_col = NULL,
  size_factors = NULL,
  verbose = TRUE
)
```

## Arguments

- counts:

  Integer matrix of raw counts (features x samples). Rownames required.

- targets:

  Data.frame with sample metadata. Rownames must match colnames of
  counts. Must include a group column.

- group_col:

  Character. Column in targets for group labels. Default: "group".

- batch_col:

  Character or NULL. Optional column name for a single blocking variable
  to include in the design formula as `~ batch_col + condition`.
  Default: NULL. **Limitation:** only one blocking variable is
  supported; multi-variable designs are not yet implemented.

- size_factors:

  NULL, a named numeric vector, or a single numeric scalar. Controls how
  per-sample size factors are set. Default: NULL.

  NULL

  :   DESeq2 estimates size factors itself via `estimateSizeFactors()`
      (median-of-ratios). This is correct when the input counts are raw
      and untouched by any prior normalization.

  Named numeric vector

  :   Names must match `rownames(targets)` (equivalently
      `colnames(counts)`); every sample must have a value. Assigned
      directly via `sizeFactors(dds) <-`, skipping DESeq2's own
      estimation entirely. Use this to hand DESeq2 an externally
      computed scale, e.g. spike-in-derived factors from
      `HORIZON_compute_spike_in_factors()`'s `size_factor_deseq2`
      column.

  Single scalar

  :   Recycled across every sample (e.g. `1`). Use this to disable
      DESeq2-level normalization entirely – appropriate when the input
      counts have already been put on a comparable scale upstream (e.g.
      by physical read downsampling) and no further per-sample
      correction should be applied.

  All values must be strictly positive (DESeq2's own requirement).

- verbose:

  Logical. Print progress. Default: TRUE.

## Value

An S3 object of class `"artemis_norm"` containing:

- dds:

  DESeqDataSet with size factors set on full dataset

- norm_counts:

  Normalized count matrix

- targets:

  The targets data.frame

- size_factors:

  Named vector of size factors per sample

- parameters:

  List of parameters used, including `size_factors_source`
  ("deseq2_estimated" or "external")

## Details

This function:

1.  Creates a DESeqDataSet from the full count matrix

2.  Sets size factors using all samples together – either estimated by
    DESeq2 or supplied via `size_factors`

3.  Stores the normalized counts and size factors

The resulting object can be passed to
[`ARTEMIS_differential_counts()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_differential_counts.md)
which will subset samples for each comparison while preserving the
global size factors.

## Why normalize first

When performing multiple DEA comparisons from the same experiment (e.g.,
A vs B, A vs C, B vs C), normalizing each comparison independently can
introduce inconsistencies. Estimating size factors once on all samples
ensures that the normalized expression values are comparable across all
comparisons.

## Externally supplied size factors

DESeq2's default (median-of-ratios) normalization assumes most features
are non-differential between samples and infers relative library scale
from the count matrix itself. When an experiment includes a spike-in
control (ATAC, CHIP/CUT&RUN, R-loop), that assumption can be exactly
what you don't want – a real, global shift in signal violates it.
`size_factors` lets you substitute a ground-truth scale (from the
spike-in) in place of DESeq2's own estimate. Do not combine externally
supplied size factors with counts that have *already* been normalized
upstream (e.g. via physical BAM downsampling to a spike-in ratio) – that
applies the same correction twice. Pick one: raw counts + `size_factors`
(spike-in derived), or already-normalized counts + `size_factors = 1`
(disables further correction).

## Examples

``` r
if (FALSE) { # \dontrun{
# Step 1: Normalize the full dataset
norm_data <- ARTEMIS_normalize_counts(counts, sample_sheet)

# Step 2: Run multiple DEA comparisons (all use same size factors)
dea_A_vs_B <- ARTEMIS_differential_counts(norm_data, "B", "A")
dea_C_vs_B <- ARTEMIS_differential_counts(norm_data, "B", "C")

# Step 3: Combine results for downstream analysis
all_genes <- ARTEMIS_select_de_genes(list(dea_A_vs_B, dea_C_vs_B))

# Spike-in normalized ATAC/CHIP: supply factors instead of letting DESeq2
# estimate its own (raw, non-downsampled counts required for this to be a
# single, correct normalization rather than a double one).
spikein_factors <- HORIZON_compute_spike_in_factors(spikein_bams)
sf <- setNames(spikein_factors$size_factor_deseq2, spikein_factors$sample_id)
norm_data <- ARTEMIS_normalize_counts(counts, sample_sheet, size_factors = sf)

# Counts already normalized upstream -- disable DESeq2-level normalization
norm_data <- ARTEMIS_normalize_counts(counts, sample_sheet, size_factors = 1)
} # }
```
