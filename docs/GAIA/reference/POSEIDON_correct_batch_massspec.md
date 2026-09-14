# Batch Correction for Mass Spectrometry Data

Removes batch effects from a `massspec_data` object's log2 intensity
matrix (`$wide`) using ComBat or
[`limma::removeBatchEffect`](https://rdrr.io/pkg/limma/man/removeBatchEffect.html).
Unlike
[`POSEIDON_correct_batch_counts()`](https://ylefol.github.io/L-Atelier/GAIA/reference/POSEIDON_correct_batch_counts.md),
this operates on data that is already log2-transformed and contains NAs
– no log/back-transform is applied, and missingness is handled via a
disposable imputed working copy used only to make the underlying
algorithm runnable.

## Usage

``` r
POSEIDON_correct_batch_massspec(
  massspec_data,
  batch_col = "PlateID",
  group_col = NULL,
  method = "combat",
  impute_method = "batchmean",
  verbose = TRUE
)
```

## Arguments

- massspec_data:

  A `massspec_data` object (from
  [`ELEUTHIA_load_massspec()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ELEUTHIA_load_massspec.md)),
  typically after
  [`HADES_filter_massspec()`](https://ylefol.github.io/L-Atelier/GAIA/reference/HADES_filter_massspec.md)/[`HADES_filter_massspec_proteins()`](https://ylefol.github.io/L-Atelier/GAIA/reference/HADES_filter_massspec_proteins.md)
  and
  [`POSEIDON_normalize_massspec()`](https://ylefol.github.io/L-Atelier/GAIA/reference/POSEIDON_normalize_massspec.md).

- batch_col:

  Character string. Column in `$sample_meta` containing batch
  information (default = `"PlateID"`).

- group_col:

  Character string or NULL. Optional column in `$sample_meta` for a
  biological covariate to preserve during correction (default = NULL, no
  covariate preserved). Samples with NA in this column (e.g. POOL
  samples, which have no biological group) are assigned their own
  `"no_group"` level rather than being dropped, so every sample is still
  corrected.

- method:

  Character string. `"combat"` (default) or `"limma"`.

- impute_method:

  Character string. How the disposable working copy used to make
  ComBat/limma runnable is filled in:

  - `"batchmean"` (default) – each missing cell is filled with its
    protein's mean WITHIN its own batch (observed cells only). Falls
    back to the protein's global mean for any protein/batch combination
    with zero observed values in that batch.

  - `"global"` – each missing cell is filled with its protein's mean
    across ALL samples, regardless of batch (the original behavior of
    this function).

  Comparing the two on the SEPSOMICS dataset (see
  `projects/SEPSOMICS/massSpec_imputation_comparison.R`) showed
  `"global"` under-corrects proteins with high missingness in a given
  batch – the placeholder pulls that batch's estimate toward the grand
  mean instead of its own, diluting the estimated batch effect.
  `"batchmean"` does not have this failure mode, but that was
  established on one dataset/filtering setup; `"global"` is kept
  available in case a different experiment's missingness pattern behaves
  differently.

- verbose:

  Logical. Print progress messages (default = TRUE).

## Value

The input `massspec_data` object with `$wide` batch-corrected. Values
that were NA before correction remain NA in the output – imputation is
used only internally to make ComBat/limma runnable; the imputed values
themselves are never written back.

## Details

DIA mass spec retains substantial missingness even after filtering, and
both ComBat and
[`limma::removeBatchEffect`](https://rdrr.io/pkg/limma/man/removeBatchEffect.html)
require a complete matrix. A temporary copy of `$wide` is imputed (per
`impute_method`) purely so the algorithm can run; the correction is
computed on that copy, but only the originally non-NA cells of the
result are written back. Cells that were NA before correction stay NA.

## Examples

``` r
if (FALSE) { # \dontrun{
ms <- HADES_filter_massspec(ms_raw, keep_types = c("SAMPLE", "POOL"))
ms <- HADES_filter_massspec_proteins(ms, max_na_fraction = 0.5)
ms$wide <- POSEIDON_normalize_massspec(ms$wide)
ms <- POSEIDON_correct_batch_massspec(ms, batch_col = "PlateID")

# Preserve a biological covariate during correction
ms <- POSEIDON_correct_batch_massspec(ms, batch_col = "PlateID", group_col = "Group")

# Use the original global-mean imputation instead
ms <- POSEIDON_correct_batch_massspec(ms, batch_col = "PlateID", impute_method = "global")
} # }
```
