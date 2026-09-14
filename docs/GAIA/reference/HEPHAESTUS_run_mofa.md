# Multi-Omics Factor Analysis via MOFA2

Trains a MOFA2 model on multiple omics views, following the same
list-of-matrices convention as
[`HEPHAESTUS_multi_convCCA()`](https://ylefol.github.io/L-Atelier/GAIA/reference/HEPHAESTUS_multi_convCCA.md)
– but unlike the sparse-CCA methods in this module, MOFA2 natively
tolerates samples missing entirely from one or more views. Internally
the views are converted to MOFA2's long-format input
(`sample, feature, view, value`) rather than its stricter matrix-input
path (which requires identical, identically-ordered sample sets across
every view), so callers can pass views with only partial sample overlap
without manually aligning/padding them first.

## Usage

``` r
HEPHAESTUS_run_mofa(
  X,
  groups = NULL,
  num_factors = 10,
  likelihoods = NULL,
  spikeslab_weights = TRUE,
  spikeslab_factors = FALSE,
  scale_views = FALSE,
  convergence_mode = "fast",
  maxiter = 1000,
  drop_factor_threshold = -1,
  seed = 42,
  use_basilisk = TRUE,
  outfile = NULL,
  verbose = TRUE
)
```

## Arguments

- X:

  Named list of matrices, one per view, samples in rows and features in
  columns (GAIA's standard convention, same as
  [`HEPHAESTUS_multi_convCCA()`](https://ylefol.github.io/L-Atelier/GAIA/reference/HEPHAESTUS_multi_convCCA.md)).
  List names become MOFA2 view names. Views do not need the same samples
  or sample order – rownames are the join key across views.

- groups:

  Optional named character vector or factor, sample ID -\> group label,
  for MOFA2's multi-group framework. Default `NULL` (single group).
  MOFA2's own documentation describes multi-group as an advanced option
  not recommended as a first approach – treat this as experimental.

- num_factors:

  Integer. Number of latent factors to fit. Default `10` (MOFA2's own
  package default). For views with few samples, check
  `variance_explained` in the result rather than assuming 10 is
  appropriate – this default is not tuned to any particular dataset.

- likelihoods:

  Optional named character vector, view name -\> `"gaussian"`,
  `"poisson"`, or `"bernoulli"`. Default `NULL` uses `"gaussian"` for
  every view. MOFA2's own guidance is to prefer transforming
  count/binary data to continuous values rather than relying on
  non-gaussian likelihoods, which give non-optimal results.

- spikeslab_weights:

  Logical. Use a spike-and-slab sparsity prior on factor weights.
  Default `TRUE` (MOFA2 default).

- spikeslab_factors:

  Logical. Use a spike-and-slab sparsity prior on factors. Default
  `FALSE` (MOFA2 default).

- scale_views:

  Logical. Scale each view to unit variance before training. Default
  `FALSE`.

- convergence_mode:

  Character. `"fast"`, `"medium"`, or `"slow"` – controls the
  ELBO-change convergence threshold. Default `"fast"`.

- maxiter:

  Integer. Maximum training iterations. Default `1000`.

- drop_factor_threshold:

  Numeric. Fraction of variance explained below which a factor is
  dropped during training; `-1` disables dropping. Default `-1` (MOFA2
  default).

- seed:

  Integer. Random seed for reproducibility. Default `42` (MOFA2
  default).

- use_basilisk:

  Logical. Let MOFA2 provision its own basilisk-managed Python
  environment for its `mofapy2` backend. Default `TRUE` – unlike GAIA's
  other basilisk-wrapped tools (CPAT, SignalP, Pfam, HOMER), no
  `.gaia_mofa_env` exists in `GAIA/R/basilisk.R`; MOFA2 manages this
  itself. Set `FALSE` only if you have manually configured a
  Python/`mofapy2` installation via
  [`reticulate::use_python()`](https://rstudio.github.io/reticulate/reference/use_python.html).

- outfile:

  Character or `NULL`. Path to save the trained model (`.hdf5`). Default
  `NULL` uses a temp file – pass an explicit path if you intend to
  reload the model later with
  [`HEPHAESTUS_load_mofa_model()`](https://ylefol.github.io/L-Atelier/GAIA/reference/HEPHAESTUS_load_mofa_model.md),
  since temp files are not guaranteed to persist across R sessions.

- verbose:

  Logical. Print progress. Default `TRUE`.

## Value

A classed `hephaestus_mofa` list:

- model:

  The raw trained `MOFA` object, for advanced use.

- factors:

  Tidy data.frame of factor scores (sample, group, factor, value).

- weights:

  Tidy data.frame of feature weights (view, feature, factor, value).

- variance_explained:

  Tidy data.frame of variance explained (group, view, factor,
  variance_explained).

- params:

  List of the training/model/data options actually used.

## Details

Field names for `model_options`/`data_options`/ `training_options` are
set from MOFA2's own
`get_default_model_options()`/`get_default_data_options()`/
`get_default_training_options()` and were not independently pinned
against a citable spec beyond MOFA2's own documentation – if a field
name changes in a future MOFA2 release, the option-setting step below
will fail loudly (unknown list element) rather than silently doing the
wrong thing.

This function is intentionally not wired into
[`MINERVA_scca_fit()`](https://ylefol.github.io/L-Atelier/GAIA/reference/MINERVA_scca_fit.md)/
[`MINERVA_cv_scca()`](https://ylefol.github.io/L-Atelier/GAIA/reference/MINERVA_cv_scca.md):
those normalize every method to a `list(W = list(...))`
canonical-weight-vector contract and search sparsity parameters via
k-fold CV, neither of which applies to a shared latent-factor model
selected by number-of-factors/ELBO.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(123)
# view_a and view_b share all 50 samples; view_c only covers a 15-sample
# subset -- MOFA2 handles this partial overlap natively.
view_a <- matrix(rnorm(50 * 100), 50, 100,
                  dimnames = list(paste0("S", 1:50), paste0("gA", 1:100)))
view_b <- matrix(rnorm(50 * 80), 50, 80,
                  dimnames = list(paste0("S", 1:50), paste0("gB", 1:80)))
view_c <- matrix(rnorm(15 * 60), 15, 60,
                  dimnames = list(paste0("S", 1:15), paste0("gC", 1:60)))

result <- HEPHAESTUS_run_mofa(
  X = list(rna = view_a, atac = view_b, cutntag = view_c),
  num_factors = 5
)
print(result)
head(result$factors)

} # }
```
