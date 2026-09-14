# Project New Samples onto a Trained MOFA2 Model

MOFA2 has no out-of-sample projection of its own:
[`predict()`](https://rdrr.io/r/stats/predict.html)/`impute()` only
reconstruct the model's own training data via its already-fitted factor
scores, and `interpolate_factors()` is MEFISTO-only (interpolates along
a trained covariate for existing training samples, not new ones). This
fills that gap: holding the trained per-view feature weights fixed, it
solves for each new sample's factor vector by ordinary least squares
against `Y_view = Z %*% t(W_view)`, after centering the new sample's
data the same way MOFA centered the training data internally
(per-feature training mean; confirmed empirically – reconstructing
training data as `Z %*% t(W)` without adding the feature mean back
leaves residuals on the order of the feature means themselves, and
centered residuals are small). A view a sample has no data for simply
drops out of that sample's fit – no imputation is attempted.

## Usage

``` r
HEPHAESTUS_project_mofa_samples(
  mofa_result,
  train_data,
  new_data,
  verbose = TRUE
)
```

## Arguments

- mofa_result:

  A `hephaestus_mofa` object from
  [`HEPHAESTUS_run_mofa()`](https://ylefol.github.io/L-Atelier/GAIA/reference/HEPHAESTUS_run_mofa.md)
  or
  [`HEPHAESTUS_load_mofa_model()`](https://ylefol.github.io/L-Atelier/GAIA/reference/HEPHAESTUS_load_mofa_model.md).

- train_data:

  The named list of view matrices (samples x features) originally passed
  as `X` to
  [`HEPHAESTUS_run_mofa()`](https://ylefol.github.io/L-Atelier/GAIA/reference/HEPHAESTUS_run_mofa.md)
  to train this model. Used only to recover each view's per-feature
  training mean for centering – MOFA2 does not expose this on the fitted
  model object.

- new_data:

  Named list of view matrices (samples x features) for the new samples
  to project. Must use the SAME normalization/transform as `train_data`
  (e.g. the same VST dispersion fit, the same batch correction) – this
  function does not check for that, only for shared feature names.
  Feature sets are intersected against `train_data` and the model's
  weights per view; extra features are silently dropped. Views may cover
  different, partially-overlapping sets of samples, same as MOFA2
  tolerates for training.

- verbose:

  Logical. Print progress. Default `TRUE`.

## Value

A matrix (new samples x factors) of projected factor scores, in the same
units as `mofa_result$factors`.

## Details

Ordinary least squares, solved independently per sample via
[`qr.solve()`](https://rdrr.io/r/base/qr.html) on that sample's stacked
available-view design. No regularization is applied – if a sample's
total usable feature count (summed across its available views) is small
relative to the number of factors, the fit is underdetermined/unstable;
this is a known limitation, not handled automatically. Samples with zero
usable features in every view are returned as all-NA rows, with a
warning.
