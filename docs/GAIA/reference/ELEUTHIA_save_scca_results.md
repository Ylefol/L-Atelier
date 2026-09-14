# Eleuthia - Save Results Functions

Functions for saving analysis results, plots, and summaries in an
organized structure. Save sCCA Analysis Results

Comprehensive save function for sCCA cross-validation analysis. Saves
plots (PNG and RDS), CV results, best model, parameter grids, and
generates a plain text summary.

## Usage

``` r
ELEUTHIA_save_scca_results(
  output_dir,
  coarse_cv_results,
  fine_cv_results,
  plots,
  method,
  X,
  pilot_results = NULL,
  prefix = "scca",
  metadata = list()
)
```

## Arguments

- output_dir:

  Character string. Directory where all results will be saved.

- coarse_cv_results:

  List. Output from MINERVA_cv_scca for coarse grid. Must contain
  \$cv_results and \$best_model.

- fine_cv_results:

  List. Output from MINERVA_cv_scca for fine grid. Must contain
  \$cv_results and \$best_model.

- plots:

  List of ggplot objects with named elements:

  - cv_marginal_full: Marginal plot with full x-axis

  - cv_marginal_zoom: Marginal plot zoomed to data (optional)

  - cv_best_scores: Bar plot of best CV scores

  - feature_weights: Multi-feature weight plot

- method:

  Character string. Method name ("ConvCCA" or "RelPMDCCA").

- X:

  List of datasets used in analysis (for metadata).

- pilot_results:

  List. Output from MINERVA_pilot_compare_methods (optional).

- prefix:

  Character string. Prefix for saved files (default = "scca").

- metadata:

  List. Additional metadata to include in summary:

  - k: Number of CV folds

  - lambda: Lambda parameter

  - nIter: Number of iterations

  - penalty: Penalty type

  - metric: CV metric used

## Value

Invisibly returns a list with paths to all saved files.

## Details

Creates the following file structure in output_dir:

    output_dir/
      ├── plots/
      │   ├── <prefix>_cv_marginal_full.png
      │   ├── <prefix>_cv_marginal_zoom.png
      │   ├── <prefix>_cv_best_scores.png
      │   └── <prefix>_feature_weights.png
      ├── plot_objects/
      │   └── <prefix>_plots.rds (all ggplot objects)
      ├── code_elements/
      │   ├── <prefix>_coarse_cv.rds
      │   ├── <prefix>_fine_cv.rds
      │   ├── <prefix>_best_model.rds
      │   ├── <prefix>_coarse_grid.csv
      │   └── <prefix>_fine_grid.csv
      └── <prefix>_summary.txt

## Examples

``` r
if (FALSE) { # \dontrun{
# After running full CV workflow
saved_files <- ELEUTHIA_save_scca_results(
  output_dir = "results/analysis_2026_01_14",
  coarse_cv_results = cv_result_coarse,
  fine_cv_results = cv_result_fine,
  plots = list(
    cv_marginal_full = p_marginal_full,
    cv_marginal_zoom = p_marginal_zoom,
    cv_best_scores = p_best_scores,
    feature_weights = p_weights
  ),
  method = "ConvCCA",
  X = sim_data$datasets,
  pilot_results = pilot_result,
  metadata = list(k = 5, lambda = 10, nIter = 100, penalty = "LASSO")
)

} # }
```
