# Plot CV Performance Comparison Between Coarse and Fine Grid Search

Creates marginal CV performance plots comparing coarse and fine grid
search results. Shows how CV score changes with each tau parameter,
demonstrating the refinement gained from fine grid search.

## Usage

``` r
AETHER_plot_cv_comparison(
  coarse_cv_results,
  fine_cv_results,
  n_datasets,
  param_names = NULL,
  metric_name = "CV Score",
  x_axis_type = "full"
)
```

## Arguments

- coarse_cv_results:

  Data frame from MINERVA_cv_scca with coarse grid. Must contain
  columns: param1, param2, ..., paramN, mean_score

- fine_cv_results:

  Data frame from MINERVA_cv_scca with fine grid. Must contain columns:
  param1, param2, ..., paramN, mean_score

- n_datasets:

  Integer. Number of datasets (tau parameters).

- param_names:

  Optional character vector of parameter names for axis labels. If NULL,
  uses "tau1", "tau2", etc. (default = NULL)

- metric_name:

  Character string for y-axis label (default = "CV Score")

- x_axis_type:

  Character string controlling x-axis range. Options: "full" = fixed 0-1
  range showing complete parameter space (default), "auto" = zoom to
  data range (fine grid region)

## Value

ggplot2 object with faceted panels showing marginal CV performance
curves for coarse and fine grids.

## Details

Creates a multi-panel plot with one panel per tau parameter. Each panel
shows:

- **Coarse grid curve**: CV performance across coarse tau values
  (averaged over other tau parameters)

- **Fine grid curve**: CV performance in refined range (averaged over
  other tau parameters)

- **Markers**: Best value from fine grid highlighted

Marginal plots are computed by:

1.  Grouping by each tau parameter

2.  Averaging CV score across all other tau values

3.  Plotting the marginal relationship

This visualization demonstrates:

- Initial broad exploration (coarse grid)

- Focused refinement (fine grid in promising region)

- Performance improvement from fine-tuning

## Examples

``` r
if (FALSE) { # \dontrun{
# After running coarse and fine CV

# Full view showing complete 0-1 parameter space
p_full <- AETHER_plot_cv_comparison(
  coarse_cv_results = coarse_cv$cv_results,
  fine_cv_results = fine_cv$cv_results,
  n_datasets = 3,
  param_names = c("RNA-seq tau", "ATAC-seq tau", "CUT&TAG tau"),
  x_axis_type = "full"
)
print(p_full)

# Zoomed view focusing on fine grid region
p_zoom <- AETHER_plot_cv_comparison(
  coarse_cv_results = coarse_cv$cv_results,
  fine_cv_results = fine_cv$cv_results,
  n_datasets = 3,
  param_names = c("RNA-seq tau", "ATAC-seq tau", "CUT&TAG tau"),
  x_axis_type = "auto"
)
print(p_zoom)

} # }
```
