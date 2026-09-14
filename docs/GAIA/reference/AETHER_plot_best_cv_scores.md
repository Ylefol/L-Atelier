# Plot Best CV Score Comparison Between Coarse and Fine Grid

Creates a simple bar chart comparing the best CV scores from coarse and
fine grid search. Designed for clarity in presentations and publications
for non-technical audiences.

## Usage

``` r
AETHER_plot_best_cv_scores(
  coarse_cv_results,
  fine_cv_results,
  metric_name = "CV Score"
)
```

## Arguments

- coarse_cv_results:

  Data frame from MINERVA_cv_scca with coarse grid. Must contain column:
  mean_score

- fine_cv_results:

  Data frame from MINERVA_cv_scca with fine grid. Must contain column:
  mean_score

- metric_name:

  Character string for y-axis label (default = "CV Score")

## Value

ggplot2 object showing bar chart of best scores from each grid.

## Details

Creates a simple two-bar comparison showing:

- Best CV score from coarse grid search

- Best CV score from fine grid search (after refinement)

This visualization clearly demonstrates whether the fine grid refinement
improved performance over the initial coarse search, without the
complexity of marginal plots.

The function also prints:

- Best score from each grid

- Absolute improvement (fine - coarse)

- Percentage improvement

## Examples

``` r
if (FALSE) { # \dontrun{
# After running coarse and fine CV
p <- AETHER_plot_best_cv_scores(
  coarse_cv_results = coarse_cv$cv_results,
  fine_cv_results = fine_cv$cv_results,
  metric_name = "Correlation"
)
print(p)

} # }
```
