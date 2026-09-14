# Visualize Batch Effects (Before/After Correction)

Creates PCA plots to visualize batch effects before and optionally after
batch correction.

## Usage

``` r
POSEIDON_plot_batch_pca(
  quant_result,
  batch_col = "batch",
  group_col = "group",
  log_transform = TRUE,
  title = "PCA - Batch Effect",
  show_labels = FALSE,
  label_size = 3
)
```

## Arguments

- quant_result:

  A list containing counts and targets.

- batch_col:

  Character string. Column name for batch (default = "batch").

- group_col:

  Character string. Column name for group (default = "group").

- log_transform:

  Logical. Log-transform counts for PCA (default = TRUE).

- title:

  Character string. Plot title (default = "PCA - Batch Effect").

- show_labels:

  Logical. Label each point with its sample ID. Uses ggrepel for
  non-overlapping placement if installed (default = FALSE).

- label_size:

  Numeric. Text size for point labels (default = 3).

## Value

A ggplot object showing PCA colored by batch and shaped by group.
