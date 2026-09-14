# Evaluate Sparse CCA Model

Computes evaluation metrics for sCCA results

## Usage

``` r
MINERVA_evaluate_scca(X_list, W_list, metric = "correlation")
```

## Arguments

- X_list:

  List of datasets (test data)

- W_list:

  List of canonical vectors

- metric:

  Metric to compute: "correlation" (default) or "both" (both computes
  correlation and sparsity metrics)

## Value

Named list of metrics
