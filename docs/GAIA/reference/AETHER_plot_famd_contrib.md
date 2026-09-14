# Plot FAMD Contribution Heatmap

Heatmap showing variable contributions to each FAMD dimension. Cleaner
alternative to bar plots for understanding which variables drive which
dimensions.

## Usage

``` r
AETHER_plot_famd_contrib(
  famd_result,
  n_dims = 5,
  var_type = "all",
  top_n = NULL,
  title = "FAMD - Variable Contributions"
)
```

## Arguments

- famd_result:

  An artemis_famd object from ARTEMIS_famd()

- n_dims:

  Integer. Number of dimensions to show. Default = 5.

- var_type:

  Character. Which variables to show: "all" (default), "quanti", or
  "quali".

- top_n:

  Integer. Only show top N contributing variables per dimension. If NULL
  (default), show all variables.

- title:

  Character. Plot title. Default = "FAMD - Variable Contributions".

## Value

A ggplot object
