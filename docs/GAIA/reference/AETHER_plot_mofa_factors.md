# Plot MOFA factor scores for two factors

Creates a scatter plot of sample scores on two MOFA factors, optionally
colored by an external metadata variable (e.g. clinical group,
timepoint).

## Usage

``` r
AETHER_plot_mofa_factors(
  mofa_result,
  factors = c(1, 2),
  color_by = NULL,
  color_label = "Group",
  point_size = 2,
  title = NULL
)
```

## Arguments

- mofa_result:

  A hephaestus_mofa object from HEPHAESTUS_run_mofa() or
  HEPHAESTUS_load_mofa_model().

- factors:

  Length-2 numeric vector, which factor indices to plot on x/y. Default:
  c(1, 2).

- color_by:

  Optional named vector/factor, sample ID -\> value, used to color
  points (e.g. a clinical metadata column keyed by sample ID). Not part
  of the hephaestus_mofa object itself – supply it separately.

- color_label:

  Legend title when color_by is supplied. Default: "Group".

- point_size:

  Point size. Default: 2.

- title:

  Plot title. If NULL, auto-generated from the factor names.

## Value

A ggplot object.

## Examples

``` r
if (FALSE) { # \dontrun{
result <- HEPHAESTUS_run_mofa(X = list(rna = view_a, atac = view_b))
AETHER_plot_mofa_factors(result, factors = c(1, 2))

# Colored by clinical group (named vector keyed by sample ID)
AETHER_plot_mofa_factors(result, factors = c(1, 2),
                          color_by = clinical_group, color_label = "Group")

} # }
```
