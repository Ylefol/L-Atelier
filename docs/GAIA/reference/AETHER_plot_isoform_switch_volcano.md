# Volcano Plot for Isoform Switch Results

Creates a four-category volcano plot (-log10 q-value vs delta Isoform
Fraction) from an `artemis_isoform_switch` object, with embedded counts
in the legend and optional gene/isoform labeling. Categories: increased
usage, decreased usage, low-dIF (significant but below the dIF cutoff),
and non-significant.

## Usage

``` r
AETHER_plot_isoform_switch_volcano(
  switch_result,
  label_col = "gene_id",
  genes_of_interest = NULL,
  show_non_sig_interest = TRUE,
  label_top_n = 0,
  q_col = "isoform_switch_q_value",
  dIF_cutoff = NULL,
  alpha = NULL,
  title = NULL,
  colors = NULL,
  point_size = 0.8,
  point_alpha = 0.7
)
```

## Arguments

- switch_result:

  An `artemis_isoform_switch` object from
  [`ARTEMIS_isoform_switch()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_isoform_switch.md).

- label_col:

  Character. Column in `isoform_results` to use for point labels.
  Default: `"gene_id"`.

- genes_of_interest:

  Character vector. Values in `label_col` to always label. Default:
  NULL.

- show_non_sig_interest:

  Logical. If FALSE, genes of interest that do not meet both thresholds
  are not labeled. Default: TRUE.

- label_top_n:

  Integer. Label the top N isoforms by `q_col` regardless of direction.
  Default: 0 (disabled).

- q_col:

  Character. Q-value column to use for significance filtering and the
  y-axis: `"isoform_switch_q_value"` (default, isoform-level) or
  `"gene_switch_q_value"` (gene-level).

- dIF_cutoff:

  Numeric or NULL. Minimum absolute dIF to call a switch
  "increased"/"decreased" usage rather than "low-dIF". Default: NULL,
  which uses `switch_result$params$dIF_cutoff`.

- alpha:

  Numeric or NULL. Significance threshold on `q_col`. Default: NULL,
  which uses `switch_result$params$alpha`.

- title:

  Character or NULL. Plot title. Default: NULL, which builds " vs " from
  `switch_result$params`.

- colors:

  Named character vector with colors for `"up"`, `"down"`, `"low_dIF"`,
  `"non_sig"`. Default uses red/blue/green/gray.

- point_size:

  Numeric. Point size. Default: 0.8.

- point_alpha:

  Numeric. Point transparency. Default: 0.7.

## Value

A ggplot object.

## Details

The four categories are:

- **increased usage**: `dIF > dIF_cutoff` AND `q_col < alpha`

- **decreased usage**: `dIF < -dIF_cutoff` AND `q_col < alpha`

- **low-dIF**: `|dIF| <= dIF_cutoff` AND `q_col < alpha`

- **non-significant**: `q_col >= alpha` or NA

A horizontal dashed line marks `alpha` on the `-log10(q_col)` scale.
Vertical dashed lines mark `±dIF_cutoff`.

## Examples

``` r
if (FALSE) { # \dontrun{
p <- AETHER_plot_isoform_switch_volcano(switch_result)

p <- AETHER_plot_isoform_switch_volcano(switch_result, label_top_n = 10,
                                        genes_of_interest = c("ENSG00000141510"))

p <- AETHER_plot_isoform_switch_volcano(switch_result,
                                        q_col = "gene_switch_q_value")
} # }
```
