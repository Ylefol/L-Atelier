# Per-Gene Isoform Switch Plot

Thin convenience wrapper around
[`IsoformSwitchAnalyzeR::switchPlot()`](https://rdrr.io/pkg/IsoformSwitchAnalyzeR/man/switchPlot.html)
that fills in `condition1`/`condition2`/`dIFcutoff`/`alphas` from an
`artemis_isoform_switch` object's `$params`, so only a gene (or isoform)
needs to be supplied.

## Usage

``` r
AETHER_plot_isoform_switch_gene(
  switch_result,
  gene_id = NULL,
  isoform_id = NULL,
  IF_cutoff = 0.05,
  dIF_cutoff = NULL,
  alphas = NULL,
  ...
)
```

## Arguments

- switch_result:

  An `artemis_isoform_switch` object from
  [`ARTEMIS_isoform_switch()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_isoform_switch.md).

- gene_id:

  Character. A single `gene_id` to plot, matching the `gene_id` column
  of `switch_result$isoform_results` (not a gene name/symbol). Mutually
  exclusive with `isoform_id`.

- isoform_id:

  Character vector. One or more isoform IDs (from the same gene) to
  plot, matching `switch_result$isoform_results$isoform_id`. Alternative
  to `gene_id`.

- IF_cutoff:

  Numeric. Minimum isoform-fraction contribution to gene expression (in
  at least one condition) required for an isoform to be plotted. This is
  `IsoformSwitchAnalyzeR`'s own display filter, unrelated to
  `switch_result$params$dIF_cutoff`. Default: `0.05`.

- dIF_cutoff:

  Numeric or NULL. dIF cutoff used to annotate increased/decreased usage
  on the transcript plot. Default: NULL, which uses
  `switch_result$params$dIF_cutoff`.

- alphas:

  Numeric vector of length two, or NULL. Q-value thresholds for "*" and
  "*\*\*" significance annotation, respectively. Default: NULL, which
  uses `c(switch_result$params$alpha, 0.001)`.

- ...:

  Additional arguments passed to
  [`IsoformSwitchAnalyzeR::switchPlot()`](https://rdrr.io/pkg/IsoformSwitchAnalyzeR/man/switchPlot.html)
  (e.g. `rescaleTranscripts`, `logYaxis`, `localTheme`,
  `additionalArguments`).

## Value

Invisibly returns the result of
[`IsoformSwitchAnalyzeR::switchPlot()`](https://rdrr.io/pkg/IsoformSwitchAnalyzeR/man/switchPlot.html).
As with the underlying function, the composite plot is drawn directly to
the current graphics device as a side effect — it is not a single
combinable ggplot object.

## Details

`condition1`/`condition2` are taken from
`switch_result$params$reference`/`$experiment` respectively — this
matches how
[`ARTEMIS_isoform_switch()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_isoform_switch.md)
builds `comparisonsToMake` internally (`condition_1 = reference`,
`condition_2 = experiment`), so orientation is guaranteed consistent
with the rest of the `artemis_isoform_switch` object (e.g. the sign of
`dIF` in
[`AETHER_plot_isoform_switch_volcano()`](https://ylefol.github.io/L-Atelier/GAIA/reference/AETHER_plot_isoform_switch_volcano.md)).

`gene_id`/`isoform_id` are validated against
`switch_result$isoform_results` before calling `switchPlot()`, so a typo
or a gene removed during `preFilter()` fails with a clear message rather
than an opaque error from IsoformSwitchAnalyzeR.

## Examples

``` r
if (FALSE) { # \dontrun{
AETHER_plot_isoform_switch_gene(switch_result, gene_id = "ENSG00000141510")

AETHER_plot_isoform_switch_gene(switch_result,
                                isoform_id = c("ENST00000269305", "ENST00000504290"))
} # }
```
