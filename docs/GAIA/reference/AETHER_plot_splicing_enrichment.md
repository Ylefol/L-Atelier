# Splicing Event Enrichment Plot

Point-range plot of the proportion of gain-vs-loss events per splicing
type (`propUp`, with confidence interval), from
`switch_result$splicing_enrichment`
([`IsoformSwitchAnalyzeR::extractSplicingEnrichment()`](https://rdrr.io/pkg/IsoformSwitchAnalyzeR/man/extractSplicingEnrichment.html)
output). Unlike
[`AETHER_plot_splicing_summary()`](https://ylefol.github.io/L-Atelier/GAIA/reference/AETHER_plot_splicing_summary.md),
this compares the two isoforms within each switch directly and is
significance-tested per type
([`prop.test()`](https://rdrr.io/r/stats/prop.test.html),
FDR-corrected).

## Usage

``` r
AETHER_plot_splicing_enrichment(
  switch_result,
  min_events = 10,
  title = NULL,
  colors = NULL
)
```

## Arguments

- switch_result:

  An `artemis_isoform_switch` object with splicing classification
  already run
  ([`ARTEMIS_isoform_switch_splicing()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_isoform_switch_splicing.md)).

- min_events:

  Integer. Minimum total events (`nUp + nDown`) required for a splicing
  type to be plotted. Default: 10, matching
  `extractSplicingEnrichment()`'s own `minEventsForPlotting` default.

- title:

  Character or NULL. Plot title. Default: NULL, which builds " vs " from
  `switch_result$params`.

- colors:

  Named character vector with colors for `"TRUE"` and `"FALSE"`
  significance. Default uses red/gray.

## Value

A ggplot object.

## Details

`propUp` is the fraction of gain/loss events (per type) that are a
*gain* – e.g. for `IR`, the fraction of IR gain-or-loss events that are
a gain of intron retention in the isoform used more (higher usage,
positive dIF) in `switch_result$params$experiment` relative to
`switch_result$params$reference` – `IsoformSwitchAnalyzeR`'s
`isoformUpregulated` convention, same as
[`AETHER_plot_splicing_summary()`](https://ylefol.github.io/L-Atelier/GAIA/reference/AETHER_plot_splicing_summary.md).
A dashed reference line at 0.5 marks "no directional preference"; points
whose confidence interval excludes 0.5 are the ones
[`prop.test()`](https://rdrr.io/r/stats/prop.test.html) calls
significant (FDR-corrected, `Significant` column). Point size is scaled
by the total number of events (`nUp + nDown`) so sparsely-supported
types are visually de-emphasized.

## Examples

``` r
if (FALSE) { # \dontrun{
p <- AETHER_plot_splicing_enrichment(switch_result)
p <- AETHER_plot_splicing_enrichment(switch_result, min_events = 5)
} # }
```
