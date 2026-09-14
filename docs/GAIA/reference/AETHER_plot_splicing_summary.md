# Splicing Event Summary Bar Plot

Bar plot of the total number of alternative splicing events per type
(IR, A5, A3, ATSS, ATTS, ES, MES, MEE), split by whether the event
occurs in the isoform used *more* or *less* in the experiment group,
from `switch_result$splicing_summary`
([`IsoformSwitchAnalyzeR::extractSplicingSummary()`](https://rdrr.io/pkg/IsoformSwitchAnalyzeR/man/extractSplicingSummary.html)
output).

## Usage

``` r
AETHER_plot_splicing_summary(
  switch_result,
  count_by = c("genes", "isoforms"),
  title = NULL,
  colors = NULL
)
```

## Arguments

- switch_result:

  An `artemis_isoform_switch` object with splicing classification
  already run
  ([`ARTEMIS_isoform_switch_splicing()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_isoform_switch_splicing.md)).

- count_by:

  Character. `"genes"` (default) plots `nrGenesWithConsequences`;
  `"isoforms"` plots `nrIsoWithConsequences`.

- title:

  Character or NULL. Plot title. Default: NULL, which builds " vs " from
  `switch_result$params`.

- colors:

  Named character vector with colors for `"more"` and `"less"`. Default
  uses red/blue.

## Value

A ggplot object.

## Details

This is the total count of events per type – it does not indicate
whether the events belong to significant switches specifically (see
[`AETHER_plot_splicing_enrichment()`](https://ylefol.github.io/L-Atelier/GAIA/reference/AETHER_plot_splicing_enrichment.md)
for the significance-tested gain/loss comparison). Bars are grouped by
`AStype` and dodged by direction (isoform used more vs. less in the
experiment group).

**Direction**: "isoform used more/less" follows
`IsoformSwitchAnalyzeR`'s `isoformUpregulated`/ `isoformDownregulated`
convention from `analyzeSwitchConsequences()` – "more" is the isoform
with higher usage in `switch_result$params$experiment` (positive dIF,
relative to `switch_result$params$reference`), not an arbitrary pairwise
label. So e.g. a tall "IR in isoform used more" bar means many genes
gain intron retention specifically in the experiment-dominant isoform.

## Examples

``` r
if (FALSE) { # \dontrun{
p <- AETHER_plot_splicing_summary(switch_result)
p <- AETHER_plot_splicing_summary(switch_result, count_by = "isoforms")
} # }
```
