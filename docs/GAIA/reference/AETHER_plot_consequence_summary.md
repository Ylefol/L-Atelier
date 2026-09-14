# Functional Consequence Summary Bar Plot

Bar plot of the number of switches with each functional consequence type
(domain gain/loss, coding potential, signal peptide, subcellular
location, IDR, etc.), faceted by feature. Built on
[`IsoformSwitchAnalyzeR::extractConsequenceSummary()`](https://rdrr.io/pkg/IsoformSwitchAnalyzeR/man/extractConsequenceSummary.html),
called live on `switch_result$switch_list` – this is a cheap tabulation
over the consequence annotations already computed by
[`ARTEMIS_isoform_switch_consequences()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_isoform_switch_consequences.md),
so nothing needs to be pre-stored under a separate field the way
`$splicing_summary`/ `$splicing_enrichment` are.

## Usage

``` r
AETHER_plot_consequence_summary(
  switch_result,
  count_by = c("genes", "isoforms"),
  consequences_to_analyze = "all",
  title = NULL,
  fill_color = "#3A3A3A"
)
```

## Arguments

- switch_result:

  An `artemis_isoform_switch` object that has already been run through
  [`ARTEMIS_isoform_switch_consequences()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_isoform_switch_consequences.md).

- count_by:

  Character. `"genes"` (default) plots `nrGenesWithConsequences`;
  `"isoforms"` plots `nrIsoWithConsequences`.

- consequences_to_analyze:

  Character vector. Passed through to `extractConsequenceSummary()`'s
  `consequencesToAnalyze`. Default: `"all"` (every consequence type
  present in `switch_result$switch_list`).

- title:

  Character or NULL. Plot title. Default: NULL, which builds " vs " from
  `switch_result$params`.

- fill_color:

  Character. Bar fill color. Default: `"#3A3A3A"`.

## Value

A ggplot object.

## Details

Unlike
[`AETHER_plot_splicing_summary()`](https://ylefol.github.io/L-Atelier/GAIA/reference/AETHER_plot_splicing_summary.md),
consequence categories don't reduce to a clean two-way "more/less" split
– a feature like "Domains identified" has gain/loss/switch outcomes,
"Coding potential" has coding/non-coding, "Sub cell location" has
gain/loss/switch, etc. – so bars are a single color and faceted by
`featureCompared` (one panel per feature) with a free x scale, rather
than dodged by direction. Panels only appear for consequence types that
were actually annotated (i.e. the external tool files supplied to
[`ARTEMIS_isoform_switch_consequences()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_isoform_switch_consequences.md)).

**Direction**: every `switchConsequence` label (e.g. "Domain gain",
"Domain loss", "Transcript is Noncoding") describes the isoform
*upregulated in `switch_result$params$experiment`* relative to the
isoform used more in `switch_result$params$reference` –
`IsoformSwitchAnalyzeR`'s own convention (`?analyzeSwitchConsequences`:
"switchConsequence ... a short description of the features of the
upregulated isoform", where "upregulated" always means higher usage in
condition_2/experiment). So e.g. a tall "Domain gain" bar means many
genes have the experiment-dominant isoform gaining a domain the
reference-dominant isoform didn't have.

## Examples

``` r
if (FALSE) { # \dontrun{
p <- AETHER_plot_consequence_summary(switch_result)
p <- AETHER_plot_consequence_summary(switch_result, count_by = "isoforms")
} # }
```
