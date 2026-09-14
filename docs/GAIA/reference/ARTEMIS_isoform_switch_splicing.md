# Classify Alternative Splicing Events for Isoform Switches

Classifies each switching isoform's splicing relative to the gene's
hypothetical pre-mRNA (all exons of the gene concatenated), into
concrete event categories: intron retention (IR), exon skipping (ES),
multiple exon skipping (MES), mutually exclusive exons (MEE),
alternative 5'/3' splice sites (A5/A3), and alternative transcription
start/end sites (ATSS/ATTS). This is a distinct layer from
[`ARTEMIS_isoform_switch_consequences()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_isoform_switch_consequences.md):
consequences describe functional impact (domain loss, NMD, etc.), this
describes splicing mechanism – and needs no external tool, only the
GTF-derived exon structure already present after
[`ARTEMIS_isoform_switch()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_isoform_switch.md).

## Usage

``` r
ARTEMIS_isoform_switch_splicing(
  switch_result,
  only_switching_genes = TRUE,
  alpha = NULL,
  dIF_cutoff = NULL,
  verbose = TRUE
)
```

## Arguments

- switch_result:

  An `artemis_isoform_switch` object that has already been run through
  [`ARTEMIS_isoform_switch_consequences()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_isoform_switch_consequences.md).
  Required because
  [`IsoformSwitchAnalyzeR::extractSplicingSummary()`](https://rdrr.io/pkg/IsoformSwitchAnalyzeR/man/extractSplicingSummary.html)/
  `extractSplicingEnrichment()` are documented as expecting
  `analyzeSwitchConsequences()` to have already been run on the
  `switchAnalyzeRlist`.

- only_switching_genes:

  Logical. Restrict classification to genes with a significant isoform
  switch (per `alpha`/`dIF_cutoff`) rather than all genes in the
  `switchAnalyzeRlist`. Default: `TRUE`.

- alpha:

  Numeric or NULL. Significance threshold. Default: NULL, which uses
  `switch_result$params$alpha`.

- dIF_cutoff:

  Numeric or NULL. Minimum absolute dIF to consider an isoform
  switching. Default: NULL, which uses
  `switch_result$params$dIF_cutoff`.

- verbose:

  Logical. Print progress. Default: `TRUE`.

## Value

The input `artemis_isoform_switch` object with three additional fields:

- splicing_events:

  Data.frame from `sar$AlternativeSplicingAnalysis` – one row per
  isoform_id with the count of each splice event type and the genomic
  coordinates of the affected region(s)

- splicing_summary:

  Data.frame from `extractSplicingSummary()` – total number of each
  event type across all switches/genes

- splicing_enrichment:

  Data.frame from `extractSplicingEnrichment()` – for each event type,
  whether gain vs loss is enriched among the switches (FDR-corrected
  [`prop.test()`](https://rdrr.io/r/stats/prop.test.html))

`switch_list` is also updated with the added
`AlternativeSplicingAnalysis` slot.
