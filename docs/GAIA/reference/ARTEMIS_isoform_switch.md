# Artemis - Isoform Switch Analysis Functions

Differential isoform usage analysis using IsoformSwitchAnalyzeR.
Identifies genes where the relative usage of transcript isoforms changes
between conditions, independently of overall gene expression changes.
Differential Isoform Usage (Isoform Switch) Analysis

Tests for differential isoform usage between two conditions using
IsoformSwitchAnalyzeR. Identifies genes where the relative proportion of
transcripts changes between conditions, independently of total gene
expression changes.

Complements
[`ARTEMIS_differential_counts()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_differential_counts.md):
DE analysis detects changes in total gene output; isoform switch
analysis detects changes in which form of the gene is produced. A gene
can show switching without overall DE, and vice versa.

## Usage

``` r
ARTEMIS_isoform_switch(
  counts,
  targets,
  reference,
  experiment,
  group_col = "group",
  gtf_path,
  method = c("DEXSeq", "satuRn"),
  alpha = 0.05,
  dIF_cutoff = 0.1,
  min_gene_expression = 1,
  min_transcript_expression = 1,
  fix_stringtie_annotation = FALSE,
  detect_unwanted_effects = TRUE,
  verbose = TRUE
)
```

## Arguments

- counts:

  Numeric matrix of transcript counts (transcripts x samples). Rownames
  must be transcript IDs matching those in `gtf_path`. Raw integer
  counts (e.g. from IsoQuant or FLAMES) and estimated counts (e.g. from
  Salmon/kallisto) are both accepted.

- targets:

  Data.frame with sample metadata. Rownames must match colnames of
  `counts`. Must include a group column.

- reference:

  Character. The reference/baseline group label (denominator).

- experiment:

  Character. The experimental/comparison group label (numerator).
  Positive dIF means higher isoform usage in experiment.

- group_col:

  Character. Column in targets containing group labels. Default:
  `"group"`.

- gtf_path:

  Character. Path to GTF annotation file. Transcript IDs in the GTF must
  match rownames of `counts`.

- method:

  Character. Statistical testing method: `"DEXSeq"` (default, negative
  binomial) or `"satuRn"` (quasi-binomial; replaces DRIMSeq in
  IsoformSwitchAnalyzeR \>= 2.x).

- alpha:

  Numeric. FDR significance threshold. Default: `0.05`.

- dIF_cutoff:

  Numeric. Minimum absolute delta Isoform Fraction required. A switch
  must change isoform usage by at least this amount to be reported.
  Default: `0.1` (10\\ is not an established universal standard — verify
  for your context.

- min_gene_expression:

  Numeric. Minimum mean expression per gene. Default: `1`.

- min_transcript_expression:

  Numeric. Minimum mean expression per transcript. Default: `1`.

- fix_stringtie_annotation:

  Logical. Passed to `importRdata()`'s `fixStringTieAnnotationProblem`.
  This ISAS feature re-derives `gene_id` for de novo/StringTie-style
  assemblies (which assign synthetic gene IDs) by mapping transcripts
  back to a reference gene via genomic overlap — but it does so by
  substituting a gene *name*, not a stable gene ID, and gene names are
  not guaranteed unique across the genome (paralogs, pseudogenes,
  readthrough/antisense genes can share a base symbol). Confirmed
  empirically: enabling this on a clean, directly annotated reference
  GTF replaced the majority of `gene_id` values with gene symbols and
  produced spurious "same gene_id on different chromosomes" collisions,
  causing genes to be dropped entirely. Default: `FALSE` — only set
  `TRUE` if `gtf_path` is itself a StringTie/discovery-mode annotation
  with synthetic gene IDs that need rescuing.

- detect_unwanted_effects:

  Logical. Passed to `importRdata()`'s `detectUnwantedEffects`. When
  `TRUE` (IsoformSwitchAnalyzeR's own default, kept here), an SVA-style
  surrogate-variable search runs automatically and any detected
  covariates are added to the design matrix and corrected for before
  testing – see `@details` below for why this can silently erase power
  at low replicate counts. Default: `TRUE`.

- verbose:

  Logical. Print progress. Default: `TRUE`.

## Value

An S3 object of class `"artemis_isoform_switch"` containing:

- switch_list:

  Full `switchAnalyzeRlist` (native ISAS format; pass directly to
  IsoformSwitchAnalyzeR functions for advanced use)

- isoform_results:

  Data.frame of per-isoform results: gene_id, gene_name, isoform_id,
  condition_1, condition_2, dIF, isoform_switch_q_value,
  gene_switch_q_value, iso_significant. Sorted by q-value. `gene_name`
  comes straight from the GTF via `importRdata()` (NA if the GTF had no
  `gene_name` attribute for that gene) – the same source
  `consequence_summary$gene_name` uses, so the two tables share a common
  gene-naming column without needing a separate lookup/join.

- summary:

  Top-level switch count summary from `extractSwitchSummary()`

- n_switches:

  Number of genes with at least one significant isoform switch

- params:

  List of parameters used

## Details

One warning from IsoformSwitchAnalyzeR is expected and not a bug: "Using
row.names as isoform_id" simply confirms it used the counts matrix
rownames as transcript IDs (this function's design).

If you see "gene_ids or isoform_ids were not unique... removed N
gene_id" with `fix_stringtie_annotation = FALSE`, that reflects genuine
gene_ids spanning more than one chromosome in `gtf_path` — inspect your
annotation. If it appears with `fix_stringtie_annotation = TRUE`, it is
likely the gene-symbol substitution artifact described above rather than
a real annotation problem.

When `method = "satuRn"`, this function attaches SummarizedExperiment
internally (no-op if already attached). This works around a confirmed
upstream bug: `isoformSwitchTestSatuRn()` calls the unqualified
`rowData()` without importing it in IsoformSwitchAnalyzeR's NAMESPACE
(present in 2.2.0, the current Bioconductor release). Since this
function calls IsoformSwitchAnalyzeR via `::` rather than
[`library()`](https://rdrr.io/r/base/library.html), its Depends chain is
never attached, so `rowData()` would otherwise error with "could not
find function".

`method = "satuRn"` also disables
[`satuRn::testDTU()`](https://rdrr.io/pkg/satuRn/man/testDTU.html)'s
built-in p-value calibration diagplots (a Frequency histogram and a
Density plot, one pair per contrast). Those are satuRn's own
diagnostics, drawn as a side effect, not part of this package's plotting
layer (AETHER).

This function does *not* support user-supplied confounder/covariate
columns (e.g. known batch, sex, age). The design matrix built internally
only contains `sampleID` and `condition` columns. IsoformSwitchAnalyzeR
itself supports additional cofactor columns on the design matrix (taken
into account automatically by both `isoformSwitchTestDEXSeq()` and
`isoformSwitchTestSatuRn()`; see its vignette section "How to handle
confounding effects (including batches)"), but that path is not
currently exposed by this wrapper. If your data has *known* confounders,
correct for them prior to calling this function rather than relying on
it here.

Separately, `importRdata()` (called internally) runs its own automatic
*unwanted-effect* detection (an SVA-style surrogate-variable search,
controlled by `detect_unwanted_effects` below) independent of the
user-supplied covariates described above. When it finds candidate
surrogate variables, it adds them to the design matrix and
batch-corrects expression estimates before testing – this happens even
though no covariate columns were passed in. At low replicate counts
(e.g. n = 3/group), this detection is unreliable: with few samples,
group is often the dominant source of variance, so estimated surrogate
variables can end up correlated with `condition` itself, and regardless
of correlation, each one consumes a residual degree of freedom that a
small design can't spare. This can turn a real effect into zero
significant results even though a pre-test guesstimate (printed during
import) suggested many candidate genes. Consider
`detect_unwanted_effects = FALSE` for low-replicate designs.
