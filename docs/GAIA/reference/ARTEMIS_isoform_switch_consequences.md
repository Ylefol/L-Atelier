# Analyze Functional Consequences of Isoform Switches

Adds biological consequence annotations to isoform switch results.
Integrates outputs from external tools (CPAT, SignalP, Pfam/HMMER,
DeepTMHMM, DeepLoc 2.1, IUPred2A) and identifies functional differences
between switching isoform pairs: changes in coding potential, NMD
sensitivity, protein domains, signal peptides, membrane topology,
subcellular localization, intrinsically disordered regions, and intron
retention.

## Usage

``` r
ARTEMIS_isoform_switch_consequences(
  switch_result,
  cpat_file = NULL,
  signalp_file = NULL,
  pfam_file = NULL,
  deeptmhmm_file = NULL,
  deeploc2_file = NULL,
  iupred2a_file = NULL,
  cpat_cutoff = 0.725,
  verbose = TRUE
)
```

## Arguments

- switch_result:

  An `artemis_isoform_switch` object from
  [`ARTEMIS_isoform_switch()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_isoform_switch.md).

- cpat_file:

  Character or NULL. Path to CPAT output file for coding potential
  annotation. Default: NULL (skipped).

- signalp_file:

  Character or NULL. Path to SignalP output file for signal peptide
  annotation. Default: NULL (skipped).

- pfam_file:

  Character or NULL. Path to Pfam/HMMER output file for protein domain
  annotation. Default: NULL (skipped).

- deeptmhmm_file:

  Character or NULL. Path to a DeepTMHMM `TMRs.gff3` result file (e.g.
  from
  [`APOLLO_run_deeptmhmm()`](https://ylefol.github.io/L-Atelier/GAIA/reference/APOLLO_run_deeptmhmm.md))
  for cell-membrane topology annotation. Default: NULL (skipped).

- deeploc2_file:

  Character or NULL. Path to a DeepLoc 2.1 result file (e.g. from
  [`APOLLO_run_deeploc2()`](https://ylefol.github.io/L-Atelier/GAIA/reference/APOLLO_run_deeploc2.md),
  already reformatted to IsoformSwitchAnalyzeR's expected 13-column
  layout) for subcellular localization annotation. Default: NULL
  (skipped).

- iupred2a_file:

  Character or NULL. Path to a combined IUPred2A result file (e.g. from
  [`APOLLO_run_iupred2a()`](https://ylefol.github.io/L-Atelier/GAIA/reference/APOLLO_run_iupred2a.md))
  for intrinsically disordered region (IDR) and ANCHOR2 binding-site
  annotation. Default: NULL (skipped).

- cpat_cutoff:

  Numeric. Coding probability cutoff for CPAT. Human default is 0.725,
  mouse is 0.44. No universal standard — verify for your organism. Only
  used if `cpat_file` is provided. Default: `0.725`.

- verbose:

  Logical. Print progress. Default: `TRUE`.

## Value

The input `artemis_isoform_switch` object with two additional fields:

- consequence_summary:

  Data.frame from `sar$switchConsequence` — one row per isoform-switch
  pair with annotated consequence types

- switch_list:

  Updated switchAnalyzeRlist with consequence annotations

## Details

At least one of `cpat_file`, `signalp_file`, `pfam_file`,
`deeptmhmm_file`, `deeploc2_file`, or `iupred2a_file` must be provided.
Structural consequences (intron retention, ORF sequence similarity, NMD
sensitivity) are always included and require no external tools, only ORF
annotations from the initial
[`ARTEMIS_isoform_switch()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_isoform_switch.md)
call; intron retention is classified internally via
[`IsoformSwitchAnalyzeR::analyzeIntronRetention()`](https://rdrr.io/pkg/IsoformSwitchAnalyzeR/man/analyzeAlternativSplicing.html)
before consequence analysis runs.

When `deeptmhmm_file` is supplied, this adds the `"isoform_topology"`
consequence type (differences in predicted intracellular/transmembrane/
extracellular topology between switching isoforms) – it does not add the
related `extracellular_region_count`/`intracellular_region_count`/
`extracellular_region_length`/`intracellular_region_length` types
`analyzeSwitchConsequences()` also supports; call
[`IsoformSwitchAnalyzeR::analyzeSwitchConsequences()`](https://rdrr.io/pkg/IsoformSwitchAnalyzeR/man/analyzeSwitchConsequences.html)
directly on `switch_result$switch_list` if those are needed.

When `deeploc2_file` is supplied, this adds the `"sub_cell_location"`
consequence type (differences in predicted subcellular localization
between switching isoforms). Unlike `signalp_file`, there is no "has
data" tracking needed here – `analyzeDeepLoc2()` always populates
`sub_cell_location` (with `NA` for isoforms with no called
localization), it never skips setting the slot the way
`analyzeSignalP()` can.

When `iupred2a_file` is supplied, this adds both `"IDR_identified"` and
`"IDR_type"` consequence types – `IDR_type` (plain IDR vs.
`IDR_w_binding_region`) comes for free once ANCHOR2 binding-site data is
present, which
[`APOLLO_run_iupred2a()`](https://ylefol.github.io/L-Atelier/GAIA/reference/APOLLO_run_iupred2a.md)
always includes.

External tools must be run independently before calling this function.
See the IsoformSwitchAnalyzeR vignette for expected file formats and
instructions for running each tool.
