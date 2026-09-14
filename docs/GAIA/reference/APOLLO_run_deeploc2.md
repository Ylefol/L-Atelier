# Run DeepLoc 2.1 Subcellular Localization Prediction

Runs DeepLoc 2.1 (subcellular localization + membrane association
prediction) against the amino acid FASTA from
[`APOLLO_extract_isoform_sequences()`](https://ylefol.github.io/L-Atelier/GAIA/reference/APOLLO_extract_isoform_sequences.md),
inside an isolated basilisk Python environment (`.gaia_deeploc2_env`).

## Usage

``` r
APOLLO_run_deeploc2(
  switch_result,
  output_dir,
  aa_fasta = NULL,
  model = c("Fast", "Accurate"),
  verbose = TRUE
)
```

## Arguments

- switch_result:

  An `artemis_isoform_switch` object, normally one that has already been
  through
  [`APOLLO_extract_isoform_sequences()`](https://ylefol.github.io/L-Atelier/GAIA/reference/APOLLO_extract_isoform_sequences.md).

- output_dir:

  Character. Directory for DeepLoc 2.1's raw output and the reformatted
  result file (created if absent).

- aa_fasta:

  Character or `NULL`. Amino acid FASTA path. Default `NULL` uses
  `switch_result$sequence_files$aa_fasta`.

- model:

  Character. `"Fast"` (ESM1b, default – higher throughput) or
  `"Accurate"` (ProtT5, ~32GB memory, downloaded on first use). Passed
  to DeepLoc 2.1's `-m` option.

- verbose:

  Logical. Print progress and DeepLoc 2.1's own console output. Default:
  `TRUE`.

## Value

Invisibly, the path to a reformatted result file compatible with
`IsoformSwitchAnalyzeR::analyzeDeepLoc2(pathToDeepLoc2resultFile = ...)`.

## Details

DeepLoc 2.1's own CLI output has 18 columns: `Protein_ID`,
`Localizations`, `Signals`, `"Membrane types"`, the 10 localization
probability columns, and 4 membrane-type probability columns
(`Peripheral`/`Transmembrane`/`Lipid anchor`/`Soluble` – a multi-label
membrane-type feature DeepLoc 2.1 added on top of DeepLoc 2.0).
[`IsoformSwitchAnalyzeR::analyzeDeepLoc2()`](https://rdrr.io/pkg/IsoformSwitchAnalyzeR/man/analyzeDeepLoc2.html)
(confirmed by inspecting its source directly) does an exact `ncol == 13`
/ exact-column-name check against the older DeepLoc 2.0-era 13-column
format, with no "Membrane types" column and no membrane-type
probabilities – this predates DeepLoc 2.1's 2024 publication, and the
currently installed IsoformSwitchAnalyzeR (2.2.0, current Bioconductor
release as of 2026-07) has not been updated for the newer format. Worth
checking for a newer IsoformSwitchAnalyzeR release before assuming this
reformatting step will always be needed. This function reformats DeepLoc
2.1's raw `results_*.csv` by dropping `"Membrane types"` and the 4
membrane-type probability columns, keeping the 10 original localization
columns unchanged (same names/ values), before returning the path.

DeepLoc 2.1 is DTU Health Tech-distributed (not on public PyPI/conda),
so `.gaia_deeploc2_env` only provisions a base Python interpreter –
DeepLoc 2.1 itself (and its own torch/fair-esm/transformers/
pytorch_lightning dependencies, pulled in automatically via its own
`setup.py`) must be installed once, manually, into that environment
(`pip install .` from the downloaded package directory, using this
environment's own pip; see
`basilisk::obtainEnvironmentPath(.gaia_deeploc2_env)`). This function
checks for the installed `deeploc2` executable and stops with setup
guidance if it's missing, rather than attempting to install it
automatically (same pattern as
[`APOLLO_run_signalp()`](https://ylefol.github.io/L-Atelier/GAIA/reference/APOLLO_run_signalp.md)).
