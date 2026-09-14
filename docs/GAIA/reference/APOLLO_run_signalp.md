# Run SignalP 6 Signal Peptide Prediction

Runs SignalP 6 against the amino acid FASTA from
[`APOLLO_extract_isoform_sequences()`](https://ylefol.github.io/L-Atelier/GAIA/reference/APOLLO_extract_isoform_sequences.md),
inside an isolated basilisk Python environment (`.gaia_signalp_env`).

## Usage

``` r
APOLLO_run_signalp(
  switch_result,
  output_dir,
  aa_fasta = NULL,
  organism = c("euk", "other"),
  mode = c("fast", "slow", "slow-sequential"),
  verbose = TRUE
)
```

## Arguments

- switch_result:

  An `artemis_isoform_switch` object, normally one that has already been
  through
  [`APOLLO_extract_isoform_sequences()`](https://ylefol.github.io/L-Atelier/GAIA/reference/APOLLO_extract_isoform_sequences.md).

- output_dir:

  Character. Directory for SignalP's output (created if absent). SignalP
  writes `prediction_results.txt` here.

- aa_fasta:

  Character or `NULL`. Amino acid FASTA path. Default `NULL` uses
  `switch_result$sequence_files$aa_fasta`.

- organism:

  Character. `"euk"` (eukaryote, default) or `"other"`, passed to
  SignalP's `--organism`.

- mode:

  Character. SignalP 6 prediction mode: `"fast"` (default), `"slow"`, or
  `"slow-sequential"` (higher accuracy, much slower).

- verbose:

  Logical. Print progress and SignalP's own console output. Default:
  `TRUE`.

## Value

Invisibly, the path to `prediction_results.txt`, directly compatible
with `ARTEMIS_isoform_switch_consequences(signalp_file = ...)` –
confirmed by inspecting `analyzeSignalP()`'s source, which natively
parses SignalP 6's output format (detected via a `"SignalP-6"` header
line), no reformatting needed.

## Details

SignalP 6 is license-gated (DTU Health Tech) and not on public PyPI, so
`.gaia_signalp_env` only provisions a base Python interpreter – SignalP
6 itself must be installed once, manually, into that environment
(download the licensed wheel, then install with the environment's own
pip; see `basilisk::obtainEnvironmentPath(.gaia_signalp_env)` to locate
it). Follow DTU's installation instructions in full, including the
separate model-weights copy step. This function checks for the installed
`signalp6` executable and stops with setup guidance if it's missing,
rather than attempting to install it automatically.
