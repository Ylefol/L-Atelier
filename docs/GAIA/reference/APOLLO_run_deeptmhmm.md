# Run DeepTMHMM Topology Prediction

Runs DeepTMHMM (transmembrane helix / signal peptide / cell-membrane
topology prediction) against the amino acid FASTA from
[`APOLLO_extract_isoform_sequences()`](https://ylefol.github.io/L-Atelier/GAIA/reference/APOLLO_extract_isoform_sequences.md),
inside an isolated basilisk Python environment (`.gaia_deeptmhmm_env`).

## Usage

``` r
APOLLO_run_deeptmhmm(
  switch_result,
  deeptmhmm_dir,
  output_dir,
  aa_fasta = NULL,
  verbose = TRUE
)
```

## Arguments

- switch_result:

  An `artemis_isoform_switch` object, normally one that has already been
  through
  [`APOLLO_extract_isoform_sequences()`](https://ylefol.github.io/L-Atelier/GAIA/reference/APOLLO_extract_isoform_sequences.md).

- deeptmhmm_dir:

  Character. Path to the locally installed, academically licensed
  DeepTMHMM package directory – a flat checkout containing `predict.py`,
  `utils.py`, the `experiments/` folder, and the five
  `deeptmhmm_cv_*.model` weight files.

- output_dir:

  Character. Directory `predict.py` writes its output to. Must **not**
  already exist – `predict.py` itself checks for this and exits with an
  error if `output_dir` is already present (unlike
  [`APOLLO_run_cpat()`](https://ylefol.github.io/L-Atelier/GAIA/reference/APOLLO_run_cpat.md)/[`APOLLO_run_signalp()`](https://ylefol.github.io/L-Atelier/GAIA/reference/APOLLO_run_signalp.md)/
  [`APOLLO_run_pfam()`](https://ylefol.github.io/L-Atelier/GAIA/reference/APOLLO_run_pfam.md),
  this function does not pre-create it).

- aa_fasta:

  Character or `NULL`. Amino acid FASTA path. Default `NULL` uses
  `switch_result$sequence_files$aa_fasta`.

- verbose:

  Logical. Print progress and DeepTMHMM's own console output. Default:
  `TRUE`.

## Value

Invisibly, the path to the `TMRs.gff3` result file, directly compatible
with
`IsoformSwitchAnalyzeR::analyzeDeepTMHMM(pathToDeepTMHMMresultFile = ...)`
/ `ARTEMIS_isoform_switch_consequences(deeptmhmm_file = ...)`.

## Details

DeepTMHMM is academically licensed (DTU Health Tech / BioLib) and not on
public PyPI/bioconda channels, so `.gaia_deeptmhmm_env` only provisions
a base Python 3.8 interpreter – `predict.py`'s own dependencies
(PyTorch, per its bundled `requirements.txt`) must be installed once,
manually, into that environment following the license holder's own
`README.txt` before calling this function. This function checks for
`predict.py` in `deeptmhmm_dir` and stops with setup guidance if it's
missing, rather than attempting to install it automatically (same
pattern as
[`APOLLO_run_signalp()`](https://ylefol.github.io/L-Atelier/GAIA/reference/APOLLO_run_signalp.md)).

Unlike CPAT/SignalP/Pfam, `predict.py` is invoked from within its own
package directory rather than as an installed executable on `PATH`: it
loads its five `deeptmhmm_cv_*.model` weight files and its own
`utils`/`experiments.tmhmm3.tm_util` modules via bare relative
paths/imports (confirmed by reading `predict.py` directly), so this
function temporarily changes the working directory to `deeptmhmm_dir`
for the duration of the call (restored via
[`on.exit()`](https://rdrr.io/r/base/on.exit.html), including on error).

`predict.py`'s own `requirements.txt` pins `torch==1.5.0+cu92` (a CUDA
9.2-era wheel) – confirm this actually resolves on your hardware/driver
stack before assuming it as fixed; a newer CPU-only or different
CUDA-version torch build may be needed instead. GPU vs CPU execution is
auto-detected by `predict.py` itself (`torch.cuda.is_available()`) and
is not controlled by this wrapper.
