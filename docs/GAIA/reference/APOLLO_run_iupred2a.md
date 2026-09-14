# Run IUPred2A Intrinsically Disordered Region Prediction

Runs IUPred2A (intrinsically disordered region + ANCHOR2 binding-site
prediction) against the amino acid FASTA from
[`APOLLO_extract_isoform_sequences()`](https://ylefol.github.io/L-Atelier/GAIA/reference/APOLLO_extract_isoform_sequences.md),
inside an isolated basilisk Python environment (`.gaia_iupred2a_env`).

## Usage

``` r
APOLLO_run_iupred2a(
  switch_result,
  iupred2a_dir,
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

- iupred2a_dir:

  Character. Path to the downloaded IUPred2A repository's `iupred2a/`
  subdirectory – a flat checkout containing `iupred2a.py`,
  `iupred2a_lib.py`, and the `data/` folder of energy
  matrices/histograms.

- output_dir:

  Character. Directory the combined result file (and per-isoform split
  input FASTAs) are written to (created if absent).

- aa_fasta:

  Character or `NULL`. Amino acid FASTA path. Default `NULL` uses
  `switch_result$sequence_files$aa_fasta`.

- verbose:

  Logical. Print progress, including a running isoform counter – this
  tool is invoked once per isoform, not once for the whole FASTA.
  Default: `TRUE`.

## Value

Invisibly, the path to a combined result file compatible with
`IsoformSwitchAnalyzeR::analyzeIUPred2A(pathToIUPred2AresultFile = ...)`.

## Details

`iupred2a.py` is strictly single-sequence: its own `read_seq()` does not
parse FASTA records at all – it discards `>` header lines and
concatenates every remaining line into one string, so handing it a
multi-isoform FASTA directly would silently merge every isoform's
sequence into one meaningless combined sequence rather than erroring
(confirmed by reading `iupred2a_lib.py` directly). This function
therefore splits `aa_fasta` into one temporary single-sequence file per
isoform (under `output_dir/split_input/`) and invokes
`iupred2a.py -a <seqfile> long` once per isoform (`-a` for ANCHOR2
binding-site prediction – always enabled, since `analyzeIUPred2A()`'s
default `annotateBindingSites=TRUE` requires the 4-column ANCHOR2 output
or it errors; `long` is the disorder mode IUPred2A's own webserver
instructions specify as standard).

`iupred2a.py`'s raw stdout has no sequence identifier in it at all (just
a citation banner, a column header, and per-residue rows) – but
`analyzeIUPred2A()`'s parser (confirmed by inspecting its source
directly) expects a combined file where each isoform's block is preceded
by its own `>isoform_id` line, matching the batch format the IUPred2A
**webserver** produces (which is what that importer was actually built
against, not the local script's native single-run output). This function
reconstructs that shape: after each per-isoform run, it prepends the
isoform's own FASTA header before concatenating all blocks into one
combined result file.

`iupred2a.py` computes its own data directory as the location of the
script itself (`os.path.dirname(os.path.realpath(__file__))`), and
Python inserts a script's own directory at `sys.path[0]` on direct
invocation – so unlike
[`APOLLO_run_deeptmhmm()`](https://ylefol.github.io/L-Atelier/GAIA/reference/APOLLO_run_deeptmhmm.md),
no [`setwd()`](https://rdrr.io/r/base/getwd.html) into `iupred2a_dir` is
needed; the script is called by its absolute path from any working
directory.

`iupred2a.py` requires no external Python libraries (confirmed by
reading its README and source directly), so `.gaia_iupred2a_env` only
provisions a bare Python interpreter – no manual pip-install step is
needed, unlike
[`APOLLO_run_signalp()`](https://ylefol.github.io/L-Atelier/GAIA/reference/APOLLO_run_signalp.md)/[`APOLLO_run_deeptmhmm()`](https://ylefol.github.io/L-Atelier/GAIA/reference/APOLLO_run_deeptmhmm.md)/
[`APOLLO_run_deeploc2()`](https://ylefol.github.io/L-Atelier/GAIA/reference/APOLLO_run_deeploc2.md).
