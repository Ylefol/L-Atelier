# Run CPAT Coding Potential Prediction

Runs CPAT (Coding-Potential Assessment Tool) against the nucleotide
FASTA from
[`APOLLO_extract_isoform_sequences()`](https://ylefol.github.io/L-Atelier/GAIA/reference/APOLLO_extract_isoform_sequences.md),
inside an isolated basilisk Python environment (`.gaia_cpat_env`). No
manual Python/conda setup is required beyond having internet access on
first use (basilisk builds the environment then).

## Usage

``` r
APOLLO_run_cpat(
  switch_result,
  hexamer_file,
  logit_model_file,
  output_dir,
  nt_fasta = NULL,
  verbose = TRUE
)
```

## Arguments

- switch_result:

  An `artemis_isoform_switch` object, normally one that has already been
  through
  [`APOLLO_extract_isoform_sequences()`](https://ylefol.github.io/L-Atelier/GAIA/reference/APOLLO_extract_isoform_sequences.md).

- hexamer_file:

  Character. Path to a CPAT hexamer frequency table (species-specific;
  CPAT publishes prebuilt tables for human/mouse/fly/ zebrafish, or
  generate one with CPAT's `make_hexamer_tab.py`).

- logit_model_file:

  Character. Path to the matching CPAT logit model file (`.RData`).

- output_dir:

  Character. Directory for CPAT's raw output and the reformatted result
  file (created if absent).

- nt_fasta:

  Character or `NULL`. Nucleotide FASTA path. Default `NULL` uses
  `switch_result$sequence_files$nt_fasta`.

- verbose:

  Logical. Print progress and CPAT's own console output. Default:
  `TRUE`.

## Value

Invisibly, the path to a reformatted result file compatible with
`ARTEMIS_isoform_switch_consequences(cpat_file = ...)`.

## Details

CPAT's own output columns
(`ID, mRNA, ORF_strand, ORF_frame, ORF_start, ORF_end, ORF, Fickett, Hexamer, Coding_prob`)
do not match what
[`IsoformSwitchAnalyzeR::analyzeCPAT()`](https://rdrr.io/pkg/IsoformSwitchAnalyzeR/man/analyzeCPAT.html)
expects to parse (a 5- or 8-column legacy format) – confirmed by
inspecting `analyzeCPAT()`'s source directly, and a known point of
confusion in the community (e.g. <https://www.biostars.org/p/9531584/>).
This function reformats CPAT's raw output into the exact column
set/order `analyzeCPAT()` requires before returning the path, so no
manual reformatting is needed downstream.
