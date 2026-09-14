# Apollo - Isoform Switch Consequence Tool Wrappers

Runs external protein/transcript annotation tools (CPAT, SignalP, Pfam
via `pfam_scan.pl`, DeepTMHMM) against isoforms from an
`artemis_isoform_switch` object, producing result files consumable by
[`ARTEMIS_isoform_switch_consequences()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_isoform_switch_consequences.md).
All four run inside isolated basilisk environments (see
`GAIA/R/basilisk.R`) – CPAT and SignalP are Python tools; Pfam is
packaged via bioconda's `pfam_scan` recipe, which basilisk installs the
same way regardless of it being a Perl/native-binary tool rather than a
Python one; DeepTMHMM is an academically licensed PyTorch model run from
a user-supplied local checkout rather than an installed executable.
Extract Isoform Nucleotide/Amino Acid Sequences

Extracts nucleotide and amino acid FASTA sequences for the isoforms in
an `artemis_isoform_switch` object, via
[`IsoformSwitchAnalyzeR::extractSequence()`](https://rdrr.io/pkg/IsoformSwitchAnalyzeR/man/extractSequence.html).
This is the required first step before running any of
[`APOLLO_run_cpat()`](https://ylefol.github.io/L-Atelier/GAIA/reference/APOLLO_run_cpat.md),
[`APOLLO_run_signalp()`](https://ylefol.github.io/L-Atelier/GAIA/reference/APOLLO_run_signalp.md),
or
[`APOLLO_run_pfam()`](https://ylefol.github.io/L-Atelier/GAIA/reference/APOLLO_run_pfam.md)
– the nucleotide FASTA feeds CPAT, the amino acid FASTA feeds SignalP
and Pfam.

## Usage

``` r
APOLLO_extract_isoform_sequences(
  switch_result,
  genome_object,
  output_dir,
  only_switching_genes = TRUE,
  verbose = TRUE
)
```

## Arguments

- switch_result:

  An `artemis_isoform_switch` object from
  [`ARTEMIS_isoform_switch()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_isoform_switch.md).

- genome_object:

  A `BSgenome` object matching the genome used to build `switch_result`
  (e.g.
  [`BSgenome.Hsapiens.UCSC.hg38::Hsapiens`](https://rdrr.io/pkg/BSgenome.Hsapiens.UCSC.hg38/man/package.html)
  after loading that annotation package). This is a hard requirement of
  [`IsoformSwitchAnalyzeR::extractSequence()`](https://rdrr.io/pkg/IsoformSwitchAnalyzeR/man/extractSequence.html)
  – a plain genome FASTA path is not accepted.

- output_dir:

  Character. Directory the FASTA files are written to (created if
  absent).

- only_switching_genes:

  Logical. If `TRUE` (default), only extract sequences for genes already
  called significant in `switch_result` (using its own
  `alpha`/`dIF_cutoff`). Set `FALSE` to extract for every tested gene.

- verbose:

  Logical. Print progress. Default: `TRUE`.

## Value

The input `artemis_isoform_switch` object with two additional fields:

- sequence_files:

  List with `nt_fasta` and `aa_fasta` paths, used as defaults by
  [`APOLLO_run_cpat()`](https://ylefol.github.io/L-Atelier/GAIA/reference/APOLLO_run_cpat.md),
  [`APOLLO_run_signalp()`](https://ylefol.github.io/L-Atelier/GAIA/reference/APOLLO_run_signalp.md),
  and
  [`APOLLO_run_pfam()`](https://ylefol.github.io/L-Atelier/GAIA/reference/APOLLO_run_pfam.md).

- switch_list:

  Updated switchAnalyzeRlist with the sequences also cached internally
  (`ntSequence`/`aaSequence`).
