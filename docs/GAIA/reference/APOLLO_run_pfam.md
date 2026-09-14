# Run Pfam Protein Domain Annotation

Runs `pfam_scan.pl` (a Perl wrapper around HMMER's `hmmscan`) against
the amino acid FASTA from
[`APOLLO_extract_isoform_sequences()`](https://ylefol.github.io/L-Atelier/GAIA/reference/APOLLO_extract_isoform_sequences.md),
inside an isolated basilisk environment (`.gaia_pfam_env`).
`pfam_scan.pl` is packaged on bioconda as `pfam_scan`, which pulls in
HMMER and its own Perl module dependencies as part of the conda recipe –
so no manual PfamScan.tar.gz/CPAN setup is needed.

## Usage

``` r
APOLLO_run_pfam(
  switch_result,
  pfam_dir,
  output_dir,
  aa_fasta = NULL,
  nproc = 4,
  verbose = TRUE
)
```

## Arguments

- switch_result:

  An `artemis_isoform_switch` object, normally one that has already been
  through
  [`APOLLO_extract_isoform_sequences()`](https://ylefol.github.io/L-Atelier/GAIA/reference/APOLLO_extract_isoform_sequences.md).

- pfam_dir:

  Character. Path to a directory containing the Pfam-A HMM database,
  pre-pressed with `hmmpress` (i.e. containing `Pfam-A.hmm`,
  `Pfam-A.hmm.dat`, and the `hmmpress`- generated index files).

- output_dir:

  Character. Directory for Pfam's output (created if absent).

- aa_fasta:

  Character or `NULL`. Amino acid FASTA path. Default `NULL` uses
  `switch_result$sequence_files$aa_fasta`.

- nproc:

  Integer. Number of CPUs for `pfam_scan.pl`'s `-cpu` option. Default:
  `4`.

- verbose:

  Logical. Print progress and `pfam_scan.pl`'s own console output.
  Default: `TRUE`.

## Value

Invisibly, the path to the Pfam result file, directly compatible with
`ARTEMIS_isoform_switch_consequences(pfam_file = ...)`.

## Details

[`IsoformSwitchAnalyzeR::analyzePFAM()`](https://rdrr.io/pkg/IsoformSwitchAnalyzeR/man/analyzePFAM.html)
expects `pfam_scan.pl`'s specific output format (confirmed against the
package's own bundled example result file), not raw
`hmmscan --domtblout` output – this is why the wrapper runs
`pfam_scan.pl` rather than calling `hmmscan` directly, even though the
latter would be sufficient to scan against Pfam-A on its own.
