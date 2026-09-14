# Load HOMER results from an existing output directory

Parses HOMER findMotifsGenome.pl output without re-running. Useful for
loading previously computed results.

## Usage

``` r
APOLLO_load_homer_results(homer_dir, verbose = TRUE)
```

## Arguments

- homer_dir:

  Path to HOMER output directory (containing knownResults.txt).

- verbose:

  Logical. Print loading messages. Default TRUE.

## Value

S3 object of class "homer_motif" with components:

- known:

  Data.frame of known motif enrichment results

- denovo:

  Data.frame of de novo motifs (NULL if not available)

- output_dir:

  Path to HOMER output directory

- metadata:

  List with n_target_seqs, n_background_seqs
