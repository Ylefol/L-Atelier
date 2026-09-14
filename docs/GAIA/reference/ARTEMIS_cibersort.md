# Run CIBERSORT immune cell deconvolution

Estimates cell type proportions from bulk gene expression data using
nu-support vector regression. Wraps the original CIBERSORT source code
with a GAIA-compatible interface that accepts R objects directly.

## Usage

``` r
ARTEMIS_cibersort(
  mixture,
  cibersort_path,
  sig_matrix,
  output_dir = "data",
  mixture_name = "mixture",
  perm = 100,
  QN = TRUE,
  verbose = TRUE
)
```

## Arguments

- mixture:

  Expression matrix (genes as rows, samples as columns) or a file path
  to a tab-delimited mixture file. Gene names should be in rownames (or
  first column for files). Can also be an artemis_norm object, in which
  case norm_counts is extracted automatically. Regardless of input type,
  the mixture is always reordered to match `sig_matrix`'s gene order
  (see Details) and rewritten to
  `<output_dir>/cibersort/<mixture_name>.txt` – a file-path input is
  read and rewritten too, not used in place.

- cibersort_path:

  Path to the original CIBERSORT.R source file.

- sig_matrix:

  Path to the signature matrix file (e.g., LM22.txt for 22 immune cell
  types).

- output_dir:

  Directory where the mixture file and CIBERSORT results will be saved.
  A "cibersort" subdirectory is created within it. Default:
  "data/cibersort".

- mixture_name:

  Character. Base name for the saved mixture file (without extension).
  Default: "mixture". The file is saved as
  `<output_dir>/cibersort/<mixture_name>.txt`.

- perm:

  Number of permutations for p-value calculation. Default 100

- QN:

  Logical. Apply quantile normalization. Default TRUE for microarray
  data. Set to FALSE for RNA-seq data.

- verbose:

  Logical. Print progress messages. Default TRUE.

## Value

S3 object of class "artemis_cibersort" with components:

- proportions:

  Matrix of cell type proportions (samples x cell types). Rows sum to 1.

- p_values:

  Numeric vector of per-sample p-values (NULL if perm=0)

- correlations:

  Numeric vector of Pearson correlations per sample

- rmse:

  Numeric vector of RMSE per sample

- raw:

  The full raw CIBERSORT output matrix

- metadata:

  List with parameters and run info

## Details

CIBERSORT uses nu-SVR to deconvolve bulk expression into cell type
proportions using a reference signature matrix. The standard LM22 matrix
contains signatures for 22 human immune cell types.

For RNA-seq data, set QN=FALSE. Quantile normalization is designed for
microarray data and can distort RNA-seq distributions.

The mixture file is saved to `<output_dir>/cibersort/<mixture_name>.txt`
for reproducibility and re-use. CIBERSORT's side-effect output
(CIBERSORT-Results.txt) is also written to the same directory.

**Gene order:** CIBERSORT.R's own gene-intersection step filters the
signature matrix and mixture down to their shared genes but does not
sort one to match the other, and the downstream SVR regression pairs
them up positionally (no gene-name join). If the two input files use
different row orders – true unless someone deliberately pre-sorts both
identically, e.g. LM22.txt is alphabetical but many other signature
matrices are not – this silently misaligns genes between the reference
and the sample and can wreck the fit (near-zero correlation, RMSE at or
above the standardized target's own SD) while looking like a
data-quality problem rather than a misalignment bug.
`ARTEMIS_cibersort()` reorders the mixture to match `sig_matrix`'s gene
order before writing it out, so this can't happen regardless of the
input files' original ordering.

## Examples

``` r
if (FALSE) { # \dontrun{
result <- ARTEMIS_cibersort(
  mixture        = norm_counts,
  cibersort_path = "data/CIBERSORT.R",
  sig_matrix     = "data/LM22.txt",
  output_dir     = "results",
  mixture_name   = "my_experiment",
  perm           = 100,
  QN             = FALSE  # RNA-seq
)

# View proportions
head(result$proportions)

# From an artemis_norm object
result <- ARTEMIS_cibersort(
  mixture        = my_norm,
  cibersort_path = "data/CIBERSORT.R",
  sig_matrix     = "data/LM22.txt",
  output_dir     = "results",
  QN             = FALSE
)

} # }
```
