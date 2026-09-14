# Run decoupleR activity inference

Core function for inferring biological activities (TF or pathway) from
gene expression data using a specified statistical method.

## Usage

``` r
ARTEMIS_run_decoupler(
  mat,
  network,
  method,
  minsize = 5,
  norm_method = "log2cpm",
  verbose = TRUE,
  ...
)
```

## Arguments

- mat:

  Numeric matrix of gene expression data (genes x samples), an
  `artemis_norm` object from
  [`ARTEMIS_normalize_counts()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_normalize_counts.md),
  or raw integer count matrix. Rownames must be gene identifiers
  matching the network.

- network:

  Prior knowledge network (from APOLLO_get\_\* functions or custom).
  Must have columns: 'source' (TF/pathway), 'target' (gene), and 'mor'
  or 'weight'.

- method:

  Character. Statistical method to use. REQUIRED - no default. Options:
  "ulm", "mlm", "viper", "wsum", "wmean", "udt", "mdt", "ora", "gsva",
  "aucell".

- minsize:

  Integer. Minimum number of targets per source to include. Default: 5.

- norm_method:

  Character. Normalization to apply to raw counts before inference.
  `"log2cpm"` (default): log2(CPM + 1). `"none"`: pass matrix through
  as-is (use when already normalized). If `mat` is an `artemis_norm`
  object, DESeq2-normalized counts are extracted and log2-transformed
  automatically, regardless of this parameter.

- verbose:

  Logical. Print progress messages. Default: TRUE.

- ...:

  Additional arguments passed to the specific method function.

## Value

A list with class "decoupler_result" containing:

- activities:

  Matrix of activity scores (sources x samples)

- pvalues:

  Matrix of p-values (if available for method)

- results_long:

  Full results in long format (from decoupleR)

- method:

  Method used

- network_info:

  Information about the network used

- n_sources:

  Number of sources (TFs/pathways) with results

- n_samples:

  Number of samples

## Details

This function requires an explicit method choice to ensure users
understand which statistical approach they are using. For method
recommendations, use
[`ARTEMIS_decoupler_list_methods()`](https://ylefol.github.io/L-Atelier/GAIA/reference/ARTEMIS_decoupler_list_methods.md).

DecoupleR methods expect log-normalized continuous expression values,
not raw counts. Use `norm_method = "log2cpm"` for raw counts, or pass an
`artemis_norm` object to reuse DESeq2 size factors from a prior
normalization.

## Examples

``` r
if (FALSE) { # \dontrun{
# Get network
network <- APOLLO_get_collectri()

# Run ULM (recommended)
result <- ARTEMIS_run_decoupler(expr_matrix, network, method = "ulm")

# Run MLM
result <- ARTEMIS_run_decoupler(expr_matrix, network, method = "mlm")

} # }
```
