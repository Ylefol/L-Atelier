# Infer transcription factor activity

Convenience wrapper for TF activity inference. Fetches the appropriate
network and runs activity inference.

## Usage

``` r
ARTEMIS_infer_tf_activity(
  mat,
  method,
  database = "collectri",
  organism = "human",
  dorothea_levels = c("A", "B", "C", "D", "E"),
  minsize = 5,
  verbose = TRUE,
  ...
)
```

## Arguments

- mat:

  Numeric matrix of gene expression data. Genes as rows, samples as
  columns.

- method:

  Character. Statistical method to use. REQUIRED - no default.

- database:

  Character. TF-target database to use. Options: "collectri" (default),
  "dorothea".

- organism:

  Character. Organism: "human" or "mouse". Default: "human".

- dorothea_levels:

  Character vector. If database = "dorothea", which confidence levels to
  include. Default: c("A", "B", "C", "D", "E").

- minsize:

  Integer. Minimum targets per TF. Default: 5.

- verbose:

  Logical. Print progress. Default: TRUE.

- ...:

  Additional arguments passed to ARTEMIS_run_decoupler().

## Value

A decoupler_result object with TF activity scores.

## Examples

``` r
if (FALSE) { # \dontrun{
# Using CollecTRI (default)
tf_activity <- ARTEMIS_infer_tf_activity(expr_matrix, method = "ulm")

# Using DoRothEA with high-confidence regulons
tf_activity <- ARTEMIS_infer_tf_activity(
  expr_matrix,
  method = "mlm",
  database = "dorothea",
  dorothea_levels = c("A", "B")
)

} # }
```
