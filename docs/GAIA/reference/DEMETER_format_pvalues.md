# Format p-values for display

Converts p-values to formatted strings, optionally in scientific
notation. Non-significant values can be replaced with empty strings or
NA.

## Usage

``` r
DEMETER_format_pvalues(
  pvalues,
  format = "scientific",
  digits = 2,
  sig_threshold = 0.05,
  hide_ns = FALSE,
  ns_string = ""
)
```

## Arguments

- pvalues:

  Numeric vector or matrix of p-values.

- format:

  Output format: "scientific" (e.g., "1.23e-04"), "stars" (\*, \*\*,
  \*\*\*), "decimal" (e.g., "0.001"), or "hybrid" (decimal if \>= 0.001,
  scientific otherwise). Default: "scientific".

- digits:

  Number of significant digits. Default: 2.

- sig_threshold:

  P-value threshold for significance. Values above this become empty
  strings (if hide_ns = TRUE). Default: 0.05.

- hide_ns:

  Logical. Replace non-significant p-values with empty string. Default:
  FALSE.

- ns_string:

  String to use for non-significant values when hide_ns = TRUE. Default:
  "".

## Value

Character vector or matrix of formatted p-values (same shape as input).

## Examples

``` r
if (FALSE) { # \dontrun{
pvals <- c(0.001, 0.023, 0.15, 0.0001)
DEMETER_format_pvalues(pvals)
DEMETER_format_pvalues(pvals, format = "stars")
DEMETER_format_pvalues(pvals, hide_ns = TRUE)

} # }
```
