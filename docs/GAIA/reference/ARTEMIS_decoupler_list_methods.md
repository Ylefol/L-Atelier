# List available decoupleR methods

Prints information about available statistical methods for activity
inference.

## Usage

``` r
ARTEMIS_decoupler_list_methods(recommended_only = FALSE)
```

## Arguments

- recommended_only:

  Logical. If TRUE, only show recommended methods. Default: FALSE.

## Value

Character vector of method names (invisibly).

## Examples

``` r
if (FALSE) { # \dontrun{
ARTEMIS_decoupler_list_methods()
ARTEMIS_decoupler_list_methods(recommended_only = TRUE)

} # }
```
