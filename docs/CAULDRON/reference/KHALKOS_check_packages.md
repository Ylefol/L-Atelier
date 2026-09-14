# Check optional package availability

Tests whether each package in `packages` can be loaded with
[`requireNamespace()`](https://rdrr.io/r/base/ns-load.html). Missing
packages are reported together with an installation hint.

## Usage

``` r
KHALKOS_check_packages(packages, install_hint = NULL, error = TRUE)
```

## Arguments

- packages:

  Character vector of package names to check.

- install_hint:

  Character scalar. Custom install instructions printed in the
  error/warning. When `NULL` (default), a generic
  [`BiocManager::install()`](https://bioconductor.github.io/BiocManager/reference/install.html)
  call is suggested.

- error:

  Logical. If `TRUE` (default), throw an error for missing packages. If
  `FALSE`, issue a warning instead.

## Value

Invisibly `TRUE` if all packages are available; invisibly `FALSE`
otherwise (only when `error = FALSE`).
