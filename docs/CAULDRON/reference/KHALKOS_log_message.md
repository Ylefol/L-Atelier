# Verbose-aware logging helper

Emits a message, warning, or error depending on `level`. At the `"info"`
level the message is only printed when `verbose = TRUE`.

## Usage

``` r
KHALKOS_log_message(msg, verbose = TRUE, level = c("info", "warning", "error"))
```

## Arguments

- msg:

  Character. Message text.

- verbose:

  Logical. Print info-level messages. Default `TRUE`.

- level:

  Character. One of `"info"` (default), `"warning"`, or `"error"`.
  `"error"` always throws; `"warning"` always warns.

## Value

Invisibly `NULL`.
