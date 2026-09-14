# Installation

``` r
install.packages("devtools")
devtools::install("path/to/L-Atelier/CYAN")
library(CYAN)
```

R \>= 4.3 is required. `basilisk` and `reticulate` are in `Imports` and
installed automatically.

No manual conda setup is needed for RNA editing:
[`CYAN_run_reditools()`](https://ylefol.github.io/L-Atelier/CYAN/reference/CYAN_run_reditools.md)
runs REDItools2 through a `basilisk`-managed Python environment (Python
3.10.14, `pysam`, `sortedcontainers`, `psutil`, `netifaces`, pinned
versions) that `basilisk` creates automatically on first use. REDItools2
itself is vendored in `inst/python/reditools.py` (the original script,
unmodified).
