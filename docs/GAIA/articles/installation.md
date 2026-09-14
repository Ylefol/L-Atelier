# Installation

GAIA is an R package installed directly from the repository. R \>= 4.1.0
is required.

``` r
install.packages("devtools")
devtools::install("path/to/L-Atelier/GAIA")
library(GAIA)
```

## Optional dependencies

Several modules have optional dependencies loaded on demand (listed in
`Suggests`). Key ones:

- **Bioconductor** — `GenomicRanges`, `IRanges`, `ChIPseeker`,
  `clusterProfiler`, `Rsubread`, `rtracklayer`, etc. Install via
  [`BiocManager::install()`](https://bioconductor.github.io/BiocManager/reference/install.html).
- **Batch correction** — `sva` (ComBat), `limma` (removeBatchEffect).
- **Annotation** — `org.Hs.eg.db`, `AnnotationDbi`.
- **CIBERSORT** — Source code must be obtained separately from the
  authors. See the
  [Reference](https://ylefol.github.io/L-Atelier/GAIA/reference/index.md)
  page.
- **Quantile normalization** — `preprocessCore`.

If Bioconductor packages are unexpectedly unavailable after
installation, run
[`BiocManager::valid()`](https://bioconductor.github.io/BiocManager/reference/valid.html)
to identify and fix any version mismatches.
