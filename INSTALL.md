# Installing L-Atelier

This page covers installing HORIZON, GAIA and CAULDRON. For what each package does, see the
[README](README.md). Each package also has its own README with the external tools it calls.

## 1. Requirements

| | |
|---|---|
| **R** | >= 4.4.0. The toolkit targets the **current Bioconductor release only** (older releases are not supported). Bioconductor 3.24 (November 2026) needs R >= 4.6.0, so install R 4.6.x. |
| **Bioconductor** | Chosen automatically by `BiocManager` to match your R version. |
| **Internet access** | To CRAN, Bioconductor and `github.com`. |
| **A compiler** | Some Bioconductor packages are built from source (see [Time](#5-how-long-it-takes)). Linux: `build-essential` (gcc, g++, gfortran). macOS: Xcode command line tools and gfortran. Windows: [Rtools](https://cran.r-project.org/bin/windows/Rtools/). |

### Installing R

Keeping your old R and adding the new one side by side is the least disruptive route.
[`rig`](https://github.com/r-lib/rig) does this on Linux, macOS and Windows:

```bash
rig add release      # installs the newest R
rig list             # both versions are listed
rig default 4.6.2    # choose the one to use
```

Each R version gets its own package library, so nothing from the old R carries over and nothing
breaks. You will reinstall your packages once for the new R.

> **Install packages as yourself, not with `sudo`.** With a current R installed this way,
> packages go to a personal library in your home directory (for example
> `~/R/x86_64-pc-linux-gnu-library/4.6`) and need no administrator rights. A package installed
> from `sudo R` goes to *root's* library, which your own R does not search; it then looks installed
> (`BiocManager` says "already up to date") but `library()` cannot find it.

## 2. Install

From R (4.4.0 or newer):

```r
source("https://raw.githubusercontent.com/Ylefol/L-Atelier/master/install.R")
```

The installer checks your R and Bioconductor versions, lets you choose which packages to install,
installs them with their required dependencies, offers optional dependencies by feature (see
[section 4](#4-optional-features)), checks that Bioconductor is consistent, and reports what it
could not do for you.

### Non-interactive install

Set the options before sourcing the script (useful for servers and scripts):

```r
options(latelier.packages = c("HORIZON", "GAIA"),   # any of HORIZON, GAIA, CAULDRON
        latelier.optional = "none")                 # "none", "all", or group names (section 4)
source("https://raw.githubusercontent.com/Ylefol/L-Atelier/master/install.R")
```

| Option | Meaning | Default |
|---|---|---|
| `latelier.packages` | packages to install | ask |
| `latelier.optional` | optional dependency groups | ask |
| `latelier.ref` | Git branch, tag or commit to install from | `"master"` |
| `latelier.upgrade` | passed to `remotes` | `"never"` |

To keep a log of the whole install, run it with `Rscript` and `tee`:

```bash
Rscript -e 'options(latelier.packages = c("HORIZON","GAIA"), latelier.optional = "none"); source("install.R")' 2>&1 | tee install.log
```

### Installing one package by hand

```r
install.packages("remotes")
remotes::install_github("Ylefol/L-Atelier", subdir = "GAIA", dependencies = TRUE)
```

Bioconductor dependencies resolve only if `BiocManager` is installed and your repositories point
at Bioconductor (`BiocManager::repositories()`), which is what `install.R` sets up. The installer
is the recommended route.

## 3. What the installer does not do

- It does **not** install system libraries.
- It does **not** create conda environments (HORIZON needs them for its external tools).
- It does **not** fetch software that must be obtained from its authors (for example CIBERSORT; see
  `GAIA/README.md`).

### System libraries (Linux)

On a Linux machine with no development libraries, compiling some R packages needs a few system
libraries. Our CI installs these three before running the installer on Ubuntu:

```bash
sudo apt-get install -y libcurl4-openssl-dev libssl-dev libxml2-dev
```

**This list is not proven to be minimal or complete.** It is what the CI job installs, and the CI
runner image already contains many other libraries and build tools. The install has been run on
Ubuntu 24.04 with those tools already present. It has **not** been tried on a bare machine or
container. If an install fails on a missing header file (the error names it, for example
`xxx.h: No such file or directory`), install the matching `-dev` package and retry. If you hit one,
please open an issue so this table can be completed.

### HORIZON external tools

HORIZON calls command-line tools from a conda environment. RNA-seq quality control, alignment and
counting use R packages and need no conda environment. The chromatin pipeline (ATAC, ChIP,
CUT&TAG) needs a `horizon_cli` environment with samtools, bedtools, macs3, bowtie2 and deeptools.
The PARSE single-cell steps need a separate environment. Exact commands are in
[`HORIZON/README.md`](HORIZON/README.md) (steps 2 to 4).

## 4. Optional features

Required (`Imports`) dependencies are always installed. Everything else is grouped by feature and
installed only if you ask for it. Choose group names for `latelier.optional`:

| Group | For | Packages |
|---|---|---|
| `enrichment` | GAIA, CAULDRON | gprofiler2, msigdbr, clusterProfiler, fgsea, rrvgo, org.Hs.eg.db, ... |
| `salmon` | HORIZON | txdbmaker, GenomicFeatures, Biostrings, Rsamtools |
| `chromatin` | GAIA | ChIPseeker, EnrichedHeatmap, BSgenome, ... |
| `batch_norm` | GAIA | sva, limma, preprocessCore, vsn, RUVSeq, imputeLCMD, mice |
| `integration` | GAIA | MOFA2, ConsensusClusterPlus, PMA, grpreg, FactoMineR, ... |
| `network_viz` | GAIA | igraph, ggraph, visNetwork, plotly, ggalluvial, ... |
| `signaling` | GAIA | decoupleR, OmnipathR, CARNIVAL, lpSolve |
| `isoform` | GAIA | IsoformSwitchAnalyzeR, reticulate |
| `sc_core` | CAULDRON | Seurat, harmony, zellkonverter, DropletUtils, celldex, ... |
| `sc_advanced` | CAULDRON | monocle3, CellChat, BPCells (GitHub only) |
| `sc_stats` | CAULDRON | speckle, orthogene, fgsea, gprofiler2, msigdbr |

Notes, as of the last test (November 2026):

- **`signaling` is currently blocked**: `decoupleR` is missing from the Bioconductor 3.24
  repositories, so it cannot be installed from Bioconductor until Bioconductor restores it.
- **Only the installs without optional groups have been tested** (see below). The optional groups
  have not been tested, and `sc_advanced` (GitHub-only packages that compile from source) is the
  most likely to need attention.
- A function that needs an optional package it cannot find stops with a message naming the package
  to install (for example `imputeLCMD` for mass spec imputation, which also needs `pcaMethods`
  from Bioconductor).

## 5. How long it takes

Measured on Ubuntu 24.04, a fast connection, R 4.6.2 and Bioconductor 3.24, with a fresh package
library and no optional groups:

| Install | Time |
|---|---|
| HORIZON + GAIA | about 9 minutes |
| CAULDRON (after HORIZON and GAIA) | about 11 minutes |

CI runs of the HORIZON + GAIA install (the whole job, including setting up R and running the
integration tests, on GitHub-hosted runners) took about 13 minutes on Ubuntu, 13 on Windows and 7
on macOS.

These are single measurements, not benchmarks. Much of the time is compiling Bioconductor packages
from source (for example `Rsubread`, `Rhtslib`, `scrapper`, `scran`, `scater`).

## 6. What has been tested

| | Status |
|---|---|
| `install.R` for HORIZON + GAIA on Ubuntu, macOS and Windows | run in CI on every push; passes |
| HORIZON -> GAIA hand-off (count files -> aggregation -> differential analysis) | run in CI after the install; passes |
| Package checks (`R CMD check`) for HORIZON, GAIA | Ubuntu, macOS, Windows, on Bioconductor 3.24 |
| Package check for CAULDRON | Ubuntu only |
| Cold install of HORIZON, GAIA and CAULDRON on a fresh R 4.6.2 library | Ubuntu 24.04, one machine |
| Bare machine / container, other Linux distributions, optional groups | **not tested** |

## 7. Troubleshooting

- **A package is installed but `library()` cannot find it.** Check which R you are running and
  where it looks:
  ```r
  R.version.string
  .libPaths()
  find.package("the_package")
  ```
  A different R version than you expect (an IDE may use the system R) or an install made with
  `sudo` are the usual causes. Install as your normal user, with the R you will actually use.
- **IDE console problems.** Some IDE consoles fail on downloads or `system2()` calls where the same
  code works from a terminal. If an install or a download fails in a way that looks
  environmental, run it with `Rscript` from a terminal.
- **`BiocManager::valid()` reports problems.** Run `BiocManager::install()` and follow its
  suggestions; the toolkit assumes a consistent installation of the current release.
- **A newer R than the current Bioconductor supports.** `BiocManager` picks the Bioconductor
  release that matches your R; if none exists yet it tells you.
- **Installing a development branch.** Set `options(latelier.ref = "dev")` before sourcing the
  installer.
