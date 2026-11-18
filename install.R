# ==============================================================================
# L-Atelier installer
# ==============================================================================
# Installs HORIZON, GAIA and/or CAULDRON from GitHub, together with their
# required (Imports) dependencies, then offers optional dependencies by feature.
#
# Usage (interactive R session):
#   source("https://raw.githubusercontent.com/Ylefol/L-Atelier/master/install.R")
#
# Non-interactive / scripted use: set options BEFORE sourcing, e.g.
#   options(latelier.packages = c("HORIZON", "GAIA"),
#           latelier.optional = "none")   # "none", "all", or group names
#   source("install.R")
#
# Available options (all optional):
#   latelier.packages  chr   which packages to install (default: ask)
#   latelier.optional  chr   optional dependency groups (default: ask)
#   latelier.ref       chr   git ref to install from (default: "master")
#   latelier.upgrade   chr   passed to remotes (default: "never")
#
# This script does NOT: create conda environments, install system libraries, or
# fetch software that must be obtained from its authors (e.g. CIBERSORT). It
# reports what is missing and points to the relevant README instead.
# ==============================================================================

# --- Configuration ------------------------------------------------------------
LATELIER_REPO  <- "Ylefol/L-Atelier"
LATELIER_PKGS  <- c("HORIZON", "GAIA", "CAULDRON")

# Policy: the toolkit targets the CURRENT Bioconductor release only (no version
# matrix). Hard minimum is R 4.4.0 / Bioconductor 3.19, because makeTxDbFromGFF()
# was moved to txdbmaker in that release. CI runs on the current R release; the
# Bioconductor 3.24 release (Nov 2026) requires R >= 4.6.0, so warn below that.
LATELIER_MIN_R    <- "4.4.0"
LATELIER_TESTED_R <- "4.6.0"   # below this: warn, do not stop

# --- Optional dependency groups -----------------------------------------------
# Each entry: which L-Atelier package(s) it belongs to, a description shown to
# the user, and the packages to install. "github" entries are not on CRAN or
# Bioconductor and are installed with the GitHub repository tarball (no GitHub API).
# NOTE: GitHub repository names below should be verified before release.
LATELIER_OPTIONAL <- list(
  enrichment = list(
    pkgs = c("GAIA", "CAULDRON"),
    desc = "Functional enrichment (clusterProfiler, gprofiler2, fgsea, msigdbr, rrvgo, org.Hs.eg.db)",
    cran = c("gprofiler2", "msigdbr", "treemap"),
    bioc = c("clusterProfiler", "fgsea", "rrvgo", "GO.db", "org.Hs.eg.db", "AnnotationDbi")
  ),
  salmon = list(
    pkgs = "HORIZON",
    desc = "Salmon transcriptome index building (txdbmaker, GenomicFeatures, Biostrings, Rsamtools)",
    cran = character(0),
    bioc = c("txdbmaker", "GenomicFeatures", "Biostrings", "Rsamtools")
  ),
  chromatin = list(
    pkgs = "GAIA",
    desc = "Chromatin annotation and tracks (ChIPseeker, EnrichedHeatmap, txdbmaker, GenomicFeatures, BSgenome)",
    cran = character(0),
    bioc = c("ChIPseeker", "EnrichedHeatmap", "txdbmaker", "GenomicFeatures",
             "GenomeInfoDb", "Biostrings", "BSgenome", "Rsamtools")
  ),
  batch_norm = list(
    pkgs = "GAIA",
    desc = "Batch correction and normalization (sva, limma, preprocessCore, vsn, RUVSeq, imputeLCMD)",
    cran = c("imputeLCMD"),
    bioc = c("sva", "limma", "preprocessCore", "vsn", "RUVSeq")
  ),
  integration = list(
    pkgs = "GAIA",
    desc = "Multi-omics integration and ML (MOFA2, PMA, grpreg, FactoMineR, missMDA, e1071, uwot)",
    cran = c("PMA", "grpreg", "FactoMineR", "missMDA", "e1071", "uwot", "cluster", "nnet"),
    bioc = c("MOFA2")
  ),
  network_viz = list(
    pkgs = "GAIA",
    desc = "Network and interactive visualization (igraph, ggraph, visNetwork, plotly, pheatmap, patchwork)",
    cran = c("igraph", "ggraph", "visNetwork", "plotly", "htmlwidgets", "pheatmap",
             "patchwork", "ggrepel", "gridExtra", "RColorBrewer", "scales"),
    bioc = character(0)
  ),
  signaling = list(
    pkgs = "GAIA",
    desc = "Pathway/TF activity and signaling (decoupleR, OmnipathR, CARNIVAL, lpSolve)",
    cran = c("lpSolve"),
    bioc = c("decoupleR", "OmnipathR", "CARNIVAL")
  ),
  isoform = list(
    pkgs = "GAIA",
    desc = "Isoform switching (IsoformSwitchAnalyzeR; also needs basilisk/reticulate for some tools)",
    cran = c("reticulate"),
    bioc = c("IsoformSwitchAnalyzeR", "basilisk")
  ),
  sc_core = list(
    pkgs = "CAULDRON",
    desc = "Single-cell I/O and clustering (Seurat, harmony, zellkonverter, DropletUtils, bluster, celldex)",
    cran = c("Seurat", "SeuratObject", "harmony", "Rtsne", "uwot", "irlba", "igraph", "jsonlite"),
    bioc = c("zellkonverter", "DropletUtils", "LoomExperiment", "bluster", "celldex", "edgeR")
  ),
  sc_advanced = list(
    pkgs = "CAULDRON",
    desc = "Trajectory, cell-cell communication, large data (monocle3, CellChat, BPCells; GitHub-only)",
    cran = character(0),
    bioc = character(0),
    github = c("cole-trapnell-lab/monocle3", "jinworks/CellChat", "bnprks/BPCells/r")
  ),
  sc_stats = list(
    pkgs = "CAULDRON",
    desc = "Single-cell composition and enrichment (speckle, orthogene, fgsea, gprofiler2, msigdbr)",
    cran = c("gprofiler2", "msigdbr"),
    bioc = c("speckle", "orthogene", "fgsea")
  )
)

# --- Helpers ------------------------------------------------------------------
.lat_cat <- function(...) cat(..., sep = "")

.lat_ask <- function(prompt, default = "") {
  if (!interactive()) return(default)
  ans <- trimws(readline(prompt))
  if (nzchar(ans)) ans else default
}

.lat_yesno <- function(prompt, default = FALSE) {
  if (!interactive()) return(default)
  d <- if (default) "[Y/n]" else "[y/N]"
  ans <- tolower(trimws(readline(paste(prompt, d, ""))))
  if (!nzchar(ans)) return(default)
  ans %in% c("y", "yes")
}

# Install an R package from a GitHub repository by downloading the repository
# tarball. This deliberately avoids remotes::install_github(), which queries the
# GitHub API: unauthenticated API calls can fail or be rate limited on a machine
# without a GitHub token (observed in CI). Tarball downloads need no token and
# no git installation. `ref` may be a branch, tag or commit.
.lat_install_github_tarball <- function(repo, subdir = NULL, ref = "master",
                                        upgrade = "never") {
  url <- sprintf("https://github.com/%s/archive/%s.tar.gz", repo, ref)
  remotes::install_url(url, subdir = subdir, dependencies = NA,
                       upgrade = upgrade, quiet = FALSE)
}

.lat_need_pkg <- function(pkg, bioc = FALSE) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    .lat_cat("    Installing ", pkg, "...\n")
    if (bioc) BiocManager::install(pkg, ask = FALSE, update = FALSE)
    else install.packages(pkg)
  }
}

# --- Step 1: environment checks -------------------------------------------------
.lat_check_environment <- function() {
  .lat_cat("[L-ATELIER] Checking environment\n")
  rv <- getRversion()
  .lat_cat("    R version: ", as.character(rv), "\n")
  if (rv < LATELIER_MIN_R) {
    stop("R >= ", LATELIER_MIN_R, " is required (found ", rv, "). ",
         "Please upgrade R before installing.", call. = FALSE)
  }
  if (rv < LATELIER_TESTED_R) {
    .lat_cat("    WARNING: R < ", LATELIER_TESTED_R, " has not been tested with this toolkit.\n")
  }

  .lat_need_pkg("BiocManager")
  bv <- as.character(BiocManager::version())
  .lat_cat("    Bioconductor version: ", bv, "\n")
  # BiocManager picks the Bioconductor release matching the running R. If the
  # user has a pinned older release, tell them rather than silently continue.
  .lat_need_pkg("remotes")

  sysname <- Sys.info()[["sysname"]]
  .lat_cat("    Platform: ", sysname, "\n")
  if (sysname != "Linux") {
    .lat_cat("    NOTE: development and testing have been done on Ubuntu. ",
             "Check the CI status badge in the README for your platform.\n")
  }
  if (sysname == "Windows") {
    .lat_cat("    NOTE: Rtools is needed to compile packages from source ",
             "(https://cran.r-project.org/bin/windows/Rtools/).\n")
  }
  .lat_cat("\n")
  invisible(TRUE)
}

# --- Step 2: choose packages ----------------------------------------------------
.lat_choose_packages <- function() {
  chosen <- getOption("latelier.packages")
  if (!is.null(chosen)) {
    bad <- setdiff(chosen, LATELIER_PKGS)
    if (length(bad)) stop("Unknown package(s): ", paste(bad, collapse = ", "), call. = FALSE)
    return(chosen)
  }
  if (!interactive()) {
    stop("Non-interactive session: set options(latelier.packages = ...) first.", call. = FALSE)
  }
  .lat_cat("[L-ATELIER] Which packages would you like to install?\n")
  .lat_cat("    1: HORIZON   upstream processing (FASTQ -> counts/BAM/BigWig)\n")
  .lat_cat("    2: GAIA      bulk multi-omics analysis (RNA-seq, ATAC-seq, CUT&TAG)\n")
  .lat_cat("    3: CAULDRON  single-cell / single-nuclei RNA-seq analysis\n")
  .lat_cat("    4: all of the above\n")
  ans <- .lat_ask("    Enter numbers separated by commas (e.g. 1,2): ")
  idx <- suppressWarnings(as.integer(strsplit(ans, "[, ]+")[[1]]))
  idx <- idx[!is.na(idx)]
  if (4L %in% idx) return(LATELIER_PKGS)
  idx <- idx[idx %in% 1:3]
  if (!length(idx)) stop("No valid selection made.", call. = FALSE)
  LATELIER_PKGS[sort(unique(idx))]
}

# --- Step 3: install required dependencies and packages -------------------------
.lat_install_packages <- function(pkgs, ref, upgrade) {
  .lat_cat("[L-ATELIER] Installing from GitHub (", LATELIER_REPO, ", ref = ", ref, ")\n")
  # Make CRAN and Bioconductor repos visible to remotes so Imports resolve.
  old <- options(repos = BiocManager::repositories())
  on.exit(options(old), add = TRUE)

  for (p in pkgs) {
    .lat_cat("    ", p, "\n")
    tryCatch(
      .lat_install_github_tarball(LATELIER_REPO, subdir = p, ref = ref,
                                  upgrade = upgrade),
      error = function(e) {
        .lat_cat("    FAILED: ", p, "\n    ", conditionMessage(e), "\n")
      }
    )
  }
  .lat_cat("\n")
}

# --- Step 4: optional dependencies ----------------------------------------------
.lat_optional_dependencies <- function(pkgs) {
  groups <- Filter(function(g) any(g$pkgs %in% pkgs), LATELIER_OPTIONAL)
  if (!length(groups)) return(invisible(NULL))

  sel <- getOption("latelier.optional")
  if (is.null(sel)) {
    if (!interactive()) sel <- "none"
  }

  if (is.null(sel)) {
    .lat_cat("[L-ATELIER] Optional features\n")
    .lat_cat("    These are only needed for specific functions. You can skip them now;\n")
    .lat_cat("    functions will tell you if a package is missing when you call them.\n\n")
    sel <- character(0)
    for (nm in names(groups)) {
      g <- groups[[nm]]
      if (.lat_yesno(paste0("    Install '", nm, "' (", paste(g$pkgs, collapse = "/"),
                            ") - ", g$desc, "?"))) {
        sel <- c(sel, nm)
      }
    }
    .lat_cat("\n")
  } else if (identical(sel, "all")) {
    sel <- names(groups)
  } else if (identical(sel, "none")) {
    sel <- character(0)
  }

  bad <- setdiff(sel, names(groups))
  if (length(bad)) {
    .lat_cat("    Unknown optional group(s) ignored: ", paste(bad, collapse = ", "), "\n")
    sel <- setdiff(sel, bad)
  }
  if (!length(sel)) {
    .lat_cat("[L-ATELIER] No optional dependencies selected\n\n")
    return(invisible(NULL))
  }

  .lat_cat("[L-ATELIER] Installing optional dependencies\n")
  for (nm in sel) {
    g <- groups[[nm]]
    .lat_cat("    Group: ", nm, "\n")
    tryCatch({
      if (length(g$cran)) install.packages(setdiff(g$cran, rownames(installed.packages())))
      if (length(g$bioc)) BiocManager::install(g$bioc, ask = FALSE, update = FALSE)
      if (length(g$github)) {
        for (gh in g$github) {
          parts  <- strsplit(gh, "/", fixed = TRUE)[[1]]   # owner/repo[/subdir]
          subdir <- if (length(parts) > 2) paste(parts[-(1:2)], collapse = "/") else NULL
          .lat_install_github_tarball(paste(parts[1:2], collapse = "/"),
                                      subdir = subdir, ref = "HEAD",
                                      upgrade = "never")
        }
      }
    }, error = function(e) {
      .lat_cat("    FAILED group ", nm, ": ", conditionMessage(e), "\n")
    })
  }
  .lat_cat("\n")
}

# --- Step 5: external (non-R) requirements --------------------------------------
.lat_external_checks <- function(pkgs) {
  if ("HORIZON" %in% pkgs) {
    .lat_cat("[L-ATELIER] External tools (HORIZON)\n")
    has_conda <- nzchar(Sys.which("conda")) || nzchar(Sys.which("mamba")) ||
      dir.exists(path.expand("~/miniconda3")) || dir.exists(path.expand("~/anaconda3"))
    if (has_conda) {
      .lat_cat("    conda detected. The ATAC/ChIP/CUT&TAG pipeline needs a 'horizon_cli'\n")
      .lat_cat("    environment (samtools, bedtools, macs3, bowtie2, deeptools).\n")
    } else {
      .lat_cat("    conda was NOT detected. The bulk RNA-seq functions work without it, but\n")
      .lat_cat("    the chromatin pipeline needs conda (https://docs.conda.io/).\n")
    }
    .lat_cat("    Setup instructions: HORIZON/README.md ('Create the HORIZON conda environment')\n\n")
  }
  if ("GAIA" %in% pkgs) {
    .lat_cat("[L-ATELIER] Software that must be obtained separately (GAIA)\n")
    .lat_cat("    CIBERSORT source code must be requested from its authors (see GAIA/README.md).\n\n")
  }
}

# --- Step 6: verify -------------------------------------------------------------
.lat_verify <- function(pkgs) {
  .lat_cat("[L-ATELIER] Verifying installation\n")
  ok <- vapply(pkgs, function(p) {
    res <- tryCatch({ loadNamespace(p); TRUE }, error = function(e) FALSE)
    .lat_cat("    ", p, ": ", if (res) paste("OK (version", as.character(packageVersion(p)), ")") else "FAILED to load", "\n")
    res
  }, logical(1))

  .lat_cat("\n[L-ATELIER] Checking Bioconductor consistency (BiocManager::valid)\n")
  v <- tryCatch(BiocManager::valid(), error = function(e) NULL)
  if (isTRUE(v)) {
    .lat_cat("    All Bioconductor packages are consistent with release ",
             as.character(BiocManager::version()), "\n")
  } else if (!is.null(v)) {
    .lat_cat("    Some packages are out of date or too new for this Bioconductor release.\n")
    .lat_cat("    See the BiocManager::valid() output for details.\n")
  }

  if (all(ok)) {
    .lat_cat("\n[L-ATELIER] Done. Load with: library(", pkgs[1], ")\n")
  } else {
    .lat_cat("\n[L-ATELIER] Some packages failed. Scroll up for the first error;\n")
    .lat_cat("    missing system libraries are the most common cause (see INSTALL.md).\n")
  }
  invisible(ok)
}

# --- Main -----------------------------------------------------------------------
LATELIER_install <- function() {
  ref     <- getOption("latelier.ref", "master")
  upgrade <- getOption("latelier.upgrade", "never")

  .lat_cat("[L-ATELIER] Installer\n\n")
  .lat_check_environment()
  pkgs <- .lat_choose_packages()
  .lat_install_packages(pkgs, ref, upgrade)
  .lat_optional_dependencies(pkgs)
  .lat_external_checks(pkgs)
  .lat_verify(pkgs)
}

LATELIER_install()
