# ==============================================================================
# GAIA - basilisk managed environment definitions
# ==============================================================================
# Environments used by GAIA for external CLI tool wrappers: isoform-switch
# consequence annotation (CPAT coding-potential prediction, SignalP 6
# signal-peptide prediction, Pfam protein domain annotation via
# pfam_scan.pl) and motif enrichment (HOMER). basilisk manages isolated
# conda environments so users don't interact with conda/pip directly for
# these tools -- note this isn't Python-specific: basilisk is a general
# reproducible-conda-env manager (every env includes a Python interpreter
# since that's what basilisk/reticulate need to drive it, but the env can
# carry arbitrary conda packages, including native-binary/Perl tools like
# pfam_scan.pl or homer). Naming convention matches CAULDRON/CYAN:
# .gaia_<tool>_env. Increment the envname suffix when package pins change so
# basilisk creates a fresh environment rather than trying to update an
# existing one.
# ==============================================================================

# -- CPAT (coding potential) --------------------------------------------------
# Version pin is an unverified starting point (CPAT 3.0.5 is the current
# PyPI release as of this writing, pip-installable) -- confirm it resolves
# correctly on first real run and adjust if needed.
.gaia_cpat_env <- basilisk::BasiliskEnvironment(
  envname  = "gaia_cpat_1",
  pkgname  = "GAIA",
  packages = "python==3.10",
  pip      = c("CPAT==3.0.5")
)

# -- SignalP 6 (signal peptide prediction) ------------------------------------
# SignalP 6 is license-gated (DTU Health Tech) and not available on public
# PyPI/conda channels, so it cannot be declared here the way CPAT is. This
# environment only provisions the base Python interpreter; SignalP 6 itself
# must be installed manually, once, following DTU's own installation
# instructions (including copying the separately-distributed model weights)
# using this environment's own pip -- see basilisk::obtainEnvironmentPath()
# to locate it. APOLLO_run_signalp() checks for the installed executable and
# stops with setup instructions if it's missing rather than attempting to
# install it automatically.
.gaia_signalp_env <- basilisk::BasiliskEnvironment(
  envname  = "gaia_signalp_1",
  pkgname  = "GAIA",
  packages = "python==3.10"
)

# -- Pfam (protein domains, via pfam_scan.pl) ---------------------------------
# pfam_scan.pl (the Perl wrapper around HMMER's hmmscan that
# IsoformSwitchAnalyzeR::analyzePFAM() expects output from) is packaged on
# bioconda as `pfam_scan`, which declares HMMER and its own Perl module
# dependencies (incl. Bio::Pfam::Scan) as part of its own recipe -- so this
# single conda install replaces what would otherwise be a manual
# PfamScan.tar.gz + CPAN Bio::Pfam::Scan setup. basilisk requires an explicit
# version pin on every conda package (not just python) -- 1.6 is the latest
# bioconda release as of this writing (bioconda.github.io/recipes/pfam_scan).
.gaia_pfam_env <- basilisk::BasiliskEnvironment(
  envname  = "gaia_pfam_1",
  pkgname  = "GAIA",
  packages = c("python==3.10", "pfam_scan==1.6"),
  channels = c("bioconda", "conda-forge", "defaults")
)

# -- HOMER (motif enrichment) --------------------------------------------------
# findMotifsGenome.pl is packaged on bioconda as `homer`. Note: HOMER also
# requires per-genome reference data installed separately via its own
# `configureHomer.pl -install <genome>` (downloads genome FASTA/annotation
# into HOMER's data directory) -- that data-install step is independent of
# how findMotifsGenome.pl itself is obtained and still needs to be run once
# against this environment. basilisk requires an explicit version pin on
# every conda package (not just python) -- 5.1 is the latest bioconda release
# as of this writing (bioconda.github.io/recipes/homer).
.gaia_homer_env <- basilisk::BasiliskEnvironment(
  envname  = "gaia_homer_1",
  pkgname  = "GAIA",
  packages = c("python==3.10", "homer==5.1"),
  channels = c("bioconda", "conda-forge", "defaults")
)
