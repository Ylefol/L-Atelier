# ==============================================================================
# CAULDRON - basilisk Python environment definitions
# ==============================================================================
# All Python environments used by CAULDRON are declared here at the package
# level. basilisk reads these at load time to manage isolated conda
# environments — users never interact with conda or pip directly.
#
# Naming convention: .cauldron_<tool>_env_<version_suffix>
# Increment the version suffix when package pins change so that basilisk
# creates a fresh environment rather than trying to update an existing one.
# ==============================================================================

# ── CellTypist ─────────────────────────────────────────────────────────────────
# Versions from a validated working environment (Python 3.11, celltypist 1.6.3).
# celltypist pulls leidenalg, scikit-learn, and other dependencies via pip;
# pinning the three main packages is sufficient for reproducibility.
.cauldron_celltypist_env <- basilisk::BasiliskEnvironment(
    envname  = "cauldron_celltypist_1",
    pkgname  = "CAULDRON",
    packages = "python==3.11",
    pip      = c(
        "celltypist==1.6.3",
        "scanpy==1.11.1",
        "anndata==0.11.4"
    )
)
