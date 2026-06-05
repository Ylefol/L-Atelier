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

# ── scVelo ─────────────────────────────────────────────────────────────────────
# Versions matched to a validated working environment (user's scvelo conda env).
# scVelo 0.3.3 rewrote the stochastic mode internals, resolving the numpy >= 1.24
# leastsq_generalized bug present in 0.2.x. All three modes (deterministic,
# stochastic, dynamical) work with this stack.
# legacy-api-wrap is required by scVelo 0.3.x for backward-compatibility shims.
.cauldron_scvelo_env <- basilisk::BasiliskEnvironment(
    envname  = "cauldron_scvelo_3",
    pkgname  = "CAULDRON",
    packages = "python==3.12",
    pip      = c(
        "scvelo==0.3.3",
        "scanpy==1.11.1",
        "anndata==0.11.4",
        "numpy==2.2.5",
        "scipy==1.15.2",
        "matplotlib==3.10.1",
        "pandas==2.2.3",
        "scikit-learn==1.5.2",
        "legacy-api-wrap==1.4.1"   # required by scVelo 0.3.x for backward-compat shims
    )
)

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

# ── scVI ───────────────────────────────────────────────────────────────────────
# scvi-tools requires PyTorch (~1 GB); first-time setup installs the full stack.
# Versions pinned against a validated working combination for Python 3.11 /
# torch 2.2. Increment the envname suffix when pins change.
.cauldron_scvi_env <- basilisk::BasiliskEnvironment(
    envname  = "cauldron_scvi_1",
    pkgname  = "CAULDRON",
    packages = c("python==3.11", "libexpat==2.6.1"),
    pip      = c(
        "setuptools==69.5.1",  # provides pkg_resources, not bundled in basilisk envs
        "scvi-tools==1.2.0",
        "torch==2.2.2",
        "anndata==0.10.8",
        "scanpy==1.10.1",
        "numpy==1.26.4",
        "scipy==1.13.1",
        "pandas==2.2.2",
        "lightning==2.2.5"
    )
)
