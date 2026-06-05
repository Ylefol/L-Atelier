# ==============================================================================
# scvi_run.py — scVI latent embedding
# Called from TALOS_run_scvi() via reticulate::py_run_file()
# inside a basilisk-managed environment (Python 3.11, scvi-tools 1.2.0).
#
# Variables injected by R before this script runs (via py_set_attr):
#   r_counts     - genes × cells scipy CSC sparse matrix (raw integer counts)
#   r_genes      - list of gene names
#   r_cells      - list of cell barcodes
#   r_batch      - list of batch label strings (one per cell) OR None
#   r_n_latent   - int: number of latent dimensions (default 10)
#   r_n_layers   - int: number of encoder/decoder layers (default 2)
#   r_n_hidden   - int: nodes per hidden layer (default 128)
#   r_max_epochs - int: maximum training epochs (default 400)
#   r_seed       - int: random seed
#
# Sets in Python global namespace on completion:
#   result_latent - numpy array (cells × n_latent), float64
# ==============================================================================

import scvi
import anndata as ad
import numpy as np
import scipy.sparse as sp
import warnings

# Suppress lightning/scVI training progress output inside basilisk
import logging
logging.getLogger("lightning.pytorch").setLevel(logging.ERROR)
warnings.filterwarnings("ignore", category=UserWarning)

scvi.settings.seed = int(r_seed)
scvi.settings.verbosity = 0

# ------------------------------------------------------------------------------
# 1. Build AnnData (cells × genes — AnnData convention)
# r_counts arrives as genes × cells scipy CSC from R's dgCMatrix.
# Transpose to cells × genes and force writable float32 CSR.
# scVI expects non-negative integer (or float) counts; float32 is memory-
# efficient and matches scVI's internal dtype.
# ------------------------------------------------------------------------------
def _writable_csr(mat):
    m = sp.csr_matrix(mat, dtype=np.float32)
    m.data    = m.data.copy()
    m.indices = m.indices.copy()
    m.indptr  = m.indptr.copy()
    return m

counts_csr = _writable_csr(r_counts.T)   # cells × genes

adata = ad.AnnData(X=counts_csr)
adata.obs_names = list(r_cells)
adata.var_names = list(r_genes)

# ------------------------------------------------------------------------------
# 2. Attach batch column if provided
# ------------------------------------------------------------------------------
if r_batch is not None:
    adata.obs["batch"] = list(r_batch)

# ------------------------------------------------------------------------------
# 3. Setup scVI and train
# setup_anndata registers the data layout (which obs column is batch, etc.).
# SCVI() defines the model architecture; train() fits the variational posterior.
# early_stopping halts training when the ELBO on a held-out validation set
# stops improving, which usually fires well before max_epochs.
# ------------------------------------------------------------------------------
scvi.model.SCVI.setup_anndata(
    adata,
    batch_key = "batch" if r_batch is not None else None
)

model = scvi.model.SCVI(
    adata,
    n_latent = int(r_n_latent),
    n_layers = int(r_n_layers),
    n_hidden = int(r_n_hidden)
)

model.train(
    max_epochs          = int(r_max_epochs),
    early_stopping      = True,
    enable_progress_bar = False
)

# ------------------------------------------------------------------------------
# 4. Extract latent representation
# get_latent_representation() returns a cells × n_latent numpy float32 array.
# Cast to float64 for consistency with R's double precision.
# ------------------------------------------------------------------------------
result_latent = np.asarray(
    model.get_latent_representation(), dtype=np.float64
)
