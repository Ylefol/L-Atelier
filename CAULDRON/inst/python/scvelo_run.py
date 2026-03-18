# ==============================================================================
# scvelo_run.py — RNA velocity estimation via scVelo
# Called from TRIPODES_run_velocity() via reticulate::py_run_file()
# inside a basilisk-managed environment (Python 3.10, numpy==1.23.5).
#
# Variables injected by R before this script runs (via py_set_attr):
#   r_spliced     - genes × cells sparse matrix (scipy CSC from dgCMatrix)
#   r_unspliced   - genes × cells sparse matrix (scipy CSC from dgCMatrix)
#   r_genes       - list/array of gene names
#   r_cells       - list/array of cell barcodes
#   r_mode        - str: "deterministic", "stochastic", or "dynamical"
#   r_n_pcs       - int: PCs for neighbour graph construction
#   r_n_neighbors - int: k for kNN graph
#   r_pca         - numpy array (cells × n_pcs) OR None
#                   Pre-computed PCA from TALOS; ensures velocity neighbours
#                   are built in the same space as the UMAP embedding.
#   r_umap        - numpy array (cells × 2) OR None
#   r_tsne        - numpy array (cells × 2) OR None
#
# Sets in Python global namespace on completion:
#   result_velocity      - numpy array (genes × cells), float64
#   result_confidence    - numpy array (cells,), float64
#   result_length        - numpy array (cells,), float64
#   result_velocity_pca  - numpy array (cells × n_pcs) OR None
#   result_velocity_umap - numpy array (cells × 2) OR None
#   result_velocity_tsne - numpy array (cells × 2) OR None
# ==============================================================================

import scvelo as scv
import anndata as ad
import scanpy as sc
import numpy as np
import scipy.sparse as sp

# Suppress scVelo progress output inside basilisk
scv.settings.verbosity = 0

# ------------------------------------------------------------------------------
# 1. Build AnnData (cells × genes — AnnData convention)
# r_spliced / r_unspliced arrive as scipy CSC matrices, shape (genes, cells),
# backed by R memory — their numpy arrays are read-only.  scVelo does in-place
# scaling during moments computation, so we must force a writable copy of all
# three CSR component arrays before building AnnData.
# ------------------------------------------------------------------------------
def _writable_csr(mat):
    m = sp.csr_matrix(mat, dtype=np.float64)
    m.data    = m.data.copy()
    m.indices = m.indices.copy()
    m.indptr  = m.indptr.copy()
    return m

spliced_cx   = _writable_csr(r_spliced.T)
unspliced_cx = _writable_csr(r_unspliced.T)

adata = ad.AnnData(X=spliced_cx.copy())
adata.obs_names = list(r_cells)
adata.var_names = list(r_genes)

adata.layers["spliced"]   = spliced_cx
adata.layers["unspliced"] = unspliced_cx

# ------------------------------------------------------------------------------
# 2. Inject pre-computed PCA if provided
# scVelo's pp.moments will use adata.obsm["X_pca"] for neighbour construction,
# ensuring velocity neighbours are consistent with the UMAP embedding.
# ------------------------------------------------------------------------------
if r_pca is not None:
    adata.obsm["X_pca"] = np.asarray(r_pca, dtype=np.float64)

# ------------------------------------------------------------------------------
# 3. Normalise spliced layer
# Required for the kNN graph underlying moments computation.
# Unspliced is not normalised here — scVelo uses raw unspliced for velocity fit.
# ------------------------------------------------------------------------------
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# ------------------------------------------------------------------------------
# 4. Build kNN graph then compute moments
# sc.pp.neighbors is called explicitly (scVelo 0.3.x deprecates implicit graph
# construction inside moments).  If PCA was provided, use it as the embedding
# so the velocity graph is consistent with the UMAP/TALOS neighbourhood.
# ------------------------------------------------------------------------------
if r_pca is not None:
    sc.pp.neighbors(adata, use_rep="X_pca",
                    n_pcs=int(r_n_pcs), n_neighbors=int(r_n_neighbors))
else:
    sc.pp.neighbors(adata, n_pcs=int(r_n_pcs), n_neighbors=int(r_n_neighbors))

scv.pp.moments(adata, n_pcs=None, n_neighbors=None)

# ------------------------------------------------------------------------------
# 5. Velocity estimation
# ------------------------------------------------------------------------------
scv.tl.velocity(adata, mode=str(r_mode))

# ------------------------------------------------------------------------------
# 6. Velocity graph (transition probability matrix between cells)
# ------------------------------------------------------------------------------
scv.tl.velocity_graph(adata)

# ------------------------------------------------------------------------------
# 7. Per-cell confidence and length scores
# ------------------------------------------------------------------------------
scv.tl.velocity_confidence(adata)

# ------------------------------------------------------------------------------
# 8. Project velocity onto all provided embeddings
# velocity_embedding is cheap — the expensive steps are already done above.
# ------------------------------------------------------------------------------
if r_pca is not None:
    scv.tl.velocity_embedding(adata, basis="pca")

if r_umap is not None:
    adata.obsm["X_umap"] = np.asarray(r_umap, dtype=np.float64)
    scv.tl.velocity_embedding(adata, basis="umap")

if r_tsne is not None:
    adata.obsm["X_tsne"] = np.asarray(r_tsne, dtype=np.float64)
    scv.tl.velocity_embedding(adata, basis="tsne")

# ------------------------------------------------------------------------------
# 9. Extract results — all returned to R as dense numpy arrays
# Velocity is transposed back to genes × cells (R convention).
# ------------------------------------------------------------------------------
vel_layer = adata.layers["velocity"]
result_velocity = np.asarray(
    vel_layer.todense() if sp.issparse(vel_layer) else vel_layer,
    dtype=np.float64
).T   # shape: (genes, cells)

result_confidence = np.asarray(
    adata.obs["velocity_confidence"].values, dtype=np.float64
)

result_length = np.asarray(
    adata.obs["velocity_length"].values, dtype=np.float64
)

result_velocity_pca = (
    np.asarray(adata.obsm["velocity_pca"], dtype=np.float64)
    if (r_pca is not None and "velocity_pca" in adata.obsm)
    else None
)

result_velocity_umap = (
    np.asarray(adata.obsm["velocity_umap"], dtype=np.float64)
    if (r_umap is not None and "velocity_umap" in adata.obsm)
    else None
)

result_velocity_tsne = (
    np.asarray(adata.obsm["velocity_tsne"], dtype=np.float64)
    if (r_tsne is not None and "velocity_tsne" in adata.obsm)
    else None
)
