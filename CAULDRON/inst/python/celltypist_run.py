# ==============================================================================
# celltypist_run.py — CellTypist annotation
# Called from PANDORA_annotate_celltypist() via reticulate::py_run_file()
# inside a basilisk-managed environment.
#
# Variables injected by R before this script runs (via reticulate::py$):
#   r_mat_t           - cells × genes sparse matrix (scipy csc_matrix)
#   r_genes           - list/array of gene names
#   r_cells           - list/array of cell barcodes
#   r_model_name      - model filename, e.g. "Mouse_Whole_Brain.pkl"
#   r_majority_voting - bool: consensus label per over-cluster
#   r_force_update    - bool: re-download model even if cached
#
# Sets in Python global namespace on completion:
#   result_df - pandas DataFrame indexed by cell barcode
#               always contains "predicted_labels"
#               also contains "conf_score" when majority_voting=False
# ==============================================================================

import celltypist
import celltypist.models as ct_models
import anndata as ad
import scanpy as sc

# Download / use cached model
ct_models.download_models(force_update=r_force_update, model=r_model_name)

# Build AnnData object from the sparse matrix passed from R
# reticulate converts dgCMatrix -> scipy csc_matrix automatically;
# AnnData expects cells x genes
adata = ad.AnnData(X=r_mat_t)
adata.obs_names = list(r_cells)
adata.var_names = list(r_genes)

# Normalise: 10,000 counts per cell then log1p — required by CellTypist
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Annotate
predictions = celltypist.annotate(
    adata,
    model=r_model_name,
    majority_voting=r_majority_voting
)

# predicted_labels is a pandas DataFrame; R retrieves it via py_to_r()
result_df = predictions.predicted_labels
