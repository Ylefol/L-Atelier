# ==============================================================================
# celltypist_markers.py — CellTypist classifier-driving genes
# Called from PANDORA_celltypist_markers() via reticulate::py_run_file()
# inside a basilisk-managed environment.
#
# Unlike celltypist_run.py (which classifies cells against a model),
# this inspects an already-trained model's linear classifier directly.
# extract_top_markers() returns the genes with the largest coefficients for
# a given class -- the genes CellTypist's model actually relied on to define
# that cell type in its reference atlas, not a statistical DE test run on
# the current dataset.
#
# Variables injected by R before this script runs (via reticulate::py_set_attr):
#   r_model_name    - model filename, e.g. "Mouse_Whole_Brain.pkl"
#   r_cell_types    - list of cell type labels (must match model.cell_types)
#   r_top_n         - int, markers per cell type
#   r_only_positive - bool, positive-coefficient markers only
#   r_force_update  - bool, re-download model even if cached
#
# Sets in Python global namespace on completion:
#   result - dict, cell type label -> list of marker genes
# ==============================================================================

import celltypist.models as ct_models

ct_models.download_models(force_update=r_force_update, model=r_model_name)
model = ct_models.Model.load(r_model_name)

missing = [ct for ct in r_cell_types if ct not in model.cell_types]
if missing:
    raise ValueError(
        "Cell type(s) not found in model '%s': %s" % (r_model_name, missing)
    )

result = {
    ct: list(model.extract_top_markers(ct, top_n=r_top_n, only_positive=r_only_positive))
    for ct in r_cell_types
}
