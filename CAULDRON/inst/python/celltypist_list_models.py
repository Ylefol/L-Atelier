# ==============================================================================
# celltypist_list_models.py — CellTypist model registry fetch
# Called from PANDORA_list_celltypist_models() via reticulate::py_run_file()
# inside a basilisk-managed environment.
#
# Requires internet access (model index is hosted by the Sanger Institute).
#
# Sets in Python global namespace on completion:
#   result_df - pandas DataFrame of available model descriptions
# ==============================================================================

import celltypist.models as ct_models

result_df = ct_models.models_description()
