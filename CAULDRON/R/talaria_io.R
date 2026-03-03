# ==============================================================================
# TALARIA - Import, Export & Format Conversion
# CAULDRON: Single-Cell Analysis Toolkit
# ==============================================================================
# The winged sandals forged by Hephaestus for Hermes. TALARIA carries data
# in and out of the toolkit — pure transport, bridging external formats and
# the internal SingleCellExperiment object model.
#
# Core responsibilities:
#   - Load 10x Genomics MEX directories (barcodes/features/matrix)
#   - Load H5 and H5AD (AnnData) files
#   - Load Loom files
#   - Convert between SingleCellExperiment and Seurat objects
#   - Export results (count matrices, metadata, reduced dims)
#   - Report generation
#
# All functions prefixed: TALARIA_
# ==============================================================================
