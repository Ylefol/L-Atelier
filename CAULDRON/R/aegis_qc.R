# ==============================================================================
# AEGIS - Quality Control & Filtering
# CAULDRON: Single-Cell Analysis Toolkit
# ==============================================================================
# The divine shield forged by Hephaestus. AEGIS protects the dataset from
# bad data — filtering low-quality cells, doublets, and ambient RNA
# contamination before any analysis proceeds.
#
# Core responsibilities:
#   - Per-cell QC metric computation (nUMI, nGenes, MT%, ribo%)
#   - Adaptive threshold filtering (MAD-based)
#   - Doublet detection (scDblFinder)
#   - Ambient RNA removal (SoupX / DecontX)
#   - QC visualisation (violin plots, scatter plots, knee plots)
#   - Sample-level QC summaries
#
# All functions prefixed: AEGIS_
# ==============================================================================
