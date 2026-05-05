"""
parse_dge_filter.py
-------------------
Filters a raw velocity AnnData (from parse_velocity.py) to the cells and
genes present in a PARSE split-pipe DGE directory, then adds cell and gene
metadata from cell_metadata.csv and all_genes.csv.

Called from HORIZON_parse_DGE_filter() via system2().

Usage:
    python parse_dge_filter.py \
        --velocity_h5ad <path to raw adata_vel.h5ad> \
        --cell_metadata <path to DGE_filtered/cell_metadata.csv> \
        --all_genes     <path to DGE_filtered/all_genes.csv> \
        --output_file   <path for filtered output h5ad>
"""

import argparse
import os
import sys
import pandas as pd
import anndata as ad

# ---------------------------------------------------------------------------
# Arguments
# ---------------------------------------------------------------------------

parser = argparse.ArgumentParser(description="Filter velocity AnnData using PARSE DGE output")
parser.add_argument("--velocity_h5ad", required=True,
                    help="Path to the raw velocity AnnData from parse_velocity.py")
parser.add_argument("--cell_metadata", required=True,
                    help="Path to cell_metadata.csv from the DGE directory")
parser.add_argument("--all_genes",     required=True,
                    help="Path to all_genes.csv from the DGE directory")
parser.add_argument("--output_file",   required=True,
                    help="Output path for the filtered AnnData")
args = parser.parse_args()

# ---------------------------------------------------------------------------
# Load AnnData
# ---------------------------------------------------------------------------

print(f"[parse_dge_filter] Loading velocity AnnData: {args.velocity_h5ad}")
adata = ad.read_h5ad(args.velocity_h5ad)
print(f"    Cells: {adata.n_obs}  |  Genes: {adata.n_vars}")

# ---------------------------------------------------------------------------
# Load and prepare cell metadata
# ---------------------------------------------------------------------------

print(f"[parse_dge_filter] Loading cell metadata: {args.cell_metadata}")
cell_meta = pd.read_csv(args.cell_metadata)

if 'bc_wells' not in cell_meta.columns:
    print("ERROR: 'bc_wells' column not found in cell_metadata.csv")
    print(f"  Available columns: {list(cell_meta.columns)}")
    sys.exit(1)

cell_meta = cell_meta.set_index('bc_wells')
print(f"    Cells in DGE metadata: {len(cell_meta)}")

# ---------------------------------------------------------------------------
# Load and prepare gene metadata
# ---------------------------------------------------------------------------

print(f"[parse_dge_filter] Loading gene metadata: {args.all_genes}")
gene_meta = pd.read_csv(args.all_genes)

if 'gene_name' not in gene_meta.columns:
    print("ERROR: 'gene_name' column not found in all_genes.csv")
    print(f"  Available columns: {list(gene_meta.columns)}")
    sys.exit(1)

# Strip any stray quotes from gene names (tscp CSV parsing artefact)
gene_meta['gene_name'] = gene_meta['gene_name'].str.strip('"')
gene_meta = gene_meta.set_index('gene_name')
n_before_dedup = len(gene_meta)
gene_meta = gene_meta[~gene_meta.index.duplicated(keep='first')]
if len(gene_meta) < n_before_dedup:
    print(f"    Dropped {n_before_dedup - len(gene_meta)} duplicate gene name(s)")
print(f"    Genes in DGE metadata: {len(gene_meta)}")

# ---------------------------------------------------------------------------
# Filter cells
# ---------------------------------------------------------------------------

# obs_names after concat carry the split-pipe run suffix (__s1, __s2 …) and
# match bc_wells in the combined cell_metadata.csv directly
n_before  = adata.n_obs
cell_mask = adata.obs_names.isin(cell_meta.index)
adata     = adata[cell_mask].copy()
print(f"[parse_dge_filter] Cell filter: kept {adata.n_obs} / {n_before}")

if adata.n_obs == 0:
    print("ERROR: No cells remained after filtering.")
    print("  Check that obs_names in the velocity h5ad match bc_wells in cell_metadata.csv.")
    print(f"  Example obs_names : {list(adata.obs_names[:3])}")
    print(f"  Example bc_wells  : {list(cell_meta.index[:3])}")
    sys.exit(1)

# ---------------------------------------------------------------------------
# Filter genes
# ---------------------------------------------------------------------------

n_before  = adata.n_vars
gene_mask = adata.var_names.isin(gene_meta.index)
adata     = adata[:, gene_mask].copy()
print(f"[parse_dge_filter] Gene filter:  kept {adata.n_vars} / {n_before}")

# ---------------------------------------------------------------------------
# Add metadata
# ---------------------------------------------------------------------------

# Cell metadata — join on obs_names (bc_wells with run suffix)
adata.obs = adata.obs.join(cell_meta, how='left')
print(f"[parse_dge_filter] Added {len(cell_meta.columns)} cell metadata columns")

# Gene metadata — join on var_names (gene_name)
adata.var = adata.var.join(gene_meta, how='left')
print(f"[parse_dge_filter] Added {len(gene_meta.columns)} gene metadata columns")

# ---------------------------------------------------------------------------
# Write output
# ---------------------------------------------------------------------------

os.makedirs(os.path.dirname(os.path.abspath(args.output_file)), exist_ok=True)
adata.write(args.output_file)
print(f"[parse_dge_filter] Filtered AnnData written to: {args.output_file}")
print(f"    Final cells: {adata.n_obs}  |  Final genes: {adata.n_vars}")
