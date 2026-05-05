"""
parse_velocity.py
-----------------
Generates spliced/unspliced AnnData objects from PARSE Biosciences
tscp_assignment.csv files for RNA velocity analysis.

Called from HORIZON_run_parse_velocity() via system2(). Arguments are
passed as command-line flags by the R wrapper.

Produces a raw concatenated AnnData containing all cells across all runs
with spliced and unspliced layers. No cell/gene filtering or metadata
addition is performed here — use parse_dge_filter.py for that step.

Usage:
    python parse_velocity.py \
        --working_dir  <root output dir from split-pipe> \
        --run_ids      <comma-separated run IDs> \
        --output_file  <path for final adata_vel.h5ad>

Note: tscp_assignment.csv files must be uncompressed before running.
      HORIZON_run_parse_velocity() warns if .gz files are found.
"""

import argparse
import gc
import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse
import scvelo as scv
import anndata as ad
import scipy.io as sio
import dask.dataframe as dd
from dask.diagnostics import ProgressBar

ProgressBar().register()

# ---------------------------------------------------------------------------
# Arguments
# ---------------------------------------------------------------------------

parser = argparse.ArgumentParser(description="PARSE velocity matrix generation")
parser.add_argument("--working_dir",   required=True,
                    help="Root output directory containing per-run split-pipe results")
parser.add_argument("--run_ids",       required=True,
                    help="Comma-separated run IDs (must match subdirectory names)")
parser.add_argument("--output_file",   required=True,
                    help="Output path for the final adata_vel.h5ad")
args = parser.parse_args()

working_dir = args.working_dir
run_ids     = [r.strip() for r in args.run_ids.split(",")]
output_file = args.output_file

# Resolve tscp paths from explicit run IDs
tscp_paths = [
    os.path.join(working_dir, rid, "process", "tscp_assignment.csv")
    for rid in run_ids
]

missing = [p for p in tscp_paths if not os.path.exists(p)]
if missing:
    print("ERROR: tscp_assignment.csv not found for the following runs:")
    for p in missing:
        print(f"  {p}")
    print("Ensure split-pipe has completed and files are uncompressed (gunzip).")
    sys.exit(1)

print(f"[parse_velocity] Processing {len(tscp_paths)} run(s)")
for p in tscp_paths:
    print(f"  {p}")

# ---------------------------------------------------------------------------
# Per-run spliced/unspliced matrix generation
# ---------------------------------------------------------------------------

def generate_splice_matrices(tscp_path, run_idx):
    print(f"\n[parse_velocity] Reading {tscp_path}")
    tscp_assign_df = dd.read_csv(tscp_path, blocksize="800MB")
    tscp_assign_df = tscp_assign_df.compute()
    tscp_assign_df['gene_name'] = tscp_assign_df['gene_name'].str.strip('"')

    cell_tscp_cnts = tscp_assign_df.groupby("bc_wells").size()
    filtered_cell_dict = dict(zip(cell_tscp_cnts.index,
                                  np.zeros(len(cell_tscp_cnts))))

    def check_filtered_cell(cell_ind):
        try:
            filtered_cell_dict[cell_ind]
        except KeyError:
            return False
        else:
            return True

    genes     = tscp_assign_df.gene_name.unique()
    bcs       = cell_tscp_cnts.index
    gene_dict = dict(zip(genes, range(len(genes))))
    barcode_dict = dict(zip(bcs, range(len(bcs))))

    reads_to_keep = tscp_assign_df.bc_wells.apply(check_filtered_cell)

    print("[parse_velocity] Filtering tscp file...")
    tscp_assign_df_filt = tscp_assign_df[reads_to_keep]
    tscp_assign_df_filt["cell_index"] = tscp_assign_df_filt.bc_wells.apply(
        lambda s: barcode_dict[s])
    tscp_assign_df_filt["gene_index"] = tscp_assign_df_filt.gene_name.apply(
        lambda s: gene_dict[s])
    print(f"    Done: {tscp_assign_df_filt.shape}")

    # Spliced (exonic)
    rcv  = tscp_assign_df_filt.query("exonic").groupby(
        ["cell_index", "gene_index"]).size().reset_index().values
    rows = list(rcv[:, 0]) + [len(barcode_dict) - 1]
    cols = list(rcv[:, 1]) + [len(genes) - 1]
    vals = list(rcv[:, 2]) + [0]
    X_exonic = scipy.sparse.csr_matrix((vals, (rows, cols)))

    # Unspliced (intronic)
    rcv  = tscp_assign_df_filt.query("~exonic").groupby(
        ["cell_index", "gene_index"]).size().reset_index().values
    rows = list(rcv[:, 0]) + [len(barcode_dict) - 1]
    cols = list(rcv[:, 1]) + [len(genes) - 1]
    vals = list(rcv[:, 2]) + [0]
    X_intronic = scipy.sparse.csr_matrix((vals, (rows, cols)))

    X     = X_exonic + X_intronic
    adata = scv.AnnData(X=X)

    adata.obs = pd.DataFrame({"barcodes": bcs}, index=bcs)
    adata.var = pd.DataFrame({"gene": genes, "gene_name": genes})
    adata.var.index = genes

    adata.var_names_make_unique()
    adata.obs_names_make_unique()
    adata.layers["spliced"]   = X_exonic
    adata.layers["unspliced"] = X_intronic
    scv.utils.show_proportions(adata)

    adata.obs.index = adata.obs.index.astype(str)
    adata.var.index = adata.var.index.astype(str)
    for col in adata.obs.select_dtypes(include="string").columns:
        adata.obs[col] = adata.obs[col].astype(str)
    for col in adata.var.select_dtypes(include="string").columns:
        adata.var[col] = adata.var[col].astype(str)

    run_out = os.path.join(working_dir, f"sub{run_idx + 1}_adata.h5ad")
    adata.write(run_out)
    print(f"[parse_velocity] Saved per-run AnnData: {run_out}")
    return adata


# ---------------------------------------------------------------------------
# Run and concatenate
# ---------------------------------------------------------------------------

ad_list_sp = []
for i, tscp_path in enumerate(tscp_paths):
    ad_list_sp.append(generate_splice_matrices(tscp_path, i))

print("\n[parse_velocity] Concatenating runs...")
subs      = list(range(1, len(ad_list_sp) + 1))
ad_splice = ad.concat(ad_list_sp, keys=subs, index_unique="__s")
print(f"[parse_velocity] Total cells after concatenation: {ad_splice.n_obs}")

os.makedirs(os.path.dirname(os.path.abspath(output_file)), exist_ok=True)
ad_splice.write(output_file)
print(f"[parse_velocity] Raw velocity AnnData written to: {output_file}")
print(f"    Cells: {ad_splice.n_obs}  |  Genes: {ad_splice.n_vars}")
print("    Run HORIZON_parse_DGE_filter() to filter and add metadata.")
