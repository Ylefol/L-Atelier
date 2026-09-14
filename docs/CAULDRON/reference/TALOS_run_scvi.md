# scVI latent embedding

Trains a scVI (single-cell Variational Inference) model on the raw count
matrix and extracts the per-cell latent representation. The result is
stored in `reducedDims(sce)[["scVI"]]` and can be used as a drop-in
replacement for PCA in
[`TALOS_run_umap`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TALOS_run_umap.md)
and
[`TALOS_build_graph`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TALOS_build_graph.md)
via their `use_rep` argument.

## Usage

``` r
TALOS_run_scvi(
  sce,
  n_latent = 10L,
  n_layers = 2L,
  n_hidden = 128L,
  max_epochs = 400L,
  batch_key = NULL,
  seed = 42L,
  verbose = TRUE
)
```

## Arguments

- sce:

  A `SingleCellExperiment` with a `"counts"` assay containing raw
  integer counts.

- n_latent:

  Integer. Dimensionality of the latent space. Default `10`. Higher
  values capture more variance but may include noise; 10–20 is a typical
  range for scRNA-seq.

- n_layers:

  Integer. Number of encoder/decoder hidden layers. Default `2`.

- n_hidden:

  Integer. Number of nodes per hidden layer. Default `128`.

- max_epochs:

  Integer. Maximum training epochs. Early stopping is always enabled and
  will usually halt well before this limit. Default `400`.

- batch_key:

  Character or `NULL`. Column name in `colData(sce)` containing batch
  labels. When provided, the scVI model conditions on batch during
  training, producing a batch-corrected latent space. `NULL` (default)
  trains without batch correction.

- seed:

  Integer. Random seed for reproducibility. Default `42L`.

- verbose:

  Logical. Print a summary. Default `TRUE`.

## Value

SCE with:

- `reducedDims(sce)[["scVI"]]`:

  Cells × `n_latent` matrix of latent coordinates.

- `metadata(sce)$scvi_params`:

  List recording the training parameters for reproducibility.

## Details

scVI models counts with a negative binomial distribution and learns a
low-dimensional latent space that accounts for library size, technical
noise, and — when `batch_key` is provided — batch effects. Unlike PCA,
no pre-normalisation is required: the function reads the raw `"counts"`
assay directly.

Training is performed inside an isolated basilisk Python environment
(Python 3.11, scvi-tools 1.2.0, PyTorch 2.2). The environment is created
automatically on first use (~1 GB download).
