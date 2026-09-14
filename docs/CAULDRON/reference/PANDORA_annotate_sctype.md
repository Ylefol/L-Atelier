# Marker-based cell type annotation via scType

Scores clusters against a curated marker gene database to assign cell
type labels. Annotation is cluster-level: expression is aggregated per
cluster, scored against each cell type's positive and negative markers
(from the bundled `ScTypeDB_full.xlsx`), and the highest-scoring type is
assigned to all cells in each cluster.

## Usage

``` r
PANDORA_annotate_sctype(
  sce,
  tissue,
  cluster_col = "cluster",
  db = NULL,
  species = "human",
  ortholog_method = "homologene",
  assay_name = "logcounts",
  label_col = "sctype_label",
  verbose = TRUE
)
```

## Arguments

- sce:

  A `SingleCellExperiment` with a `"logcounts"` (or other) assay.

- tissue:

  Character. Tissue type string matching an entry in the scTypeDB
  (case-sensitive).

- cluster_col:

  Character. `colData` column with cluster assignments. Default
  `"cluster"`.

- db:

  `NULL` (bundled `ScTypeDB_full.xlsx`), a file path to a custom Excel
  or CSV, or a `data.frame` with columns `tissueType`, `cellName`,
  `geneSymbolmore1` (positive markers, comma-separated),
  `geneSymbolmore2` (negative markers, comma-separated; `NA` = none).

- species:

  Character. Target organism. `"human"` (default) uses the DB markers
  directly. Any other value (e.g. `"mouse"`, `"Mus musculus"`) triggers
  ortholog conversion via `orthogene`. Accepts any species string
  recognised by `orthogene`.

- ortholog_method:

  Character. Ortholog database passed to
  [`orthogene::convert_orthologs()`](https://rdrr.io/pkg/orthogene/man/convert_orthologs.html).
  Default `"homologene"` (offline, fast, covers common model organisms).
  Use `"gprofiler"` for broader species coverage (requires internet).

- assay_name:

  Character. Assay to score. Default `"logcounts"`.

- label_col:

  Character. `colData` column for the result. Default `"sctype_label"`.

- verbose:

  Logical. Print a summary. Default `TRUE`.

## Value

SCE with `colData(sce)[[label_col]]` populated.

## Details

The bundled ScTypeDB uses human HGNC gene symbols. For non-human data,
set `species` to the target organism and orthologous gene symbols will
be resolved via `orthogene` before scoring. The conversion report
(database used, mapping rate) is printed to allow you to judge the
reliability of the annotation before interpreting results.

Use
[`PANDORA_list_references`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/PANDORA_list_references.md)`(source = "sctype")`
to see available tissues and cell type counts before running.
