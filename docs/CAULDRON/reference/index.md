# Package index

## Talaria

Winged sandals of Hermes, pure transport — I/O: load 10x
MEX/H5/H5AD/Loom, Seurat interop, BPCells backing, export (matrices,
DE/GSEA/ORA results, interactive R Markdown widgets).

- [`TALARIA_build_de_comparison_table()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TALARIA_build_de_comparison_table.md)
  : Build a browsable table of DE results across several comparisons
- [`TALARIA_export_de()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TALARIA_export_de.md)
  : Export DE results to files
- [`TALARIA_export_de_comparison_widget()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TALARIA_export_de_comparison_widget.md)
  : Export an interactive DE cross-comparison bar-chart widget for R
  Markdown
- [`TALARIA_export_first_steps()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TALARIA_export_first_steps.md)
  : Export a first-steps analysis
- [`TALARIA_export_gene_expr_widget()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TALARIA_export_gene_expr_widget.md)
  : Export an interactive gene-expression-by-group widget for R Markdown
- [`TALARIA_export_gene_table_widget()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TALARIA_export_gene_table_widget.md)
  : Export a searchable rowData(sce) table widget for R Markdown
- [`TALARIA_export_gsea()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TALARIA_export_gsea.md)
  : Export GSEA results to files
- [`TALARIA_export_matrix()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TALARIA_export_matrix.md)
  : Export a count matrix to disk
- [`TALARIA_export_metadata()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TALARIA_export_metadata.md)
  : Export cell metadata to CSV
- [`TALARIA_export_ora()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TALARIA_export_ora.md)
  : Export ORA results to files
- [`TALARIA_from_seurat()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TALARIA_from_seurat.md)
  : Convert a Seurat object to a SingleCellExperiment
- [`TALARIA_generate_log_filename()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TALARIA_generate_log_filename.md)
  : Generate a timestamp-based log filename
- [`TALARIA_get_log_file()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TALARIA_get_log_file.md)
  : Get the current log file path
- [`TALARIA_is_logging()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TALARIA_is_logging.md)
  : Check whether logging is active
- [`TALARIA_load_10x()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TALARIA_load_10x.md)
  : Load a 10x Genomics MEX directory into a SingleCellExperiment
- [`TALARIA_load_h5()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TALARIA_load_h5.md)
  : Load a 10x Genomics HDF5 file into a SingleCellExperiment
- [`TALARIA_load_h5ad()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TALARIA_load_h5ad.md)
  : Load an AnnData H5AD file into a SingleCellExperiment
- [`TALARIA_load_loom()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TALARIA_load_loom.md)
  : Load a Loom file into a SingleCellExperiment
- [`TALARIA_load_sce()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TALARIA_load_sce.md)
  : Load a SingleCellExperiment object from disk
- [`TALARIA_save_sce()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TALARIA_save_sce.md)
  : Save a SingleCellExperiment object to disk
- [`TALARIA_start_log()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TALARIA_start_log.md)
  : Start console logging
- [`TALARIA_stop_log()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TALARIA_stop_log.md)
  : Stop console logging
- [`TALARIA_to_bpcells()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TALARIA_to_bpcells.md)
  : Convert a SingleCellExperiment counts assay to an on-disk BPCells
  matrix
- [`TALARIA_to_seurat()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TALARIA_to_seurat.md)
  : Convert a SingleCellExperiment to a Seurat object

## Aegis

Divine shield of Athena/Zeus, protects — QC: per-cell metrics, doublet
detection (scDblFinder), cell filtering.

- [`AEGIS_compute_qc_metrics()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/AEGIS_compute_qc_metrics.md)
  : Compute per-cell QC metrics
- [`AEGIS_detect_doublets()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/AEGIS_detect_doublets.md)
  : Detect doublets with scDblFinder
- [`AEGIS_filter_cells()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/AEGIS_filter_cells.md)
  : Filter low-quality cells

## Pyri

Sacred forge fire, purifies raw ore — preprocessing: normalization, HVG
selection and tuning.

- [`PYRI_normalize()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/PYRI_normalize.md)
  : Normalise a SingleCellExperiment
- [`PYRI_select_hvg()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/PYRI_select_hvg.md)
  : Select highly variable genes
- [`PYRI_tune_hvg()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/PYRI_tune_hvg.md)
  : Sweep HVG count to find a data-driven optimum

## Talos

Bronze automaton that traversed Crete — embedding & clustering: PCA,
Harmony, UMAP, tSNE, scVI, SNN graph construction/clustering, parameter
sweeps.

- [`TALOS_build_graph()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TALOS_build_graph.md)
  : Build a shared nearest-neighbour graph
- [`TALOS_cluster()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TALOS_cluster.md)
  : Graph-based clustering
- [`TALOS_run_harmony()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TALOS_run_harmony.md)
  : Harmony batch integration
- [`TALOS_run_pca()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TALOS_run_pca.md)
  : Principal component analysis
- [`TALOS_run_scvi()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TALOS_run_scvi.md)
  : scVI latent embedding
- [`TALOS_run_tsne()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TALOS_run_tsne.md)
  : tSNE embedding
- [`TALOS_run_umap()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TALOS_run_umap.md)
  : UMAP embedding
- [`TALOS_tune_k()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TALOS_tune_k.md)
  : Sweep k values for SNN graph construction
- [`TALOS_tune_pc()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TALOS_tune_pc.md)
  : Sweep number of PCA components
- [`TALOS_tune_resolution()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TALOS_tune_resolution.md)
  : Sweep clustering resolution
- [`TALOS_tune_tsne()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TALOS_tune_tsne.md)
  : Sweep tSNE perplexity
- [`TALOS_tune_umap()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TALOS_tune_umap.md)
  : Sweep UMAP parameters

## Pandora

Fashioned by Hephaestus, gifted by all gods — annotation: SingleR,
scType, CellTypist, manual label assignment, ortholog conversion.

- [`PANDORA_annotate_celltypist()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/PANDORA_annotate_celltypist.md)
  : Cell type annotation via CellTypist
- [`PANDORA_annotate_sctype()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/PANDORA_annotate_sctype.md)
  : Marker-based cell type annotation via scType
- [`PANDORA_annotate_singler()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/PANDORA_annotate_singler.md)
  : Reference-based cell type annotation via SingleR
- [`PANDORA_assign_labels()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/PANDORA_assign_labels.md)
  : Manually assign cell type labels to clusters
- [`PANDORA_celltypist_markers()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/PANDORA_celltypist_markers.md)
  : CellTypist's own classifier-driving genes for a cell type
- [`PANDORA_convert_orthologs()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/PANDORA_convert_orthologs.md)
  : Convert gene symbols to orthologs
- [`PANDORA_list_celltypist_models()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/PANDORA_list_celltypist_models.md)
  : List available CellTypist models
- [`PANDORA_list_references()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/PANDORA_list_references.md)
  : List available annotation references
- [`PANDORA_summarise_labels()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/PANDORA_summarise_labels.md)
  : Summarise cell type label composition

## Keraunos

Thunderbolts forged by Hephaestus, decisive — expression: cluster
markers, cell-level and pseudobulk (DESeq2) DE, GSEA/ORA (single and
pairwise), propeller proportion testing.

- [`KERAUNOS_contrast_groups()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/KERAUNOS_contrast_groups.md)
  : Cell-level contrast between groups (quick path)
- [`KERAUNOS_de_pseudobulk()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/KERAUNOS_de_pseudobulk.md)
  : Pseudobulk differential expression via DESeq2
- [`KERAUNOS_fetch_marker_genes()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/KERAUNOS_fetch_marker_genes.md)
  : Fetch known cell type marker genes from a curated database
- [`KERAUNOS_find_markers()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/KERAUNOS_find_markers.md)
  : Identify cluster marker genes
- [`KERAUNOS_gsea()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/KERAUNOS_gsea.md)
  : Gene set enrichment analysis via fgsea
- [`KERAUNOS_gsea_markers_pairwise()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/KERAUNOS_gsea_markers_pairwise.md)
  : Run GSEA between pairs of clusters using marker effect sizes
- [`KERAUNOS_gsea_pseudobulk()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/KERAUNOS_gsea_pseudobulk.md)
  : Run GSEA across all clusters from pseudobulk DE results
- [`KERAUNOS_list_gene_sets()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/KERAUNOS_list_gene_sets.md)
  : List available MSigDB gene set collections
- [`KERAUNOS_list_species()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/KERAUNOS_list_species.md)
  : List species available in MSigDB
- [`KERAUNOS_ora()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/KERAUNOS_ora.md)
  : Over-representation analysis via gprofiler2
- [`KERAUNOS_ora_pseudobulk()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/KERAUNOS_ora_pseudobulk.md)
  : Run ORA across all clusters from pseudobulk DE results
- [`KERAUNOS_propeller_proportions()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/KERAUNOS_propeller_proportions.md)
  : Test for differences in cell-type proportions across samples
  (propeller)
- [`KERAUNOS_rank_pairwise_markers()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/KERAUNOS_rank_pairwise_markers.md)
  : Pairwise ranked gene list from marker scores (for pairwise GSEA)
- [`KERAUNOS_score_markers()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/KERAUNOS_score_markers.md)
  : Annotation-focused cluster marker scoring
- [`KERAUNOS_summarise_de()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/KERAUNOS_summarise_de.md)
  : Summarise differential expression results

## Tripodes

Self-moving golden tripods, self-directed motion — trajectory: Monocle3
pseudotime, RNA velocity (scVelo), CytoTRACE stemness scoring.

- [`TRIPODES_assess_splicing()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TRIPODES_assess_splicing.md)
  : Assess spliced/unspliced read distribution across cells and genes
- [`TRIPODES_check_velocity_ready()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TRIPODES_check_velocity_ready.md)
  : Check whether an SCE object is ready for RNA velocity analysis
- [`TRIPODES_run_monocle()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TRIPODES_run_monocle.md)
  : Pseudotime and trajectory inference via Monocle3
- [`TRIPODES_run_velocity()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TRIPODES_run_velocity.md)
  : Estimate RNA velocity via scVelo
- [`TRIPODES_score_cytotrace_v1()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/TRIPODES_score_cytotrace_v1.md)
  : Score cells by differentiation state using CytoTRACE v1

## Aspis

Shield of Achilles, depicts everything — visualization: atlas plots,
embeddings, QC, markers, composition, velocity, pseudotime, trajectory.

- [`ASPIS_plot_atlas()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/ASPIS_plot_atlas.md)
  : Atlas visualisation: metadata tracks + colour ring + dendrogram +
  scatter
- [`ASPIS_plot_composition()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/ASPIS_plot_composition.md)
  : Plot a 100% stacked bar chart of per-sample cell-type composition
- [`ASPIS_plot_db_heatmap()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/ASPIS_plot_db_heatmap.md)
  : Marker database heatmap
- [`ASPIS_plot_elbow()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/ASPIS_plot_elbow.md)
  : Elbow plot of PCA variance explained
- [`ASPIS_plot_embedding_grid()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/ASPIS_plot_embedding_grid.md)
  : Grid of embedding plots from a parameter sweep
- [`ASPIS_plot_marker_dotplot()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/ASPIS_plot_marker_dotplot.md)
  : Marker gene dot plot
- [`ASPIS_plot_marker_heatmap()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/ASPIS_plot_marker_heatmap.md)
  : Marker gene heatmap
- [`ASPIS_plot_pca()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/ASPIS_plot_pca.md)
  : PCA scatter plot
- [`ASPIS_plot_propeller()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/ASPIS_plot_propeller.md)
  : Plot a propeller cell-type proportion test result
- [`ASPIS_plot_pseudotime()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/ASPIS_plot_pseudotime.md)
  : Plot pseudotime as a colour gradient on a 2D embedding
- [`ASPIS_plot_qc()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/ASPIS_plot_qc.md)
  : QC metric violin plots
- [`ASPIS_plot_trajectory()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/ASPIS_plot_trajectory.md)
  : Plot Monocle3 principal graph trajectory on a 2D embedding
- [`ASPIS_plot_tsne()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/ASPIS_plot_tsne.md)
  : tSNE embedding plot
- [`ASPIS_plot_umap()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/ASPIS_plot_umap.md)
  : UMAP embedding plot
- [`ASPIS_plot_velocity()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/ASPIS_plot_velocity.md)
  : Plot RNA velocity with coloured per-cell arrows and origin circles
- [`ASPIS_plot_velocity_stream()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/ASPIS_plot_velocity_stream.md)
  : Plot RNA velocity as smooth streamlines over a KDE cluster
  background

## Khalkos

Bronze, the foundational material of the forge — utilities: SCE
validation/merging, assay accessors, colour palettes, ggplot2 theme,
logging.

- [`KHALKOS_assign_metadata_colours()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/KHALKOS_assign_metadata_colours.md)
  : Assign a stable colour to every level of qualifying colData columns
- [`KHALKOS_check_packages()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/KHALKOS_check_packages.md)
  : Check optional package availability
- [`KHALKOS_default_palette()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/KHALKOS_default_palette.md)
  : Default CAULDRON colour palettes
- [`KHALKOS_get_assay()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/KHALKOS_get_assay.md)
  : Safe assay retrieval
- [`KHALKOS_log_message()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/KHALKOS_log_message.md)
  : Verbose-aware logging helper
- [`KHALKOS_merge_sce()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/KHALKOS_merge_sce.md)
  : Merge multiple SingleCellExperiment objects
- [`KHALKOS_set_assay()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/KHALKOS_set_assay.md)
  : Safe assay assignment
- [`KHALKOS_theme_cauldron()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/KHALKOS_theme_cauldron.md)
  : CAULDRON ggplot2 theme
- [`KHALKOS_validate_sce()`](https://ylefol.github.io/L-Atelier/CAULDRON/reference/KHALKOS_validate_sce.md)
  : Validate a SingleCellExperiment object
