# Shared test fixture: a small simulated SingleCellExperiment.
#
# Three cell types (A/B/C) in two samples. Each type has its own set of
# marker genes expressed more highly in that type; the remaining genes are
# background. 10 mitochondrial ("mt-") and 10 ribosomal ("Rps") genes are
# included so the QC helpers have something to find.
make_test_sce <- function(n_per_type = 40, n_marker = 15, n_background = 150,
                          seed = 1) {
  set.seed(seed)
  types <- rep(c("A", "B", "C"), each = n_per_type)
  n_cells <- length(types)

  marker_names <- unlist(lapply(c("A", "B", "C"),
                                function(t) paste0("marker", t, "_", seq_len(n_marker))))
  gene_names <- c(marker_names,
                  paste0("gene", seq_len(n_background)),
                  paste0("mt-g", 1:10),
                  paste0("Rps", 1:10))
  n_genes <- length(gene_names)

  # baseline mean per gene, then a boost for each type's own markers
  mu <- matrix(rgamma(n_genes, shape = 2, rate = 1), n_genes, n_cells)
  for (t in c("A", "B", "C")) {
    idx <- grep(paste0("^marker", t, "_"), gene_names)
    mu[idx, types == t] <- mu[idx, types == t] * 8
  }
  # cell-to-cell library size variation
  mu <- sweep(mu, 2, rlnorm(n_cells, sdlog = 0.2), `*`)

  counts <- matrix(rpois(length(mu), mu), nrow = n_genes,
                   dimnames = list(gene_names, paste0("cell", seq_len(n_cells))))

  sce <- SingleCellExperiment::SingleCellExperiment(
    assays  = list(counts = counts),
    colData = S4Vectors::DataFrame(
      cell_type = types,
      sample_id = rep(c("s1", "s2"), length.out = n_cells),
      row.names = colnames(counts)
    )
  )
  sce
}
