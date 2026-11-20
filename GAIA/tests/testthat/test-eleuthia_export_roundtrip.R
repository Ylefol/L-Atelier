# Export -> load round trips. What ELEUTHIA writes to disk must be what DEMETER
# reads back, because analyses are routinely resumed from exported results.

export_experiments <- function(experiments = c("exp1", "exp2")) {
  d <- make_test_counts()
  norm <- quiet(ARTEMIS_normalize_counts(d$counts, d$targets, verbose = FALSE))
  dea  <- quiet(ARTEMIS_differential_counts(norm, "control", "treated", verbose = FALSE))
  sel  <- ARTEMIS_select_de_genes(dea, l2fc_thresh = 1, verbose = FALSE)

  parent <- tempfile("exports_"); dir.create(parent)
  for (e in experiments) {
    quiet(ELEUTHIA_export_dea_results(file.path(parent, e), dea_result = dea,
                                      gene_selection = sel, sample_info = d$targets,
                                      save_plots = FALSE, prefix = "dea"))
  }
  list(parent = parent, dea = dea, sel = sel, sim = d)
}

test_that("an exported DEA result is written as csv, rds and metadata", {
  ex <- export_experiments("exp1"); on.exit(unlink(ex$parent, recursive = TRUE), add = TRUE)
  files <- list.files(file.path(ex$parent, "exp1"))

  expect_true(all(c("dea_all_results.csv", "dea_significant.csv", "dea_selected_genes.csv",
                    "dea_metadata.txt", "dea_dea_result.rds") %in% files))
})

test_that("load_dea returns one DEA result per experiment, identical to the original", {
  ex <- export_experiments(); on.exit(unlink(ex$parent, recursive = TRUE), add = TRUE)
  loaded <- quiet(DEMETER_load_dea(ex$parent))

  expect_named(loaded, c("exp1", "exp2"))
  expect_s3_class(loaded$exp1, "artemis_dea")
  expect_equal(loaded$exp1$results, ex$dea$results)
  expect_equal(loaded$exp1$comparison, "treated vs control")
  expect_equal(loaded$exp1$norm_counts, ex$dea$norm_counts)
})

test_that("load_dea can load a subset of experiments", {
  ex <- export_experiments(); on.exit(unlink(ex$parent, recursive = TRUE), add = TRUE)
  loaded <- quiet(DEMETER_load_dea(ex$parent, experiments = "exp2"))
  expect_named(loaded, "exp2")
})

test_that("load_dea falls back to the csv when the rds is missing", {
  ex <- export_experiments("exp1"); on.exit(unlink(ex$parent, recursive = TRUE), add = TRUE)
  unlink(file.path(ex$parent, "exp1", "dea_dea_result.rds"))

  loaded <- quiet(DEMETER_load_dea(ex$parent))
  expect_true(isTRUE(loaded$exp1$.reconstructed))
  expect_equal(loaded$exp1$results$log2FoldChange, ex$dea$results$log2FoldChange)
  expect_equal(loaded$exp1$results$feature_id, ex$dea$results$feature_id)
  expect_null(loaded$exp1$dds)          # the full DESeq object cannot be rebuilt from a csv
})

test_that("load_dea explains what is wrong with a bad directory", {
  expect_error(quiet(DEMETER_load_dea(file.path(tempdir(), "no_such_dir"))),
               "Directory not found")
  empty <- tempfile("empty_"); dir.create(empty); on.exit(unlink(empty, recursive = TRUE), add = TRUE)
  expect_error(quiet(DEMETER_load_dea(empty)), "No experiment subdirectories")
})

test_that("load_dea_genes gives the union of significant genes across experiments", {
  ex <- export_experiments(); on.exit(unlink(ex$parent, recursive = TRUE), add = TRUE)
  genes <- quiet(DEMETER_load_dea_genes(ex$parent))

  expect_setequal(genes, as.character(ex$sel))          # same genes in both experiments
  expect_equal(attr(genes, "n_per_experiment"), c(exp1 = 20L, exp2 = 20L))
  expect_named(attr(genes, "per_experiment"), c("exp1", "exp2"))
})

test_that("the three gene sources agree when thresholds match", {
  ex <- export_experiments(); on.exit(unlink(ex$parent, recursive = TRUE), add = TRUE)
  sig <- quiet(DEMETER_load_dea_genes(ex$parent, source = "significant"))
  sel <- quiet(DEMETER_load_dea_genes(ex$parent, source = "selected"))
  all <- quiet(DEMETER_load_dea_genes(ex$parent, source = "all", l2fc_thresh = 1, p_thresh = 0.05))

  expect_setequal(sig, sel)
  expect_setequal(sig, all)
})

# ---- enrichment --------------------------------------------------------------

make_enrichment_object <- function() {
  combined <- data.frame(
    module = rep(c("m1", "m2"), each = 4), source = rep(c("GO:BP", "KEGG"), 4),
    term_id = sprintf("T:%03d", 1:8), term_name = paste("term", 1:8),
    p_value = c(0.001, 0.01, 0.02, 0.03, 0.002, 0.004, 0.2, 0.5),
    term_size = c(50, 40, 30, 20, 60, 70, 80, 90), query_size = 10,
    intersection_size = c(5, 4, 3, 2, 6, 5, 4, 3), precision = 0.5, recall = 0.1,
    stringsAsFactors = FALSE)
  structure(list(combined = combined,
                 results = split(combined, combined$module),
                 parameters = list(organism = "hsapiens")),
            class = c("gost_enrichment", "list"))
}

test_that("an exported enrichment loads back from rds with full fidelity", {
  enr <- make_enrichment_object()
  out <- tempfile("enr_"); on.exit(unlink(out, recursive = TRUE), add = TRUE)
  quiet(ELEUTHIA_export_enrichment(enr, out, save_plots = FALSE, save_full_gostplot = FALSE))

  back <- quiet(DEMETER_load_enrichment(file.path(out, "enrichment.rds")))
  expect_s3_class(back, "gost_enrichment")
  expect_equal(back$combined, enr$combined)
})

test_that("an exported enrichment can be reconstructed from the csv alone", {
  enr <- make_enrichment_object()
  out <- tempfile("enr_"); on.exit(unlink(out, recursive = TRUE), add = TRUE)
  quiet(ELEUTHIA_export_enrichment(enr, out, save_plots = FALSE, save_full_gostplot = FALSE))

  back <- quiet(DEMETER_load_enrichment(file.path(out, "enrichment_combined.csv")))
  expect_s3_class(back, "gost_enrichment")
  expect_equal(nrow(back$combined), 8)
  expect_equal(back$combined$term_name, enr$combined$term_name)
  expect_equal(back$combined$p_value, enr$combined$p_value)

  # and it works with the same filter used on live results
  top <- APOLLO_filter_enrichment(back, top_n = 1)
  expect_equal(nrow(top), 4)
})

test_that("load_enrichment reports a missing file", {
  expect_error(quiet(DEMETER_load_enrichment(file.path(tempdir(), "none.rds"))),
               "File not found")
})
