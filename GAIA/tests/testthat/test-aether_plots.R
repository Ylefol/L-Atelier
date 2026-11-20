# AETHER: visualisation. The plots are built from the real output of the
# analysis functions (not hand-made tables), so a change in an analysis result
# that breaks plotting is caught here. Each plot must build without error and
# carry the data it should; how it LOOKS is not tested.

expect_builds <- function(p) {
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
}

# ---- colours -------------------------------------------------------------------

test_that("default_palette returns named hex colours for omics types and groups", {
  om <- AETHER_default_palette("omics")
  expect_true(all(c("ATACseq", "RNAseq", "CHIPseq") %in% names(om)))
  expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", om)))
  expect_false(anyDuplicated(om) > 0)

  expect_gt(length(AETHER_default_palette("groups")), 0)
  expect_error(AETHER_default_palette("zzz"), "type must be 'omics' or 'groups'")
})

test_that("generate_gradient steps from the base colour towards light or dark", {
  light <- AETHER_generate_gradient("#336699", n = 5)
  expect_length(light, 5)
  expect_equal(light[1], "#336699")                      # starts at the base colour
  expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", light)))
  # lighter each step: mean RGB increases
  expect_false(is.unsorted(colMeans(grDevices::col2rgb(light))))

  dark <- AETHER_generate_gradient("#336699", n = 3, direction = "to_dark")
  expect_equal(dark[3], "#336699")
  expect_error(AETHER_generate_gradient("#336699", n = 0), "n must be at least 1")
})

# ---- PC / metadata association -----------------------------------------------

test_that("the PC association heatmap has one cell per PC and variable", {
  x <- make_pc_data()
  pa <- quiet_msg(ARTEMIS_pc_metadata_association(x$m, x$meta, n_pcs = 5,
                                                  sample_col = "sample", verbose = FALSE))
  p <- AETHER_plot_pc_association(pa)

  expect_builds(p)
  expect_equal(nrow(p$data), 5 * 3)
  expect_named(p$data, c("PC", "Variable", "log10p"))
  # the batch effect on PC1 is the strongest cell
  top <- p$data[which.max(p$data$log10p), ]
  expect_equal(as.character(top$PC), "PC1")
  expect_equal(as.character(top$Variable), "batch")

  expect_error(AETHER_plot_pc_association(list()), "must be an artemis_pc_assoc object")
})

# ---- PART clustering plots ---------------------------------------------------

part_fixture <- function() {
  x <- make_divergence()
  pr <- part_quiet(x$mat, q = 0.25, min_size = 8, B = 50, seed = 1)
  dv <- quiet_msg(ARTEMIS_part_trajectory_divergence(pr, x$info, group_col = "group",
                                                     time_col = "timepoint",
                                                     perm = 200, seed = 1))
  list(x = x, pr = pr, dv = dv)
}

test_that("the divergence plot shows every cluster and marks the significant one", {
  f <- part_fixture()
  p <- AETHER_plot_cluster_divergence(f$dv)

  expect_builds(p)
  expect_equal(nrow(p$data), 2)
  expect_equal(sum(grepl("<", p$data$significant)), 1)   # one cluster below the p_adj cutoff
  expect_true(all(c("cluster", "divergence_stat", "p_adj", "neg_log10_p") %in% names(p$data)))
})

test_that("cluster mean, group-mean and trajectory plots build from a PART result", {
  f <- part_fixture()

  means <- AETHER_plot_cluster_means(f$pr, f$x$info)
  expect_builds(means)
  expect_equal(nrow(means$data), 2 * 4 * 2)           # clusters x timepoints x groups

  expect_builds(AETHER_plot_cluster_group_means(f$pr, f$x$info))
  expect_builds(AETHER_plot_cluster_trajectories(f$pr, f$x$info))
})

test_that("the PART heatmap draws", {
  f <- part_fixture()
  h <- AETHER_plot_part_heatmap(f$pr, f$x$info, sample_col = NULL,
                                group_col = "group", time_col = "timepoint")
  expect_s4_class(h, "HeatmapList")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(ComplexHeatmap::draw(h))
})

# ---- time-series summary -----------------------------------------------------

test_that("the time-series summary counts significant genes per comparison and direction", {
  d <- make_timeseries()
  nt <- quiet_msg(ARTEMIS_normalize_timeseries(d$counts, d$targets, verbose = FALSE))
  te <- quiet_msg(ARTEMIS_timeseries_temporal(nt, verbose = FALSE))
  p <- AETHER_plot_timeseries_summary(te)

  expect_builds(p)
  expect_equal(nrow(p$data), 4)                         # 2 comparisons x up / down
  expect_named(p$data, c("label", "direction", "count", "abs_count"))
  # each comparison finds ~10 genes, all of them increases
  up <- p$data[p$data$direction == "Up", ]
  expect_true(all(up$count >= 8))
  expect_true(all(p$data$count[p$data$direction == "Down"] == 0))
})

# ---- CCA weight plots --------------------------------------------------------

test_that("feature-weight plots mark which features carry weight", {
  s <- simulate_cca()
  res <- quiet(HEPHAESTUS_convCCA(s$X[[1]], s$X[[2]], 0.15, 0.15, nIter = 100))

  p <- quiet(AETHER_plot_feature_weights(res$W_1, "sCCA"))
  expect_builds(p)
  expect_equal(nrow(p$data), 40)
  expect_named(p$data, c("feature_index", "weight", "is_nonzero"))
  expect_equal(sum(p$data$is_nonzero), sum(abs(res$W_1) > 1e-6))

  multi <- quiet(AETHER_plot_multi_feature_weights(list(res$W_1, res$W_2), "sCCA",
                                                   dataset_names = c("A", "B")))
  expect_builds(multi)
})
