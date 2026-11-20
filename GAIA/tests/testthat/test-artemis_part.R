# ARTEMIS PART clustering, trajectory divergence and PC-metadata association.

test_that("PART recovers three planted trajectory clusters", {
  x <- make_trajectories()
  pr <- part_quiet(x$mat, q = 0.25, min_size = 8, B = 50, seed = 1)

  expect_s3_class(pr, "artemis_part")
  expect_equal(pr$n_clusters, 3)
  expect_equal(pr$n_outliers, 0)
  expect_setequal(names(pr$clusters), rownames(x$mat))

  # every planted cluster maps onto exactly one PART cluster, and vice versa
  tab <- table(pr$clusters[rownames(x$mat)], x$planted)
  expect_equal(sum(tab > 0), 3)
  expect_equal(sort(as.vector(tab[tab > 0])), c(15, 20, 25))
})

test_that("PART is reproducible for a given seed", {
  x <- make_trajectories()
  a <- part_quiet(x$mat, B = 30, seed = 5)
  b <- part_quiet(x$mat, B = 30, seed = 5)
  expect_identical(a$clusters, b$clusters)
})

test_that("features in clusters smaller than min_size become outliers (cluster C0)", {
  x <- make_trajectories()
  pr <- part_quiet(x$mat, min_size = 18, B = 30, seed = 1)   # the 15-feature cluster is too small

  expect_equal(pr$n_clusters, 2)
  expect_equal(pr$n_outliers, 15)
  expect_equal(sum(pr$clusters == "C0"), 15)
  expect_true(all(pr$clusters[x$planted == 3] == "C0"))
  real <- table(pr$clusters[pr$clusters != "C0"])
  expect_true(all(real >= 18))
})

test_that("part_summary prints the cluster structure", {
  x <- make_trajectories()
  pr <- part_quiet(x$mat, B = 30, seed = 1)
  out <- utils::capture.output(ARTEMIS_part_summary(pr))
  expect_gt(length(out), 0)
  expect_true(any(grepl("C1", out)))
})

# ---- trajectory divergence ---------------------------------------------------

test_that("trajectory divergence flags the cluster whose groups follow different paths", {
  x <- make_divergence()
  pr <- part_quiet(x$mat, q = 0.25, min_size = 8, B = 50, seed = 1)
  expect_equal(pr$n_clusters, 2)

  dv <- quiet_msg(ARTEMIS_part_trajectory_divergence(pr, x$info, group_col = "group",
                                                     time_col = "timepoint",
                                                     perm = 300, seed = 1))
  expect_s3_class(dv, "artemis_part_divergence")
  res <- dv$results
  expect_equal(nrow(res), 2)
  expect_named(res, c("cluster", "n_genes", "divergence_stat", "p_value", "p_adj"))

  # which PART cluster holds the diverging features?
  div_cluster <- unique(pr$clusters[grepl("^div", names(pr$clusters))])
  same_cluster <- unique(pr$clusters[grepl("^same", names(pr$clusters))])
  expect_length(div_cluster, 1)
  expect_length(same_cluster, 1)

  d <- res[res$cluster == div_cluster, ]
  s <- res[res$cluster == same_cluster, ]
  expect_lt(d$p_adj, 0.05)
  expect_gt(s$p_adj, 0.05)
  expect_gt(d$divergence_stat, 100 * s$divergence_stat)
  expect_equal(d$n_genes, 20)
})

test_that("trajectory divergence is reproducible for a given seed", {
  x <- make_divergence()
  pr <- part_quiet(x$mat, B = 30, seed = 1)
  f <- function() quiet_msg(ARTEMIS_part_trajectory_divergence(
    pr, x$info, group_col = "group", time_col = "timepoint", perm = 100, seed = 2))$results
  expect_identical(f(), f())
})

# ---- PC / metadata association -----------------------------------------------

test_that("the PC that carries the batch effect is associated with batch, not with noise variables", {
  x <- make_pc_data()
  pa <- quiet_msg(ARTEMIS_pc_metadata_association(x$m, x$meta, n_pcs = 5,
                                                  sample_col = "sample", verbose = FALSE))
  expect_s3_class(pa, "artemis_pc_assoc")
  expect_equal(dim(pa$pvalues), c(5L, 3L))
  expect_equal(colnames(pa$pvalues), c("batch", "age", "noise"))

  expect_lt(pa$pvalues["PC1", "batch"], 1e-4)
  expect_gt(pa$effect_sizes["PC1", "batch"], 0.5)
  expect_gt(min(pa$pvalues["PC1", c("age", "noise")]), 0.01)
  expect_equal(which.min(pa$pvalues[, "batch"]), c(PC1 = 1L))
})

test_that("categorical and numeric variables get appropriate tests", {
  x <- make_pc_data()
  pa <- quiet_msg(ARTEMIS_pc_metadata_association(x$m, x$meta, n_pcs = 3,
                                                  sample_col = "sample", verbose = FALSE))
  expect_equal(pa$test_used[["batch"]], "kruskal")
  expect_equal(pa$test_used[["noise"]], "kruskal")
  expect_equal(pa$test_used[["age"]], "spearman")
})

test_that("variance explained decreases across PCs and sums to at most 1", {
  x <- make_pc_data()
  pa <- quiet_msg(ARTEMIS_pc_metadata_association(x$m, x$meta, n_pcs = 5,
                                                  sample_col = "sample", verbose = FALSE))
  expect_length(pa$var_explained, 5)
  expect_false(is.unsorted(rev(pa$var_explained)))
  expect_lte(sum(pa$var_explained), 1)
  expect_gt(pa$var_explained[1], 0.2)    # the batch effect dominates PC1
})
