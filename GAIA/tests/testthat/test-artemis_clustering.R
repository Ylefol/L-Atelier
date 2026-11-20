# ARTEMIS clustering of mixed (numeric + categorical) clinical-type data, via
# Gower distance. The `cluster` package is a recommended package that ships with
# R, but the tests still skip if it is missing.

# Three groups separated on two numeric variables and (with 10% noise) a
# categorical one.
make_mixed_clusters <- function(seed = 1) {
  set.seed(seed)
  n <- 90
  planted <- rep(1:3, each = 30)
  lab <- c("u", "v", "w")[planted]
  flip <- runif(n) < 0.1
  lab[flip] <- sample(c("u", "v", "w"), sum(flip), replace = TRUE)
  list(planted = planted,
       df = data.frame(x1 = rnorm(n, c(0, 10, 20)[planted], 1.5),
                       x2 = rnorm(n, c(0, 10, 0)[planted], 1.5),
                       grp = factor(lab)))
}

purity <- function(clusters, planted) {
  sum(apply(table(clusters, planted), 1, max)) / length(planted)
}

cluster_quiet <- function(df, ...) quiet_msg(ARTEMIS_cluster_mixed(df, ...))

test_that("mixed clustering recovers the planted groups", {
  skip_if_not_installed("cluster")
  m <- make_mixed_clusters()
  cm <- cluster_quiet(m$df, k = 3, method = "pam")

  expect_s3_class(cm, "artemis_cluster")
  expect_equal(cm$k, 3)
  expect_length(cm$clusters, 90)
  expect_equal(sum(cm$cluster_sizes), 90)
  expect_gte(purity(cm$clusters, m$planted), 0.85)
  expect_gt(cm$silhouette_avg, 0.7)
})

test_that("hierarchical clustering gives a similar partition", {
  skip_if_not_installed("cluster")
  m <- make_mixed_clusters()
  hc <- cluster_quiet(m$df, k = 3, method = "hierarchical")
  expect_gte(purity(hc$clusters, m$planted), 0.85)
})

test_that("k can be chosen automatically, and the choice is explained", {
  skip_if_not_installed("cluster")
  m <- make_mixed_clusters()
  auto <- cluster_quiet(m$df, k = "auto", k_range = 2:6)
  expect_equal(auto$k, 3)

  ev <- quiet_msg(ARTEMIS_evaluate_k_range(m$df, k_range = 2:6))
  expect_s3_class(ev, "artemis_k_evaluation")
  expect_equal(ev$optimal_k, 3)
  expect_equal(ev$metrics$k, 2:6)
  expect_equal(ev$metrics$k[which.max(ev$metrics$silhouette_avg)], 3)
})

test_that("cluster stability is high for well-separated groups", {
  skip_if_not_installed("cluster")
  m <- make_mixed_clusters()
  cm <- cluster_quiet(m$df, k = 3)
  set.seed(1)
  st <- quiet_msg(ARTEMIS_cluster_stability(cm, n_boot = 20))

  expect_s3_class(st, "artemis_stability")
  expect_gt(st$jaccard_mean, 0.8)
  expect_length(st$jaccard_per_cluster, 3)
  expect_equal(st$n_boot, 20)
})

test_that("cluster characterisation ranks the variables that define the clusters", {
  skip_if_not_installed("cluster")
  m <- make_mixed_clusters()
  cm <- cluster_quiet(m$df, k = 3)
  ch <- quiet_msg(ARTEMIS_characterize_clusters(cm))

  expect_s3_class(ch, "artemis_characterization")
  expect_equal(ch$n_clusters, 3)
  expect_setequal(ch$results$variable, c("x1", "x2", "grp"))
  expect_true(all(ch$results$p_adjusted < 0.05))
  # numeric variables are tested with a rank test, categorical ones with chi-squared
  expect_equal(ch$results$test[ch$results$variable == "grp"], "chi-sq")
  expect_equal(ch$results$test[ch$results$variable == "x1"], "kruskal")
})
