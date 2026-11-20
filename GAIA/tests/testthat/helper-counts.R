# Shared fixture: simulated bulk count data with a known answer.
#
# 200 features, 8 samples (4 'control', 4 'treated', alternating across two
# batches). Features 1-10 are up in 'treated' (x4), features 11-20 are down
# (/4); the other 180 are unchanged. Each sample has its own sequencing depth,
# which the normalisation should recover.
make_test_counts <- function(n_features = 200, n_per_group = 4, seed = 123) {
  set.seed(seed)
  n_samples <- 2 * n_per_group
  group <- rep(c("control", "treated"), each = n_per_group)
  depth <- rep(c(1, 1.5, 0.7, 2), length.out = n_samples)   # true size factors

  base <- rgamma(n_features, shape = 2, rate = 0.02) + 20     # mean count per feature
  fc <- rep(1, n_features)
  fc[1:10]  <- 4
  fc[11:20] <- 1 / 4

  mu <- outer(base, depth)
  mu[, group == "treated"] <- mu[, group == "treated"] * fc
  counts <- matrix(rnbinom(length(mu), mu = mu, size = 20), nrow = n_features)
  rownames(counts) <- paste0("feat", seq_len(n_features))
  colnames(counts) <- paste0("S", seq_len(n_samples))

  targets <- data.frame(
    group = group,
    batch = rep(c("b1", "b2"), length.out = n_samples),
    row.names = colnames(counts),
    stringsAsFactors = FALSE
  )
  list(counts = counts, targets = targets, depth = depth,
       up = rownames(counts)[1:10], down = rownames(counts)[11:20],
       null = rownames(counts)[21:n_features])
}

quiet <- function(expr) {
  # these functions print progress with cat(); keep test output readable
  utils::capture.output(res <- force(expr))
  res
}
