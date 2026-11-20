# Shared fixtures for the PART / PC-association tests and the AETHER plot tests.

# Features with three distinct trajectories over six conditions (25 / 20 / 15).
make_trajectories <- function(seed = 3) {
  set.seed(seed)
  planted <- rep(1:3, c(25, 20, 15))
  pat <- rbind(c(0, 1, 2, 3, 4, 5), c(5, 4, 3, 2, 1, 0), c(0, 3, 0, 3, 0, 3))
  mat <- t(sapply(planted, function(k) pat[k, ] + rnorm(6, sd = 0.3)))
  dimnames(mat) <- list(sprintf("f%02d", seq_along(planted)), paste0("t", 1:6))
  list(mat = mat, planted = planted)
}

part_quiet <- function(...) quiet_msg(ARTEMIS_part(..., verbose = FALSE))

make_divergence <- function(seed = 7) {
  set.seed(seed)
  groups <- rep(c("ctrl", "trt"), each = 12)
  tps    <- rep(rep(1:4, each = 3), 2)
  n <- 24
  ids <- sprintf("s%02d", seq_len(n))
  info <- data.frame(sample_id = ids, group = groups, timepoint = tps,
                     row.names = ids, stringsAsFactors = FALSE)
  noisy <- function(ctrl, trt, nf) {
    t(sapply(seq_len(nf), function(i) ifelse(groups == "ctrl", ctrl(tps), trt(tps)) + rnorm(n, sd = 0.3)))
  }
  same     <- noisy(function(t) t, function(t) t,     20)   # identical in both groups
  diverges <- noisy(function(t) t, function(t) 5 - t, 20)   # trt falls while ctrl rises
  mat <- rbind(same, diverges)
  dimnames(mat) <- list(c(sprintf("same%02d", 1:20), sprintf("div%02d", 1:20)), ids)
  list(mat = mat, info = info)
}

make_pc_data <- function(seed = 2) {
  set.seed(seed)
  m <- matrix(rnorm(60 * 30), 60, 30)
  batch <- rep(c("b1", "b2"), each = 15)
  m[1:20, batch == "b2"] <- m[1:20, batch == "b2"] + 3      # a strong batch effect
  meta <- data.frame(sample = paste0("x", 1:30), batch = batch, age = rnorm(30),
                     noise = sample(c("u", "v"), 30, replace = TRUE),
                     stringsAsFactors = FALSE)
  dimnames(m) <- list(paste0("g", 1:60), meta$sample)
  list(m = m, meta = meta)
}
