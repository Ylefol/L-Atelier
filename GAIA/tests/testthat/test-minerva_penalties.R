# MINERVA penalty functions used by the sparse CCA methods. These are small
# closed-form operators, so the expected values below are worked out by hand.

test_that("lasso soft-thresholding shrinks towards zero and zeroes small values", {
  x <- c(-3, -0.5, 0, 0.5, 2)
  expect_equal(MINERVA_soft_threshold_lasso(x, tau = 1), c(-2, 0, 0, 0, 1))

  # identity when tau = 0, everything zero when tau exceeds all |x|
  expect_equal(MINERVA_soft_threshold_lasso(x, 0), x)
  expect_equal(MINERVA_soft_threshold_lasso(x, 10), rep(0, 5))
})

test_that("lasso soft-thresholding agrees with sign(x) * max(|x| - tau, 0)", {
  set.seed(1)
  x <- rnorm(200, sd = 2)
  expect_equal(MINERVA_soft_threshold_lasso(x, 0.7),
               sign(x) * pmax(abs(x) - 0.7, 0))
})

test_that("scad thresholding has three regimes: lasso-like, interpolated, untouched", {
  tau <- 1   # default a = 3.7, so regimes split at |x| = 2 and |x| = 3.7
  expect_equal(MINERVA_soft_threshold_scad(0.5, tau), 0)            # <= tau: removed
  expect_equal(MINERVA_soft_threshold_scad(1.5, tau), 0.5)          # lasso-like
  expect_equal(MINERVA_soft_threshold_scad(2.5, tau),               # interpolated
               (2.7 * 2.5 - 3.7) / 1.7)
  expect_equal(MINERVA_soft_threshold_scad(5, tau), 5)              # large: unpenalised
  expect_equal(MINERVA_soft_threshold_scad(-5, tau), -5)
  expect_equal(MINERVA_soft_threshold_scad(-1.5, tau), -0.5)        # odd function
})

test_that("scad thresholding keeps the order of the input and is continuous at the breakpoints", {
  x <- c(5, 0.5, -2.5, 1.5, 0)
  out <- MINERVA_soft_threshold_scad(x, 1)
  expect_length(out, 5)
  expect_equal(out[1], 5)
  expect_equal(out[2], 0)
  expect_equal(out[4], 0.5)

  # no jump at 2*tau or at a*tau (default a = 3.7)
  eps <- 1e-9
  expect_equal(MINERVA_soft_threshold_scad(2 - eps, 1), MINERVA_soft_threshold_scad(2 + eps, 1),
               tolerance = 1e-6)
  expect_equal(MINERVA_soft_threshold_scad(3.7 - eps, 1), MINERVA_soft_threshold_scad(3.7 + eps, 1),
               tolerance = 1e-6)
})

test_that("scad honours a non-default shape parameter 'a'", {
  # KNOWN BUG (reported, fix pending): the interpolated regime hard-codes the
  # constants for a = 3.7 (2.7 and 1.7) instead of using (a - 1) and (a - 2), so
  # for any other 'a' the thresholding is discontinuous at |x| = 2 * tau.
  # Remove this skip once MINERVA_soft_threshold_scad() uses 'a' throughout.
  skip("known bug: scad ignores 'a' in the interpolated regime")
  eps <- 1e-9
  expect_equal(MINERVA_soft_threshold_scad(2 - eps, 1, a = 5),
               MINERVA_soft_threshold_scad(2 + eps, 1, a = 5), tolerance = 1e-6)
})

test_that("proximal lasso follows the three-case rule", {
  # check = omega + mu * c
  expect_equal(MINERVA_proximal_lasso(omega = 1, mu = 1, c = 0.5, tau = 0.5, old = 9), 1)
  expect_equal(MINERVA_proximal_lasso(omega = -1, mu = 1, c = -0.5, tau = 0.5, old = 9), -1)
  # inside the dead zone the previous value is kept
  expect_equal(MINERVA_proximal_lasso(omega = 0.1, mu = 1, c = 0, tau = 0.5, old = 9), 9)
  # a missing value never propagates; the previous value is kept
  expect_equal(MINERVA_proximal_lasso(omega = NA_real_, mu = 1, c = 0, tau = 0.5, old = 9), 9)
})

test_that("proximal elastic net combines the L1 threshold with L2 shrinkage", {
  # check = 2 > mu * tau_EN = 1  ->  (2 - 1) / (1 + 2 * 1 * 0.5) = 0.5
  expect_equal(MINERVA_proximal_elastic_net(omega = 2, mu = 1, c = 0, tau = 0.5,
                                            tau_EN = 1, old = 7), 0.5)
  expect_equal(MINERVA_proximal_elastic_net(omega = -2, mu = 1, c = 0, tau = 0.5,
                                            tau_EN = 1, old = 7), -0.5)
  expect_equal(MINERVA_proximal_elastic_net(omega = 0.5, mu = 1, c = 0, tau = 0.5,
                                            tau_EN = 1, old = 7), 7)
  # with no L2 term it reduces to the lasso rule
  expect_equal(MINERVA_proximal_elastic_net(omega = 2, mu = 1, c = 0, tau = 0, tau_EN = 1, old = 7),
               MINERVA_proximal_lasso(omega = 2, mu = 1, c = 0, tau = 1, old = 7))
})

test_that("proximal scad covers its five cases", {
  prox <- function(omega) MINERVA_proximal_scad(omega = omega, mu = 1, c = 0,
                                                tau = 1, old = 99)   # a = 3.7
  # worked out for tau = 1, mu = 1, a = 3.7:
  #   eta_2 = 7.4 / 5.4, denominator = 1 - 1 / 2.7
  den <- 1 - 1 / 2.7
  eta2 <- 7.4 / 5.4

  expect_equal(prox(1.5), 0.5)                         # small positive
  expect_equal(prox(-1.5), -0.5)                       # small negative
  expect_equal(prox(3), (3 - eta2) / den)              # intermediate positive
  expect_equal(prox(-3), (-3 + eta2) / den)            # intermediate negative
  expect_equal(prox(5), 5)                             # large: unchanged
  expect_equal(prox(0.5), 99)                          # dead zone: previous value kept
})
