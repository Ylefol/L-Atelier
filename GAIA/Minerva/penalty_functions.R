###############################################################################
########### Penalty Functions for Sparse CCA ###########
###############################################################################

#' Apply Soft-Thresholding with LASSO Penalty
#'
#' @description Soft-thresholding operator for LASSO (L1) penalty.
#' Shrinks values toward zero by threshold tau.
#'
#' @param x Input value or vector
#' @param tau Threshold parameter (sparsity level)
#'
#' @return Soft-thresholded value(s)
#'
#' @export
MINERVA_soft_threshold_lasso <- function(x, tau) {
    u.new <- abs(x) - tau
    u.new <- (u.new + abs(u.new)) / 2
    u.new <- u.new * sign(x)
    return(as.numeric(u.new))
}

#' Apply SCAD Penalty Thresholding
#'
#' @description Smoothly Clipped Absolute Deviation (SCAD) penalty.
#' Treats small, medium, and large coefficients differently:
#' - Small coefficients: shrunk like LASSO
#' - Medium coefficients: partially shrunk
#' - Large coefficients: no penalty
#'
#' @param x Input value or vector
#' @param tau Threshold parameter
#' @param a SCAD shape parameter (default = 3.7, from Fan & Li 2001)
#'
#' @return SCAD-thresholded value(s)
#'
#' @export
MINERVA_soft_threshold_scad <- function(x, tau, a = 3.7) {
    names(x) <- 1:length(x)

    # First part: small coefficients (LASSO-like)
    vx1 <- x[abs(x) <= 2*tau]
    u.new1 <- abs(vx1) - tau
    u.new1 <- (u.new1 + abs(u.new1)) / 2
    u.new1 <- u.new1 * sign(vx1)

    # Second part: medium coefficients (partial shrinkage)
    vx2 <- x[abs(x) > 2*tau & abs(x) <= a*tau]
    u.new2 <- (2.7*vx2 - sign(vx2)*a*tau) / 1.7

    # Third part: large coefficients (no penalty)
    vx3 <- x[abs(x) > a*tau]
    u.new3 <- vx3

    # Combine and restore original order
    u.new <- c(u.new1, u.new2, u.new3)
    u.new <- u.new[names(x)]
    u.new <- as.numeric(u.new)

    return(u.new)
}

#' Apply Proximal Operator for LASSO (ADMM Framework)
#'
#' @description Element-wise proximal operator for LASSO in ADMM.
#' Used in RelPMDCCA algorithm.
#'
#' @param omega Proximal gradient term
#' @param mu Step size parameter
#' @param c Cross-product term (gradient of correlation)
#' @param tau Sparsity tuning parameter
#' @param old Previous iteration value (fallback if update fails)
#'
#' @return Updated value
#'
#' @export
MINERVA_proximal_lasso <- function(omega, mu, c, tau, old) {
    check <- omega + mu * c
    condition_1 <- check > mu * tau
    condition_2 <- check < -mu * tau

    if (!is.na(condition_1) && condition_1) {
        updateW <- check - mu * tau
    } else if (!is.na(condition_2) && condition_2) {
        updateW <- check + mu * tau
    } else {
        updateW <- old
    }

    return(updateW)
}

#' Apply Proximal Operator for SCAD (ADMM Framework)
#'
#' @description Element-wise proximal operator for SCAD in ADMM.
#' Used in RelPMDCCA algorithm.
#'
#' @param omega Proximal gradient term
#' @param mu Step size parameter
#' @param c Cross-product term
#' @param tau Sparsity tuning parameter
#' @param a SCAD shape parameter (default = 3.7)
#' @param old Previous iteration value (fallback)
#'
#' @return Updated value
#'
#' @export
MINERVA_proximal_scad <- function(omega, mu, c, tau, a = 3.7, old) {
    eta_1 <- -1 / (2 * (a - 1))
    eta_2 <- (2 * a * tau) / (2 * (a - 1))
    check <- omega + mu * c

    # Five cases for SCAD proximal operator
    condition_1 <- (0 < check - mu*tau) && (check - mu*tau <= tau)
    condition_2 <- (-tau <= check + mu*tau) && (check + mu*tau < 0)
    condition_3 <- (tau < (check - mu*eta_2)/(1 + 2*mu*eta_1)) &&
                   ((check - mu*eta_2)/(1 + 2*mu*eta_1) <= a*tau)
    condition_4 <- (-a*tau <= (check + mu*eta_2)/(1 + 2*mu*eta_1)) &&
                   ((check + mu*eta_2)/(1 + 2*mu*eta_1) <= -tau)
    condition_5 <- abs(check) > a*tau

    if (!is.na(condition_1) && condition_1) {
        updateW <- check - mu*tau
    } else if (!is.na(condition_2) && condition_2) {
        updateW <- check + mu*tau
    } else if (!is.na(condition_3) && condition_3) {
        updateW <- (check - mu*eta_2) / (1 + 2*mu*eta_1)
    } else if (!is.na(condition_4) && condition_4) {
        updateW <- (check + mu*eta_2) / (1 + 2*mu*eta_1)
    } else if (!is.na(condition_5) && condition_5) {
        updateW <- check
    } else {
        updateW <- old
    }

    return(updateW)
}

#' Apply Proximal Operator for Elastic-Net (ADMM Framework)
#'
#' @description Element-wise proximal operator for Elastic-Net penalty.
#' Combines L1 (LASSO) and L2 (Ridge) penalties.
#'
#' @param omega Proximal gradient term
#' @param mu Step size parameter
#' @param c Cross-product term
#' @param tau L2 penalty parameter (Ridge component)
#' @param tau_EN L1 penalty parameter (LASSO component)
#' @param old Previous iteration value (fallback)
#'
#' @return Updated value
#'
#' @export
MINERVA_proximal_elastic_net <- function(omega, mu, c, tau, tau_EN, old) {
    check <- omega + mu * c
    condition_1 <- check > mu * tau_EN
    condition_2 <- check < -mu * tau_EN

    if (!is.na(condition_1) && condition_1) {
        updateW <- (check - mu * tau_EN) / (1 + 2 * mu * tau)
    } else if (!is.na(condition_2) && condition_2) {
        updateW <- (check + mu * tau_EN) / (1 + 2 * mu * tau)
    } else {
        updateW <- old
    }

    return(updateW)
}
