#' TBASS - Bayesian MARS with a Student's t likelihood
#'
#' Wrapper around \code{gbass()} for a Student's t error model.
#'
#' @param X An \eqn{N \times p} numeric matrix of predictor variables.
#' @param y A numeric response vector of length \eqn{N}.
#' @param df Degrees of freedom. Default is \code{5}.
#' @param ... Additional arguments passed to \code{gbass()}.
#'
#' @details
#' Uses an inverse-gamma latent-scale representation through a GIG prior on the
#' local variance components. If \code{df > 2}, the scale is set to
#' \code{(df - 2) / df} so that the marginal variance matches \code{w}.
#'
#' @export
tbass <- function(X, y, df = 5, ...) {
  v_prior <- list(type = "GIG", p = -df / 2, a = 0, b = df)
  scale <- ifelse(df > 2, (df - 2) / df, 1)
  md <- gbass(X, y, v_prior = v_prior, scale = scale, ...)
  md$df <- df
  class(md) <- c("tbass", "gbass")
  md
}


#' HBASS - Bayesian MARS with Horseshoe or Strawderman-Berger likelihood
#'
#' Wrapper around \code{gbass()} using generalized beta prime latent-scale priors.
#'
#' @param X An \eqn{N \times p} numeric matrix of predictor variables.
#' @param y A numeric response vector of length \eqn{N}.
#' @param likelihood Character string specifying the likelihood family.
#'   Use \code{"h"} for Horseshoe or \code{"sb"} for Strawderman-Berger.
#' @param ... Additional arguments passed to \code{gbass()}.
#'
#' @details
#' The Horseshoe and Strawderman-Berger likelihoods can be expressed using
#' generalized beta prime latent-scale distributions. This function provides a
#' convenient wrapper for these likelihoods. For additional flexibility, use
#' \code{gbass()} directly.
#'
#' @export
hbass <- function(X, y, likelihood = "h", ...) {
  w_prior <- list(type = "GBP", p = 1, a = 1 / 2, b = 1 / 2, prop_sigma = 0.5)

  if (likelihood == "h") {
    v_prior <- list(type = "GBP", p = 1, a = 1 / 2, b = 1 / 2, prop_sigma = 0.5)
  } else if (likelihood == "sb") {
    v_prior <- list(type = "GBP", p = 1 / 2, a = 1, b = 1 / 2, prop_sigma = 0.5)
  } else {
    stop("likelihood not recognized. See documentation.")
  }

  scale <- var(y)
  md <- gbass(X, y, w_prior = w_prior, v_prior = v_prior, scale = scale, ...)
  class(md) <- c("hbass", "gbass")
  md
}

bass_shoe <- hbass


#' QBASS - Bayesian MARS with an asymmetric Laplace likelihood
#'
#' Wrapper around \code{gbass()} for quantile regression.
#'
#' @param X An \eqn{N \times p} numeric matrix of predictor variables.
#' @param y A numeric response vector of length \eqn{N}.
#' @param q Quantile of interest. Default is \code{0.5} for median regression.
#' @param prop_sigma_v Proposal standard deviation for local variance updates.
#'   This is retained for backward compatibility.
#' @param w_prior Prior for the global variance factor.
#' @param ... Additional arguments passed to \code{gbass()}.
#'
#' @details
#' Performs quantile regression for quantile \code{q} using the asymmetric
#' Laplace representation. For many quantiles, fitting separate models in
#' parallel may be convenient.
#'
#' @export
qbass <- function(X, y, q = 0.5, prop_sigma_v = 0.25,
                  w_prior = list(type = "GIG", p = 0, a = 0, b = 0, prop_sigma = 0.2), ...) {
  v_prior <- list(type = "GIG", p = 1, a = 2, b = 0, prop_sigma = prop_sigma_v)
  scale <- 2 / (q * (1 - q))
  m_beta <- 1 / q - 1 / (1 - q)
  s_beta <- 0

  md <- gbass(
    X, y,
    v_prior = v_prior,
    w_prior = w_prior,
    scale = scale,
    m_beta = m_beta,
    s_beta = s_beta,
    ...
  )
  md$q <- q
  class(md) <- c("qbass", "gbass")
  md
}







