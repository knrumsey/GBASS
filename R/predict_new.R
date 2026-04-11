#' Predict method for GBASS objects
#'
#' Returns posterior draws of either:
#'   (i) the linear predictor / mean surface, or
#'   (ii) the full posterior predictive distribution.
#'
#' @param object an object of class "gbass" (including subclasses like
#'   "tbass", "qbass", "nwbass")
#' @param newdata a matrix of predictor variables. Defaults to training inputs.
#' @param mcmc.use optional vector indicating which posterior draws to use.
#' @param predictive logical. If TRUE, return posterior predictive draws.
#'   If FALSE, return draws of the linear predictor.
#' @param bias_correct logical. Ignored unless predictive = FALSE.
#'   If TRUE, return the posterior mean response rather than just the linear predictor.
#'
#' @details
#' If predictive = FALSE and bias_correct = FALSE, this returns draws of B(x)a.
#'
#' If predictive = FALSE and bias_correct = TRUE, this returns draws of
#'   B(x)a + E(error | posterior draw),
#' i.e. the mean response under the fitted GBASS error model.
#'
#' If predictive = TRUE, this returns posterior predictive draws by simulating
#' a fresh latent local variance v_new and Gaussian draw for each posterior sample.
#'
#' For qbass objects, bias_correct = TRUE is usually not what you want, because
#' qbass is typically being used for quantile regression rather than mean regression.
#'
#' Currently, posterior predictive draws are implemented for GIG-based models.
#' If the fitted object uses a GBP prior for v, predictive = TRUE will stop.
#'
#' @return a matrix with rows corresponding to posterior draws and columns
#'   corresponding to rows of newdata.
#'
#' @export
predict.gbass <- function(object, newdata = NULL, mcmc.use = NULL,
                          predictive = TRUE, bias_correct = FALSE) {

  if (is.null(newdata)) {
    newdata <- object$X
  }
  if (is.null(dim(newdata))) {
    newdata <- matrix(newdata, ncol = ncol(object$X))
  }

  N <- nrow(newdata)

  if (is.null(mcmc.use)) {
    if (!is.null(object$M)) {
      mcmc.use <- seq_along(object$M)
    } else if (!is.null(object$nbasis)) {
      mcmc.use <- seq_along(object$nbasis)
    } else {
      stop("Could not determine available MCMC draws from object.")
    }
  }

  if (predictive && isTRUE(bias_correct)) {
    warning("bias_correct is ignored when predictive = TRUE.")
  }

  if (inherits(object, "qbass") && !predictive && isTRUE(bias_correct)) {
    warning("bias_correct=TRUE is usually not what you want for qbass/quantile regression.")
  }

  # small helper to build basis matrix at newdata for one posterior draw
  build_basis_gbass <- function(object, draw_index, newdata) {
    tX <- t(newdata)
    basis_curr <- object$lookup[unlist(object$basis[[draw_index]])]
    B_curr <- matrix(1, nrow = nrow(newdata), ncol = length(basis_curr) + 1)

    if (length(basis_curr) > 0) {
      for (m in seq_along(basis_curr)) {
        thetaB <- basis_curr[[m]]
        B_curr[, m + 1] <- makeBasis(thetaB$s, thetaB$u, thetaB$t, tX, 1)
      }
    }
    B_curr
  }

  # helper: draw from the prior for a new local variance factor v_new
  draw_v_new <- function(object, draw_index, n) {

    # nwbass: assume intended NW model from paper:
    # v ~ GIG(-1/2, gamma^2, 1)
    if (inherits(object, "nwbass")) {
      if (is.null(object$gamma)) {
        stop("nwbass object does not contain gamma draws.")
      }
      gam <- object$gamma[draw_index]
      return(replicate(n, rgig2(-1/2, gam^2, 1)))
    }

    # generic gbass / tbass / qbass path
    vp <- object$v_prior

    if (is.null(vp$type)) {
      stop("object$v_prior$type is missing.")
    }

    if (vp$type == "GIG") {
      return(replicate(n, rgig2(vp$p, vp$a, vp$b)))
    }

    stop("predictive=TRUE is currently implemented only for GIG-based v_prior objects.")
  }

  # helper: mean of v under the prior, for bias correction
  mean_v_prior <- function(object, draw_index) {

    if (inherits(object, "nwbass")) {
      if (is.null(object$gamma)) {
        stop("nwbass object does not contain gamma draws.")
      }
      gam <- object$gamma[draw_index]
      # For GIG(-1/2, gamma^2, 1), E(v) = 1/gamma
      return(1 / gam)
    }

    vp <- object$v_prior

    if (vp$type == "GIG") {
      return(mu_gig(vp$p, vp$a, vp$b))
    }

    stop("bias_correct=TRUE currently requires a GIG-based v_prior (or nwbass).")
  }

  ndraw <- length(mcmc.use)
  res <- matrix(NA_real_, nrow = ndraw, ncol = N)

  jj <- 1L
  for (i in mcmc.use) {

    B_curr <- build_basis_gbass(object, i, newdata)
    eta <- as.numeric(B_curr %*% object$a[[i]])

    # draw-specific parameters
    w_i <- object$w[i]
    beta_i <- if (!is.null(object$beta)) object$beta[i] else 0
    scale_i <- if (!is.null(object$scale)) object$scale else 1

    if (!predictive) {
      if (isTRUE(bias_correct)) {
        mu_v_i <- mean_v_prior(object, i)
        res[jj, ] <- eta + sqrt(w_i) * beta_i * mu_v_i
      } else {
        res[jj, ] <- eta
      }
    } else {
      v_new <- draw_v_new(object, i, N)
      z_new <- rnorm(N)
      eps_new <- sqrt(w_i) * (beta_i * v_new + sqrt(scale_i * v_new) * z_new)
      res[jj, ] <- eta + eps_new
    }

    jj <- jj + 1L
  }

  return(res)
}
