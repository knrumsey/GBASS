#' Generalized Bayesian MARS with a Normal-Wald likelihood
#'
#' Fits a generalized BMARS model under a Normal-Wald likelihood. This provides
#' flexible nonlinear regression with a unimodal, potentially skewed error model.
#'
#' @param X An \eqn{N \times p} numeric matrix of predictor variables. A numeric
#'   vector is treated as a single-column matrix.
#' @param y A numeric response vector of length \eqn{N}.
#' @param w_prior A named list specifying the prior for the global variance
#'   component. See Details.
#' @param maxInt Integer giving the maximum degree of interaction in spline basis
#'   functions. Defaults to \code{3}.
#' @param maxBasis Maximum number of basis functions.
#' @param npart Minimum number of nonzero points required for a proposed basis
#'   function. Defaults to \code{min(20, 0.1 * N)}.
#' @param nmcmc Total number of MCMC iterations.
#' @param nburn Number of initial MCMC iterations discarded as burn-in.
#' @param thin Thinning interval for retained draws.
#' @param moveProbs A length-3 vector giving probabilities for birth, death, and
#'   mutation moves.
#' @param a_tau Prior hyperparameter for \code{tau}.
#' @param b_tau Prior hyperparameter for \code{tau}. Defaults to \code{N / 2}.
#' @param a_lambda Prior hyperparameter for \code{lambda}.
#' @param b_lambda Prior hyperparameter for \code{lambda}.
#' @param m_beta Prior mean for \code{beta}.
#' @param s_beta Prior standard deviation for \code{beta}.
#' @param lag_beta Number of initial iterations for which \code{beta} is fixed at
#'   \code{m_beta}. This is often used to stabilize early sampling.
#' @param m_gamma Prior mean for \code{gamma}.
#' @param s_gamma Prior standard deviation for \code{gamma}.
#' @param scale Fixed variance-scale parameter. Defaults to \code{1}.
#' @param Iw0 Vector of positive nominal weights for interaction order in proposed
#'   basis functions. Must have length \code{maxInt}.
#' @param Zw0 Vector of positive nominal weights for variable selection in proposed
#'   basis functions. Must have length \code{ncol(X)}.
#' @param verbose Logical; should progress be printed?
#'
#' @details
#' The latent local-scale prior is fixed internally to the Normal-Wald form
#' \eqn{v_i \sim \mathrm{GIG}(-1/2, \gamma^2, 1)}. Unlike \code{gbass()},
#' \code{nwbass()} does not expose a user-specified \code{v_prior}.
#'
#' The \code{w_prior} list should contain:
#' \enumerate{
#'   \item \code{type}: either \code{"GIG"} or \code{"GBP"}.
#'   \item \code{p}, \code{a}, \code{b}: hyperparameters for the prior.
#'   \item \code{lower_bound}: optional lower bound for \code{w}. For backward
#'     compatibility, \code{lb} is also accepted.
#'   \item \code{prop_sigma}: proposal standard deviation on the log scale for
#'     Metropolis updates of \code{w}. This is only used when \code{w} is updated
#'     by Metropolis-Hastings, such as the \code{"GBP"} case.
#'   \item \code{adapt}, \code{adapt_delay}, \code{adapt_thin}: optional controls
#'     for adaptive Metropolis updates of \code{w} when applicable.
#' }
#'
#' Retained draws are taken at iterations
#' \code{nburn + 1, nburn + 1 + thin, ...}. Thus \code{nburn} is interpreted as
#' the number of initial iterations discarded as burn-in.
#'
#' @return
#' An object of class \code{c("nwbass", "gbass")} containing posterior draws and
#' fitted model information.
#'
#' @import Matrix
#' @export
#'
#' @examples
#' # Simple example
#' # n <- 200
#' # X <- lhs::maximinLHS(n, 2)
#' # y <- 6 * apply(X, 1, function(x) sin(pi * x[1]) + x[2]^2)
#' # mod <- nwbass(X, y, nmcmc = 2000, nburn = 1000)
nwbass <- function(
    X, y,
    w_prior = list(type = "GIG", p = -0.1, a = 0, b = 0.1),
    maxInt = 3, maxBasis = 1000, npart = NULL,
    nmcmc = 10000, nburn = 9000, thin = 1,
    moveProbs = rep(1 / 3, 3),
    a_tau = 1 / 2, b_tau = NULL,
    a_lambda = 1, b_lambda = 1,
    m_beta = 0, s_beta = 0, lag_beta = 1001,
    m_gamma = 1, s_gamma = 0,
    scale = 1,
    Iw0 = rep(1, maxInt), Zw0 = rep(1, ncol(X)),
    verbose = TRUE
) {
  # ---------------------------------------------------------------------------
  # Input checks and normalization
  # ---------------------------------------------------------------------------
  if (is.null(dim(X))) {
    X <- matrix(X, ncol = 1)
  } else {
    X <- as.matrix(X)
  }

  if (!is.numeric(X)) {
    stop("X must be numeric.")
  }
  if (!is.numeric(y)) {
    stop("y must be numeric.")
  }
  if (anyNA(X) || anyNA(y)) {
    stop("X and y must not contain missing values.")
  }
  if (nrow(X) != length(y)) {
    stop("nrow(X) and length(y) should match.")
  }
  if (!is.numeric(nmcmc) || length(nmcmc) != 1L || nmcmc < 1 || nmcmc != as.integer(nmcmc)) {
    stop("nmcmc must be a positive integer.")
  }
  if (!is.numeric(nburn) || length(nburn) != 1L || nburn < 0 || nburn != as.integer(nburn)) {
    stop("nburn must be a nonnegative integer.")
  }
  if (nburn >= nmcmc) {
    stop("nburn must be strictly less than nmcmc.")
  }
  if (!is.numeric(thin) || length(thin) != 1L || thin < 1 || thin != as.integer(thin)) {
    stop("thin must be a positive integer.")
  }
  if (!is.numeric(maxBasis) || length(maxBasis) != 1L || maxBasis < 1 || maxBasis != as.integer(maxBasis)) {
    stop("maxBasis must be a positive integer.")
  }
  if (!is.numeric(maxInt) || length(maxInt) != 1L || maxInt < 1 || maxInt != as.integer(maxInt)) {
    stop("maxInt must be a positive integer.")
  }
  if (!is.numeric(scale) || length(scale) != 1L || scale <= 0) {
    stop("scale must be a positive numeric scalar.")
  }
  if (!is.numeric(moveProbs) || length(moveProbs) != 3L || any(moveProbs < 0) || sum(moveProbs) <= 0) {
    stop("moveProbs must be a numeric vector of length 3 with nonnegative entries and positive sum.")
  }
  moveProbs <- moveProbs / sum(moveProbs)

  N <- nrow(X)
  p <- ncol(X)
  maxInt <- min(maxInt, p)

  if (length(Iw0) != maxInt || any(Iw0 <= 0)) {
    stop("Iw0 must have length maxInt and positive entries.")
  }
  if (length(Zw0) != p || any(Zw0 <= 0)) {
    stop("Zw0 must have length ncol(X) and positive entries.")
  }

  if (!is.numeric(lag_beta) || length(lag_beta) != 1L || lag_beta < 0 || lag_beta != as.integer(lag_beta)) {
    stop("lag_beta must be a nonnegative integer.")
  }
  if (lag_beta > nburn) {
    lag_beta <- nburn
    warning("lag_beta cannot exceed nburn. Setting lag_beta = nburn.")
  }

  if (is.null(w_prior$type) || !(w_prior$type %in% c("GIG", "GBP"))) {
    stop("w_prior$type must be either 'GIG' or 'GBP'.")
  }
  if (!is.null(w_prior$lb) && is.null(w_prior$lower_bound)) {
    w_prior$lower_bound <- w_prior$lb
  }

  if (is.null(npart)) {
    npart <- min(20, 0.1 * N)
  }
  if (is.null(b_tau)) {
    b_tau <- N / 2
  }

  keep_idx <- seq.int(from = nburn + 1L, to = nmcmc, by = thin)
  nkeep <- length(keep_idx)
  keep_iter <- rep(FALSE, nmcmc)
  keep_iter[keep_idx] <- TRUE

  # Fixed latent-scale prior for the Normal-Wald model
  v_prior <- list(type = "GIG", p = -1 / 2, a = 1, b = 1)

  # Initial moments under the current NW latent structure
  mu_v <- 1
  s2_v <- 1

  if (is.null(w_prior$lower_bound)) {
    w_prior$lower_bound <- var(y) / N
  }

  # ---------------------------------------------------------------------------
  # Allocate storage
  # ---------------------------------------------------------------------------
  a_mc <- list()
  M_mc <- lam_mc <- tau_mc <- w_mc <- bet_mc <- gam_mc <- ss <- bias_mc <- s2_mc <- rep(NA_real_, nkeep)
  v_mc <- matrix(NA_real_, nrow = nkeep, ncol = N)
  lookup <- basis_mc <- list()
  basis_index <- integer(0)
  nlook <- 1L
  kk <- 1L

  # ---------------------------------------------------------------------------
  # Initialize parameters
  # ---------------------------------------------------------------------------
  M <- 0L
  tau <- (b_tau + 0.1) / (a_tau + 0.1)
  gam <- 1
  v <- rep(1, N)
  w <- var(y) / scale
  lam <- rgamma(1, a_lambda, b_lambda)
  bet <- rnorm(1, m_beta, s_beta)
  a <- mean(y)

  s_beta_hold <- s_beta
  if (lag_beta > 0L) {
    s_beta <- 0
  }

  z <- y - bet * v * sqrt(w)
  B <- matrix(1, nrow = N)
  Vinv <- Matrix::Diagonal(x = 1 / v)
  U <- solve(symchol(crossprod(B, Vinv) %*% B + scale / tau * Diagonal(M + 1L)))
  U2 <- crossprod(z / v, B) %*% U

  bias <- sqrt(w) * bet * mu_v
  s2 <- scale * w * mu_v + w * bet^2 * s2_v

  cnt1 <- cnt2 <- rep(0, 3)
  cntw <- 0
  adapt_iter <- rep(FALSE, nmcmc)

  if (w_prior$type == "GBP") {
    if (is.null(w_prior$prop_sigma)) {
      w_prior$prop_sigma <- min(max(var(y) / scale / 2, 1e-3), 5)
    }
    if (is.null(w_prior$adapt)) {
      w_prior$adapt <- TRUE
    }

    if (w_prior$adapt) {
      if (is.null(w_prior$adapt_delay)) {
        w_prior$adapt_delay <- floor(nburn / 10)
      }
      if (is.null(w_prior$adapt_thin)) {
        w_prior$adapt_thin <- 1L
      }

      w_prior$count <- 1L
      w_prior$sum1 <- log(w)
      w_prior$welford <- 0
      adapt_iter <- seq_len(nmcmc) %in% seq.int(1L, nburn, by = w_prior$adapt_thin)
    }
  }

  tX <- t(X)

  if (verbose) {
    pr <- c("MCMC iteration", 0, myTimestamp(), "nbasis:", M)
    cat(pr, "\n")
  }

  # ---------------------------------------------------------------------------
  # MCMC loop
  # ---------------------------------------------------------------------------
  for (k in seq_len(nmcmc)) {
    if (k == lag_beta) {
      s_beta <- s_beta_hold
    }

    move <- move_type(M, maxBasis, moveProbs)

    # -------------------------------------------------------------------------
    # Birth step
    # -------------------------------------------------------------------------
    if (move == "B") {
      J_cand <- sample(maxInt, 1, prob = Iw0)
      u_cand <- sample(p, J_cand, replace = FALSE, prob = Zw0)
      s_cand <- 1 - 2 * rbinom(J_cand, 1, 0.5)
      t_cand <- runif(J_cand)

      B_new <- makeBasis(s_cand, u_cand, t_cand, tX, 1)

      if (sum(B_new > 0) >= npart) {
        B_cand <- cbind(B, B_new)

        U_cand <- tryCatch(
          solve(symchol(crossprod(B_cand, Vinv) %*% B_cand + scale / tau * Diagonal(M + 2L))),
          error = function(e) FALSE
        )

        if (isFALSE(U_cand)) {
          log_accept_prob_B <- -Inf
        } else {
          U2_cand <- crossprod(z / v, B_cand) %*% U_cand
          log_accept_prob_B <-
            sum(log(abs(diag(U_cand)))) - sum(log(abs(diag(U)))) +
            1 / (2 * w * scale) * (TCP(U2_cand) - TCP(U2))

          log_accept_prob_B <-
            log_accept_prob_B +
            log(scale) / 2 - log(tau) / 2 + log(lam) +
            log(moveProbs[2]) - log(moveProbs[1]) -
            log(M + 1) - log(maxInt) - lchoose(p, J_cand) -
            log(Iw0[J_cand] / sum(Iw0)) -
            log(dmwnchBass(Zw0, u_cand))
        }
      } else {
        log_accept_prob_B <- -Inf
      }

      if (is.nan(as.numeric(log_accept_prob_B))) {
        browser()
      }

      cnt1[1] <- cnt1[1] + 1
      if (log(runif(1)) < as.numeric(log_accept_prob_B)) {
        cnt2[1] <- cnt2[1] + 1
        B <- B_cand
        U <- U_cand
        U2 <- U2_cand
        Iw0[J_cand] <- Iw0[J_cand] + 1
        Zw0[u_cand] <- Zw0[u_cand] + 1
        M <- M + 1L
        lookup[[nlook]] <- list(J = J_cand, s = s_cand, t = t_cand, u = u_cand)
        basis_index <- c(basis_index, nlook)
        nlook <- nlook + 1L
      }
    }

    # -------------------------------------------------------------------------
    # Death step
    # -------------------------------------------------------------------------
    if (move == "D") {
      m0 <- sample(M, 1)
      basis_cand <- lookup[[basis_index[m0]]]
      J_cand <- basis_cand$J
      u_cand <- basis_cand$u

      Zw_cand <- Zw0
      Zw_cand[u_cand] <- Zw_cand[u_cand] - 1
      B_cand <- B[, -(m0 + 1), drop = FALSE]

      U_cand <- tryCatch(
        solve(symchol(crossprod(B_cand, Vinv) %*% B_cand + scale / tau * Diagonal(M))),
        error = function(e) FALSE
      )

      if (isFALSE(U_cand)) {
        log_accept_prob_D <- -Inf
      } else {
        U2_cand <- crossprod(z / v, B_cand) %*% U_cand
        log_accept_prob_D <-
          sum(log(abs(diag(U_cand)))) - sum(log(abs(diag(U)))) +
          1 / (2 * w * scale) * (TCP(U2_cand) - TCP(U2))

        log_accept_prob_D <-
          log_accept_prob_D -
          log(scale) / 2 + log(tau) / 2 - log(lam) -
          log(moveProbs[2]) + log(moveProbs[1]) +
          log(M) + log(maxInt) + lchoose(p, J_cand) +
          log((Iw0[J_cand] - 1) / (sum(Iw0) - 1)) +
          log(dmwnchBass(Zw_cand, u_cand))
      }

      cnt1[2] <- cnt1[2] + 1
      if (is.nan(as.numeric(log_accept_prob_D))) {
        browser()
      }

      if (log(runif(1)) < as.numeric(log_accept_prob_D)) {
        cnt2[2] <- cnt2[2] + 1
        B <- B_cand
        U <- U_cand
        U2 <- U2_cand
        Iw0[J_cand] <- Iw0[J_cand] - 1
        Zw0 <- Zw_cand
        M <- M - 1L
        basis_index <- basis_index[-m0]
      }
    }

    # -------------------------------------------------------------------------
    # Mutation step
    # -------------------------------------------------------------------------
    if (move == "M") {
      m0 <- sample(M, 1)
      basis_cand <- lookup[[basis_index[m0]]]
      J_cand <- basis_cand$J
      u_cand <- basis_cand$u
      s_cand <- 1 - 2 * rbinom(J_cand, 1, 0.5)
      t_cand <- runif(J_cand)

      B_new <- makeBasis(s_cand, u_cand, t_cand, tX, 1)

      if (sum(B_new > 0) > npart) {
        B_cand <- B
        B_cand[, m0 + 1] <- B_new

        U_cand <- tryCatch(
          solve(symchol(crossprod(B_cand, Vinv) %*% B_cand + scale / tau * Diagonal(M + 1L))),
          error = function(e) FALSE
        )

        if (isFALSE(U_cand)) {
          log_accept_prob_M <- -Inf
        } else {
          U2_cand <- crossprod(z / v, B_cand) %*% U_cand
          log_accept_prob_M <-
            sum(log(abs(diag(U_cand)))) - sum(log(abs(diag(U)))) +
            1 / (2 * w * scale) * (sum(U2_cand^2) - sum(U2^2))
        }
      } else {
        log_accept_prob_M <- -Inf
      }

      if (is.nan(as.numeric(log_accept_prob_M))) {
        browser()
      }

      cnt1[3] <- cnt1[3] + 1
      if (log(runif(1)) < as.numeric(log_accept_prob_M)) {
        cnt2[3] <- cnt2[3] + 1
        B <- B_cand
        U <- U_cand
        U2 <- U2_cand
        lookup[[nlook]] <- list(J = J_cand, s = s_cand, t = t_cand, u = u_cand)
        basis_index[m0] <- nlook
        nlook <- nlook + 1L
      }
    }

    # -------------------------------------------------------------------------
    # Gibbs / Metropolis updates
    # -------------------------------------------------------------------------
    a <- TCP(U, U2) + sqrt(scale * w) * (U %*% rnorm(M + 1L))
    yhat <- B %*% a
    r <- y - yhat
    pen <- sum(a^2)

    lam <- rgamma(1, a_lambda + M, b_lambda + 1)
    tau <- 1 / rgamma(1, a_tau + (M + 1L) / 2, b_tau + pen / (2 * w))

    bet <- rnorm(
      1,
      (s_beta^2 * sum(r) / sqrt(w) + scale * m_beta) / (s_beta^2 * sum(v) + scale),
      sqrt(scale * s_beta^2 / (s_beta^2 * sum(v) + scale))
    )

    gam <- rnorm(
      1,
      (s_gamma^2 * N + m_gamma) / (s_gamma^2 * sum(v) + 1),
      s_gamma / sqrt(s_gamma^2 * sum(v) + 1)
    )

    # Update w
    rss_over_v <- sum(r^2 / v) / scale
    sum_r <- sum(r)
    quad_term <- w_prior$b + rss_over_v + pen / tau
    shape_shift <- (N + M + 1L) / 2

    if (w_prior$type == "GIG") {
      if (abs(bet) < 1e-9) {
        w_cand <- rgig2(
          p = w_prior$p - shape_shift,
          a = w_prior$a,
          b = quad_term
        )

        if (w_cand >= w_prior$lower_bound) {
          w <- w_cand
        }
      } else {
        w_tform <- rMHN(
          n = 1,
          alpha = N + M - 2 * w_prior$p + 1,
          beta = 0.5 * quad_term,
          gamma = bet * sum_r / scale
        )

        w <- w_tform^(-2)
      }
    } else {
      w_cand <- exp(log(w) + rnorm(1, 0, w_prior$prop_sigma))

      log_alpha_w <-
        (w_prior$p * w_prior$a - shape_shift) * (log(w_cand) - log(w)) -
        0.5 * (rss_over_v + pen / tau) * (1 / w_cand - 1 / w) -
        (w_prior$a + w_prior$b) *
        (log(1 + w_cand^w_prior$p) - log(1 + w^w_prior$p)) +
        (bet / scale) * sum_r * (1 / sqrt(w_cand) - 1 / sqrt(w)) +
        log(w_cand > w_prior$lower_bound)

      if (log_alpha_w > log(runif(1))) {
        cntw <- cntw + 1
        w <- w_cand
      }

      if (w_prior$adapt) {
        if (adapt_iter[k]) {
          theta <- log(w)

          n_adapt <- w_prior$count + 1L
          S1_adapt <- w_prior$sum1 + theta

          w_prior$count <- n_adapt
          w_prior$sum1 <- S1_adapt

          xbar_n <- S1_adapt / n_adapt
          xbar_n_1 <- (S1_adapt - theta) / (n_adapt - 1L)
          w_prior$welford <- w_prior$welford + (theta - xbar_n_1) * (theta - xbar_n)

          if (k > w_prior$adapt_delay) {
            var_adapt <- max(w_prior$welford / (n_adapt - 1L), 0)
            ps_adapt <- sqrt(2.38^2 * var_adapt + 1e-8)
            w_prior$prop_sigma <- min(max(ps_adapt, 1e-3), 5)
          }
        }
      }
    }

    # Update latent local scales v
    v <- rgig2.vec(
      p = v_prior$p - 1 / 2,
      a = v_prior$a + bet^2 / scale,
      b = as.numeric(v_prior$b + r^2 / (w * scale))
    )

    z <- y - bet * v * sqrt(w)
    Vinv <- Matrix::Diagonal(x = 1 / v)
    U <- solve(symchol(crossprod(B, Vinv) %*% B + scale / tau * Diagonal(M + 1L)))
    U2 <- crossprod(z / v, B) %*% U

    v_prior$a <- gam^2
    mu_v <- mu_gig(v_prior$p, v_prior$a, v_prior$b)
    s2_v <- var_gig(v_prior$p, v_prior$a, v_prior$b)

    bias <- sqrt(w) * bet * mu_v
    s2 <- scale * w * mu_v + w * bet^2 * s2_v

    if (keep_iter[k]) {
      ss[kk] <- mean((y - yhat)^2)
      M_mc[kk] <- M
      lam_mc[kk] <- lam
      tau_mc[kk] <- tau
      w_mc[kk] <- w
      bet_mc[kk] <- bet
      gam_mc[kk] <- gam
      v_mc[kk, ] <- v
      basis_mc[[kk]] <- basis_index
      a_mc[[kk]] <- a
      bias_mc[kk] <- bias
      s2_mc[kk] <- s2
      kk <- kk + 1L
    }

    if (verbose && k %% 1000 == 0) {
      pr <- c("MCMC iteration", k, myTimestamp(), "nbasis:", M)
      cat(pr, "\n")
    }
  }

  obj <- list(
    nbasis = M_mc,
    M = M_mc,
    w = w_mc,
    v = v_mc,
    tau = tau_mc,
    lam = lam_mc,
    lamb = lam_mc,
    a = a_mc,
    beta = bet_mc,
    gamma = gam_mc,
    basis = basis_mc,
    lookup = lookup,
    cnt1 = cnt1,
    cnt2 = cnt2,
    cntw = cntw,
    ss = ss,
    v_prior = v_prior,
    X = X,
    y = y,
    scale = scale,
    s2 = s2_mc,
    bias = bias_mc
  )

  class(obj) <- c("nwbass", "gbass")
  obj
}


#' @export
#' @param ... (for backwards compatability)
#' @rdname nwbass
nwbass2 <- function(...) {
  .Deprecated("nwbass")
  nwbass(...)
}

