#' Random generation from the Modified Half-Normal distribution
#'
#' Generates random samples from the modified half-normal distribution with density
#' proportional to x^(alpha - 1) * exp(-beta * x^2 + gamma * x) for x > 0.
#'
#' @param n Number of samples.
#' @param alpha Shape parameter. Must satisfy alpha >= 1.
#' @param beta Positive rate parameter on x^2.
#' @param gamma Real-valued linear parameter.
#'
#' @return A numeric vector of length \code{n}.
#' @export
rMHN <- function(n = 1, alpha, beta, gamma) {
  .mhn_validate_args(n = n, alpha = alpha, beta = beta, gamma = gamma)

  if (n == 0L) {
    return(numeric(0))
  }

  if (gamma <= 0) {
    return(.mhn_sample_negative_gamma(
      n = n, alpha = alpha, beta = beta, gamma = gamma
    ))
  }

  log_k1 <- .mhn_log_k1(alpha = alpha, beta = beta, gamma = gamma)
  log_k2 <- .mhn_log_k2(alpha = alpha, beta = beta, gamma = gamma)

  if (log_k2 <= log_k1) {
    .mhn_sample_positive_gamma_sqrt_gamma(
      n = n, alpha = alpha, beta = beta, gamma = gamma
    )
  } else {
    .mhn_sample_positive_gamma_normal(
      n = n, alpha = alpha, beta = beta, gamma = gamma
    )
  }
}

.mhn_validate_args <- function(n, alpha, beta, gamma) {
  if (!is.numeric(n) || length(n) != 1L || is.na(n) || n < 0 || n != as.integer(n)) {
    stop("n must be a nonnegative integer.")
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) || alpha < 1) {
    stop("alpha must be a scalar with alpha >= 1.")
  }
  if (!is.numeric(beta) || length(beta) != 1L || is.na(beta) || beta <= 0) {
    stop("beta must be a positive scalar.")
  }
  if (!is.numeric(gamma) || length(gamma) != 1L || is.na(gamma)) {
    stop("gamma must be a finite scalar.")
  }
  invisible(NULL)
}

.mhn_rejection_fill <- function(n, draw_fn, accept_fn) {
  out <- numeric(0)

  while (length(out) < n) {
    need <- n - length(out)
    proposal <- draw_fn(need)
    keep <- accept_fn(proposal)
    out <- c(out, proposal[keep])
  }

  out
}

.mhn_sample_truncated_normal_positive <- function(n, mean, sd) {
  out <- numeric(0)

  while (length(out) < n) {
    need <- n - length(out)
    x <- stats::rnorm(need, mean = mean, sd = sd)
    out <- c(out, x[x > 0])
  }

  out
}

.mhn_log_k1 <- function(alpha, beta, gamma) {
  mu <- (gamma + sqrt(gamma^2 + 8 * beta * (alpha - 1))) / (4 * beta)

  log(2 * sqrt(pi)) +
    (alpha - 1) * (log(sqrt(beta) * (alpha - 1)) - log(2 * beta * mu - gamma)) -
    (alpha - 1) +
    beta * mu^2
}

.mhn_log_k2 <- function(alpha, beta, gamma) {
  delta <- beta + (gamma^2 - gamma * sqrt(gamma^2 + 8 * alpha * beta)) / (4 * alpha)

  (alpha / 2) * log(beta) +
    lgamma(alpha / 2) +
    gamma^2 / (4 * (beta - delta)) -
    (alpha / 2) * log(delta)
}

.mhn_sample_positive_gamma_sqrt_gamma <- function(n, alpha, beta, gamma) {
  delta <- beta + (gamma^2 - sqrt(gamma^4 + 8 * alpha * beta * gamma^2)) / (4 * alpha)

  .mhn_rejection_fill(
    n = n,
    draw_fn = function(m) sqrt(stats::rgamma(m, shape = alpha / 2, rate = delta)),
    accept_fn = function(x) {
      log(stats::runif(length(x))) <=
        (-(beta - delta) * x^2 + gamma * x - gamma^2 / (4 * (beta - delta)))
    }
  )
}

.mhn_sample_positive_gamma_normal <- function(n, alpha, beta, gamma) {
  mode <- (gamma + sqrt(gamma^2 + 8 * beta * (alpha - 1))) / (4 * beta)
  sd <- 1 / sqrt(2 * beta)

  .mhn_rejection_fill(
    n = n,
    draw_fn = function(m) .mhn_sample_truncated_normal_positive(m, mean = mode, sd = sd),
    accept_fn = function(x) {
      log(stats::runif(length(x))) <=
        ((alpha - 1) * log(x / mode) + (2 * beta * mode - gamma) * (mode - x))
    }
  )
}

.mhn_sample_negative_gamma <- function(n, alpha, beta, gamma) {
  if (gamma > 0) {
    stop("Internal error: gamma must satisfy gamma <= 0.")
  }

  a <- beta
  b <- abs(gamma)

  if (alpha <= 1) {
    stop("Internal error: alpha < 1 is unsupported in this version of rMHN.")
  }

  x_mode <- (sqrt(b^2 + 8 * (alpha - 1) * a) - b) / (4 * a)
  scale <- x_mode

  new_a <- a * scale^2
  new_b <- b * scale

  shape_par <- ((new_a + new_b) / (2 * new_a + new_b)) * alpha
  rate_par <- new_a + new_b
  exponent_power <- (new_a + new_b) / (2 * new_a + new_b)

  x_scaled <- .mhn_rejection_fill(
    n = n,
    draw_fn = function(m) stats::rgamma(m, shape = shape_par, rate = rate_par)^exponent_power,
    accept_fn = function(x) {
      exp_g <- x^(1 / exponent_power)
      log(stats::runif(length(x))) <=
        (-new_a * x^2 - new_b * x + (new_a + new_b) * exp_g)
    }
  )

  x_scaled * scale
}
