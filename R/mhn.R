#' Random generation from the Modified Half-Normal distribution
#'
#' Generates random samples from the modified half-normal distribution with density
#' proportional to x^(alpha - 1) * exp(-beta * x^2 + gamma * x) for x > 0.
#'
#' @param n Number of samples.
#' @param alpha Positive shape parameter.
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

  if (alpha > 1) {
    log_k1 <- .mhn_log_k1(alpha = alpha, beta = beta, gamma = gamma)
    log_k2 <- .mhn_log_k2(alpha = alpha, beta = beta, gamma = gamma)

    if (log_k2 <= log_k1) {
      return(.mhn_sample_positive_gamma_sqrt_gamma(
        n = n, alpha = alpha, beta = beta, gamma = gamma
      ))
    }

    return(.mhn_sample_positive_gamma_normal(
      n = n, alpha = alpha, beta = beta, gamma = gamma
    ))
  }

  .mhn_sample_positive_gamma_general(
    n = n, alpha = alpha, beta = beta, gamma = gamma
  )
}

.mhn_validate_args <- function(n, alpha, beta, gamma) {
  if (!is.numeric(n) || length(n) != 1L || is.na(n) || n < 0 || n != as.integer(n)) {
    stop("n must be a nonnegative integer.")
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) || alpha <= 0) {
    stop("alpha must be a positive scalar.")
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

.mhn_sample_positive_gamma_general <- function(n, alpha, beta, gamma,
                                               tail_tol = 1e-12,
                                               ratio_cap = 0.5,
                                               max_terms = 10000L) {
  weights <- .mhn_positive_gamma_weights(
    alpha = alpha,
    beta = beta,
    gamma = gamma,
    tail_tol = tail_tol,
    ratio_cap = ratio_cap,
    max_terms = max_terms
  )

  w <- sample.int(length(weights), size = n, replace = TRUE, prob = weights) - 1L
  shape <- (alpha + w) / 2
  sqrt(stats::rgamma(n, shape = shape, rate = beta))
}

.mhn_positive_gamma_weights <- function(alpha, beta, gamma,
                                        tail_tol = 1e-12,
                                        ratio_cap = 0.5,
                                        max_terms = 10000L) {
  log_q <- numeric(max_terms)
  log_q[1] <- lgamma(alpha / 2)

  k <- 0L
  ratio <- Inf
  tail_bound <- Inf
  log_q_max <- log_q[1]

  repeat {
    k <- k + 1L
    if (k >= max_terms) {
      stop("Maximum number of mixture terms reached in .mhn_positive_gamma_weights().")
    }

    log_ratio <-
      log(gamma / sqrt(beta)) +
      lgamma((alpha + k) / 2) -
      lgamma((alpha + k - 1) / 2) -
      log(k)

    log_q[k + 1L] <- log_q[k] + log_ratio
    log_q_max <- max(log_q_max, log_q[k + 1L])

    ratio <- exp(log_ratio)

    if (ratio <= ratio_cap) {
      tail_bound <- exp(log_q[k + 1L] - log_q_max) * ratio / (1 - ratio)

      partial_sum <- sum(exp(log_q[seq_len(k + 1L)] - log_q_max))

      if (tail_bound <= tail_tol * partial_sum) {
        break
      }
    }
  }

  q <- exp(log_q[seq_len(k + 1L)] - log_q_max)
  q / sum(q)
}

.mhn_sample_negative_gamma <- function(n, alpha, beta, gamma) {
  if (gamma > 0) {
    stop("Internal error: gamma must satisfy gamma <= 0.")
  }

  a <- beta
  b <- abs(gamma)

  if (alpha <= 1) {
    shape_par <- ((a + b) / (2 * a + b)) * alpha
    rate_par <- a + b
    exponent_power <- (a + b) / (2 * a + b)

    return(.mhn_rejection_fill(
      n = n,
      draw_fn = function(m) stats::rgamma(m, shape = shape_par, rate = rate_par)^exponent_power,
      accept_fn = function(x) {
        x_base <- x^(1 / exponent_power)
        log(stats::runif(length(x))) <=
          (-a * x^2 - b * x + (a + b) * x_base)
      }
    ))
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
      x_base <- x^(1 / exponent_power)
      log(stats::runif(length(x))) <=
        (-new_a * x^2 - new_b * x + (new_a + new_b) * x_base)
    }
  )

  x_scaled * scale
}
