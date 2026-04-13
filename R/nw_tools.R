#' Calibrate a prior for the Normal-Wald gamma parameter
#'
#' Selects prior hyperparameters for \code{gamma} in the Normal-Wald model using
#' probability statements about the steepness parameter
#' \eqn{\xi = (1 + \gamma)^{-1/2}}.
#'
#' @param q1 Lower reference value for the steepness parameter. Default is
#'   \code{0.1}.
#' @param q2 Upper reference value for the steepness parameter. Default is
#'   \code{0.9}.
#' @param p1 Target probability associated with \code{q1}. Default is
#'   \code{0.5}.
#' @param p2 Target probability associated with \code{q2}. Default is
#'   \code{0.05}.
#' @param par0 Optional starting values for the numerical optimization.
#' @param lambda Optional ridge penalty used in the optimization. Default is
#'   \code{0}.
#'
#' @return A named numeric vector with entries \code{m_gamma} and
#'   \code{s_gamma}.
#'
#' @details
#' This function chooses hyperparameters for a prior on \code{gamma} by solving
#' a simple calibration problem based on the steepness parameter
#' \eqn{(1 + \gamma)^{-1/2}}. The optimization is carried out with
#' \code{\link[stats]{optim}} using the Nelder-Mead method.
#'
#' @export
#'
#' @examples
#' nw_gamma_prior()
#' nw_gamma_prior(q1 = 0.2, q2 = 0.8, p1 = 0.4, p2 = 0.1)
nw_gamma_prior <- function(q1=0.1, q2=0.9, p1=0.5, p2=0.05, par0=NULL, lambda=0){
  logit <- function(z) log(z/(1-z))
  fff <- function(xx, q1,q2,p1,p2,lambda){
    lhs1 <- 1 + pnorm((1-q1^-2-xx[1])/xx[2]) - pnorm((q1^-2-1-xx[1])/xx[2])
    lhs2 <- pnorm((q2^-2-1-xx[1])/xx[2]) - pnorm((1-q2^-2-1-xx[1])/xx[2])
    (logit(lhs1)-logit(p1))^2 + (logit(lhs2)-logit(p2))^2 + lambda*sum(xx^2)
  }
  if(is.null(par0)) par0 <- c(100, 100)
  opt <- optim(par0, fff, method="Nelder-Mead", q1=q1, q2=q2, p1=p1, p2=p2, lambda=lambda, control=list(maxit=c(50000)))
  par <- opt$par
  names(par) <- c("m_gamma", "s_gamma")
  return(par)
}

#' Method-of-moments estimation for the Normal-Wald distribution
#'
#' Estimates Normal-Wald parameters from sample moments or from a supplied
#' vector of moments.
#'
#' @param data Optional data vector. If supplied, empirical moments are computed
#'   from \code{data} and \code{stats} is ignored.
#' @param stats Optional numeric vector of length 4 containing mean, variance,
#'   skewness, and kurtosis. Entries may be \code{NA} when the corresponding
#'   parameter is fixed.
#' @param mu Location parameter, if fixed.
#' @param delta Scale parameter, if fixed.
#' @param beta Skewness parameter, if fixed.
#' @param alpha Tail parameter, if fixed.
#' @param triangle Logical; if \code{TRUE}, return only the asymmetry and
#'   steepness summary parameters.
#' @param ... Additional arguments passed to \code{\link[stats]{optim}}.
#'
#' @return
#' If \code{triangle = FALSE}, a named numeric vector with entries
#' \code{mu}, \code{delta}, \code{beta}, and \code{alpha}.
#' If \code{triangle = TRUE}, a named numeric vector with entries
#' \code{asymmetry} and \code{steepness}.
#'
#' @details
#' This function computes method-of-moments estimates for the Normal-Wald
#' distribution. If \code{data} is supplied, sample mean, variance, skewness,
#' and kurtosis are computed automatically. Otherwise, the user may provide
#' these moments directly through \code{stats}.
#'
#' If some parameters are fixed, only the remaining parameters are estimated.
#' When \code{triangle = TRUE}, the returned values are the transformed summary
#' quantities
#' \eqn{\text{steepness} = (1 + |\gamma|)^{-1/2}} and a corresponding
#' asymmetry measure.
#'
#' @export
#'
#' @examples
#' n <- 500
#' y <- rgamma(n, 3, 1.5) + rlnorm(n, 1, 0.5)
#' z <- (y - mean(y)) / sd(y)
#' skew <- mean(z^3)
#' kurt <- mean(z^4)
#'
#' nw_est_mom(stats = c(NA, NA, skew, kurt), mu = 0, delta = 1, triangle = TRUE)
#' nw_est_mom(stats = c(NA, var(y), skew, kurt), mu = 0, triangle = TRUE)
nw_est_mom <- function(data=NULL, stats=NULL,
                       mu=NA, delta=NA, beta=NA, alpha=NA,
                       triangle=FALSE, ...){
  if(!is.null(data)){
    zdata <- (data-mean(data))/sd(data)
    stats <- c(mean(data), var(data), mean(zdata^3), mean(zdata^4))
  }
  moments_requested <- which(c(is.na(mu), is.na(delta), is.na(beta), is.na(alpha)))
  moments_provided  <- which(!is.na(stats))

  if(length(moments_requested) > length(moments_provided)){
    stop("Need more moments to estimate desired parameters")
  }
  if(moments_requested[1] == 1 & moments_provided[1] != 1){
    stop("Cannot estimate mu without sample mean")
  }

  f2opt <- function(par, moments, mu=NA, delta=NA, beta=NA, alpha=NA){
    cnt <- 1
    if(is.na(mu)){
      mu <- par[cnt]
      cnt <- cnt + 1
    }
    if(is.na(delta)){
      delta <- exp(par[cnt])
      cnt <- cnt + 1
    }
    if(is.na(beta)){
      beta <- par[cnt]
      cnt <- cnt + 1
    }
    if(is.na(alpha)){
      alpha <- par[cnt]
    }
    if(alpha^2 < beta^2) return(1e7)
    gamma <- sqrt(alpha^2 - beta^2)

    m1 <- mu + delta*beta/gamma
    m2 <- delta*alpha^2/gamma^3
    m3 <- 3*beta/(alpha*sqrt(delta*gamma))
    m4 <- 3 + 3*(1 + 4*beta^2)/(alpha^2*delta*gamma)

    d1 <- d2 <- d3 <- d4 <- 0
    if(!is.na(moments[1])) d1 <- (m1 - moments[1])^2
    if(!is.na(moments[2])) d2 <- (m2 - moments[2])^2
    if(!is.na(moments[3])) d3 <- (m3 - moments[3])^2
    if(!is.na(moments[4])) d4 <- (m4 - moments[4])^2

    d1 + d2 + d3 + d4
  }

  par0 <- c(0, 0, 0, 1)[moments_requested]
  fit <- optim(par0, fn=f2opt, moments= stats, mu=mu, delta=delta, beta=beta, alpha=alpha, ...)
  if(fit$convergence != 0) warning("Optim may not have converged.")

  cnt <- 1
  if(1 %in% moments_requested){
    mu <- fit$par[cnt]
    cnt <- cnt + 1
  }
  if(2 %in% moments_requested){
    delta <- exp(fit$par[cnt])
    cnt <- cnt + 1
  }
  if(3 %in% moments_requested){
    beta <- fit$par[cnt]
    cnt <- cnt + 1
  }
  if(4 %in% moments_requested){
    alpha <- abs(fit$par[cnt])
    cnt <- cnt + 1
  }

  if(triangle){
    gamma <- sqrt(alpha^2 - beta^2)
    steep <- 1 / sqrt(1 + abs(gamma))
    asymm <- beta * steep / alpha
    res <- c(asymm, steep)
    names(res) <- c("asymmetry", "steepness")
  } else {
    res <- c(mu, delta, beta, alpha)
    names(res) <- c("mu", "delta", "beta", "alpha")
  }

  return(res)
}


#' Plot Normal-Wald shape parameters on the asymmetry-steepness triangle
#'
#' Visualizes posterior draws of the Normal-Wald shape parameters in terms of
#' asymmetry and steepness. The plotted coordinates are
#' \eqn{\text{steepness} = (1 + \gamma)^{-1/2}} and
#' \eqn{\text{asymmetry} = \beta (\beta^2 + \gamma^2)^{-1/2} \times \text{steepness}}.
#'
#' @param obj A fitted object containing components \code{beta} and \code{gamma},
#'   typically an object returned by \code{nwbass()}.
#' @param add Logical; if \code{FALSE}, initialize a new triangle plot. If
#'   \code{TRUE}, add points to the current plot.
#' @param details Logical; if \code{TRUE} and \code{add = FALSE}, add a few
#'   reference distributions to the plot.
#' @param ... Additional graphical arguments passed to \code{\link{points}}.
#'
#' @return Invisibly returns a data frame with columns \code{asymmetry} and
#'   \code{steepness}.
#'
#' @details
#' This plot provides a simple geometric summary of the shape of the
#' Normal-Wald likelihood. The upper boundary corresponds to symmetric models,
#' while points away from zero on the horizontal axis indicate asymmetry.
#'
#' @export
#'
#' @examples
#' # Suppose mod is a fitted nwbass model
#' # nw_triangle(mod)
nw_triangle <- function(obj, add = FALSE, details = FALSE, ...) {
  if (is.null(obj$beta) || is.null(obj$gamma)) {
    stop("obj must contain components 'beta' and 'gamma'.")
  }

  steepness <- (obj$gamma + 1)^(-1 / 2)
  asymmetry <- obj$beta / sqrt(obj$beta^2 + obj$gamma^2) * steepness

  if (!add) {
    plot(
      NULL,
      xlim = c(-1, 1),
      ylim = c(0, 1),
      xlab = "Asymmetry",
      ylab = "Steepness"
    )

    segments(
      x0 = c(-1, -1, 0),
      x1 = c(1, 0, 1),
      y0 = c(1, 1, 0),
      y1 = c(1, 0, 1),
      lwd = 3,
      col = adjustcolor("gray", alpha.f = 0.5)
    )

    if (details) {
      points(rep(0, 5), c(1, 0.9, 0.5, 0.1428, 0), pch = 0:4, bg = "black")

      legend(
        "bottomright",
        c("Cauchy", "t(5)", "t(10)", "t(100)", "Gaussian",
          "Lnorm(0.5)", "Lnorm(0.1)", "Lnorm(0.01)"),
        pch = c(0:4, 15:17),
        col = 1,
        bty = "n",
        cex = 1.2
      )

      cnt <- 0
      for (a in c(0.5, 0.1, 0.01)) {
        k <- exp(4 * a) + 2 * exp(3 * a) + 3 * exp(2 * a) - 6
        s <- (exp(a) + 2) * sqrt(exp(a) - 1)
        g <- 9 / (3 * k - 4 * s^2)
        p <- s * sqrt(g) / 3

        xi <- (1 + abs(g))^(-1 / 2)
        chi <- p * xi

        points(chi, xi, pch = 15 + cnt)
        cnt <- cnt + 1
      }
    }
  }

  points(asymmetry, steepness, ...)

  invisible(data.frame(asymmetry = asymmetry, steepness = steepness))
}
