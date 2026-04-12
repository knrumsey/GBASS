#' Function to select prior for gamma
#'
#' A  function to select hyperparameters for gamma in terms of steepness parameter xi = (1+gamma)^(-1/2)
#'
#' @param q1 default 0.1
#' @param q2 default 0.9
#' @param p1 default 0.5
#' @param p2 default 0.05
#' @param par0 optional starting values
#' @param par0 optional  ridge penalty for optimization
#' @return hyperparameter values for prior of gamma
#' @details Plots the posterior draws of the steepness and asymmetry parameter, defined as (1+gamma)^(-1/2) and beta/sqrt(gamma^2+beta^2)*steepness respectively.
#' These parameters are location and scale invariant. Normal (0,0) and Cauchy (0,1) occur as limiting cases.
#' @import Matrix
#' @export
#' @examples
#' #not yet
nw_gamma_prior <- function(q1=0.1, q2=0.9, p1=0.5, p2=0.05, par0=NULL, lambda=0){
  logit <- function(z) log(z/(1-z))
  fff <- function(xx, q1,q2,p1,p2,lambda){
    lhs1 <- 1 + pnorm((1-q1^-2-xx[1])/xx[2]) - pnorm((q1^-2-1-xx[1])/xx[2])
    lhs2 <- pnorm((q2^-2-1-xx[1])/xx[2]) - pnorm((1-q2^-2-xx[1])/xx[2])
    (logit(lhs1)-logit(p1))^2 + (logit(lhs2)-logit(p2))^2 + lambda*sum(xx^2)
  }
  if(is.null(par0)) par0 <- c(100, 100)
  opt <- optim(par0, fff, method="Nelder-Mead", q1=q1, q2=q2, p1=p1, p2=p2, lambda=lambda, control=list(maxit=c(50000)))
  par <- opt$par
  #return(opt)
  #if(opt$value - lambda*sum(par^2) > 0.5 ) warning("Attained a minimum of ", opt$value - lambda*sum(par^2), ", optimization may not have converged")
  names(par) <- c("m_gamma", "s_gamma")
  return(par)
}


#' Method of Moments for Normal-Wald
#'
#' A function to estimate parameters of the Normal-Wald distribution
#'
#' @param data the data vector
#' @param stats a vector of length 4 containing mean, variance, skew, kurtosis. Ignored when data is provided.
#' @param mu location parameter, if fixed
#' @param alpha tail heaviness parameter, if fixed
#' @param beta skewness parameter, if fixed
#' @param delta scale parameter, if fixed
#' @param triangle logical. When TRUE, only the steepness and asymmetry values are returned.
#' @param ... additional parameters passed to optim.
#' @return estimated parameters
#' @details Method of moments estimators for NW parameters. The stats vector can contain NA values when parameters are fixed. If the mean is to be estimated, then \code{stats[1]} must be provided.
#' @import Matrix
#' @export
#' @examples
#' n <- 500
#' y <- rgamma(n, 3, 1.5) + rlnorm(n, 1, 0.5)
#' skew <- mean(((y-mean(y))/sd(y))^3) # Sample skewness
#' kurt <- mean(((y-mean(y))/sd(y))^4) # Sample kurtosis
#'
#' nw_est_mom(stats=c(NA, NA, skew, kurt), mu=0, delta=1, triangle=TRUE)
#' nw_est_mom(stats=c(NA, var(y), skew, kurt), mu=0, triangle=TRUE)
nw_est_mom <- function(data=NULL, stats=NULL,
                       mu=NA, delta=NA, beta=NA, alpha=NA,
                       triangle=FALSE, ...){
  if(!is.null(data)){
    zdata <- (data-mean(data))/sd(data)
    stats <- c(mean(data), var(data), mean(zdata^3), mean(zdata^4))
  }
  # Figure out which moments to use
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
    if(!is.na(moments[1])){
      d1 <- (m1 - moments[1])^2
    }
    if(!is.na(moments[2])){
      d2 <- (m2 - moments[2])^2
    }
    if(!is.na(moments[3])){
      d3 <- (m3 - moments[3])^2
    }
    if(!is.na(moments[4])){
      d4 <- (m4 - moments[4])^2
    }
    return(d1 + d2 + d3 + d4)
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
    gamma <-  sqrt(alpha^2 - beta^2)
    steep <- 1/sqrt(1 + abs(gamma))
    asymm <- beta*steep/alpha
    res <- c(asymm, steep)
    names(res) <- c("asymmetry", "steepness")
  }else{
    res <- c(mu, delta, beta, alpha)
    names(res) <- c("mu", "delta", "beta", "alpha")
  }
  return(res)
}
