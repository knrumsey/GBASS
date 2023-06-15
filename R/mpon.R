#' Random generator for the "monomial perturbation of normal" distribution
#'
#' A function to generate from the mpon distribution. Uses a uniformly bounded rejection sampling scheme, adapted from Devroye (2014).
#'
#' @param n currently, n = 1. Changing this argument has no effect on the code.
#' @param alpha parameter
#' @param gamma parameter
#' @param mu parameter
#' @param max_steps maximum number of steps for Newton Raphson
#' @param tol tolerance for Newton Raphson
#' @return a random draw from the mpon distribution
#' @details Generates a RV with (unnormalized) density function f(x) = x^(alpha-1)exp(-gamma*(x-mu)^2)
#' @export
#' @examples
#' n <- 10000
#' alpha <- 5
#' gamma <- 1
#' mu <- 2
#' y <- rep(NA, n)
#' for(i in 1:n){
#'   y[i] <- rmpon(1, alpha, gamma, mu)
#' }
#' hist(y, breaks=30, freq=F)
#' c <- integrate(dmpon, lower=0, upper=20, alpha=alpha, gamma=gamma, mu=mu)$value
#' curve(dmpon(x, alpha, gamma, mu)/c, add=TRUE, lwd=2, col='blue')
rmpon <- function(n=1, alpha, gamma, mu, max_steps=10000, tol=1e-4){
  m <- mu/2 + 1/(2*gamma)*sqrt(gamma*(2*alpha+gamma*mu^2-2))
  lf <- function(xx) dmpon(xx, alpha, gamma, mu, log=TRUE)
  phi <- function(xx) exp((alpha-1)*(log(xx+m)-log(m)) - gamma*((xx+m-mu)^2 - (m-mu)^2))
  psi <- function(xx) (alpha-1)*log(xx+m) - gamma*(xx+m-mu)^2 - lf(m)
  psip <- function(xx)  (alpha-1)/(xx+m) - 2*gamma*(xx+m-mu)
  #Newton-Raphson
  s <- t <- 0.5
  #Get decent starting values
  if(-psi(t) < 1){
    while(-psi(t) < 0.1){
      t <- 3*t
      s <- 3*s
    }
  }else{
    while(-psi(t) > 10){
      t <- t/3
      s <- s/3
    }
  }
  t_flag <- s_flag <- TRUE
  cnt <- 1
  while((t_flag | s_flag) & cnt <= max_steps){
    t <- t - (psi(t)+1)/psip(t)
    s <- s + (psi(-s)+1)/psip(-s)

    if(abs(psi(t) + 1) < tol) t_flag <- FALSE
    if(abs(psi(-s) + 1) < tol) s_flag <- FALSE
    cnt <- cnt + 1
  }
  if(cnt > max_steps){
    warning("cnt reached max_steps -- may not have converged. psi(t) = ",
            round(psi(t), 7), " and psi(-s) = ", round(psi(-s),7),
            "\n\talpha = ", alpha,
            "\n\tgamma = ", gamma,
            "\n\tmu = ", mu)
  }

  if(min(s,t) < 0){
    #browser()
    if(gamma/alpha > 1e10){
      warning(paste0("Error: s and/or t is negative: alpha=", alpha, ", gamma=", gamma, ", mu=", mu), "\n gamma/alpha is large, so we return mu")
      return(mu)
    }else{
      if(alpha/gamma > 1e10){
        warning(paste0("Error: s and/or t is negative: alpha=", alpha, ", gamma=", gamma, ", mu=", mu), "\n alpha/gamma is large, so we return 1")
        return(1)
      }else{
        stop(paste0("Error: s and/or t is negative: alpha=", alpha, ", gamma=", gamma, ", mu=", mu), "\n cannot determine limiting case. Fail to return.")
      }
    }

  }
  #Devroye algorithm
  p <- 1/psip(-s)
  r <- -1/psip(t)
  tp <- t + r*psi(t)
  sp <- s + p*psi(-s)
  q <- tp + sp
  flag <- TRUE

  while(flag){
    U <- runif(1); V <- runif(1); W <- runif(1)
    if(U < q/(q+r+p)){
      X <- -sp + q*V
    }else{
      if(U < (q+r)/(q+r+p)){
        X <- tp - r*log(V)
      }else{
        X <- -sp + p*log(V)
      }
    }

    if(X > -m){
      if(X > tp){
        chi <- exp(psi(t) + psip(t)*(X-t))
      }else{
        if(X > -sp){
          chi <- 1
        }else{
          chi <- exp(psi(-s) + psip(-s)*(X+s))
        }
      }
      if(log(W) <= psi(X) - log(chi)){
        flag <- FALSE
      }
    }
  }
  return(X+m)
}


#' Un-normalized density of the "monomial perturbation of normal (mpon)" distribution
#'
#' A function to evaluate the unnormalized density of the mpon distribution.
#'
#' @param x locations for density evaluation
#' @param alpha parameter
#' @param gamma parameter
#' @param mu parameter
#' @return density of the mpon distribution
#' @details Evaluates the (unnormalized) density function f(x) = x^(alpha-1)exp(-gamma*(x-mu)^2)
#' @export
#' @examples
#' n <- 10000
#' alpha <- 5
#' gamma <- 1
#' mu <- 2
#' y <- rep(NA, n)
#' for(i in 1:n){
#'   y[i] <- rmpon(1, alpha, gamma, mu)
#' }
#' hist(y, breaks=30, freq=F)
#' c <- integrate(dmpon, lower=0, upper=20, alpha=alpha, gamma=gamma, mu=mu)$value
#' curve(dmpon(x, alpha, gamma, mu)/c, add=TRUE, lwd=2, col='blue')
dmpon <- function(x, alpha, gamma, mu, log=FALSE){
  res <- (alpha-1)*log(x) - gamma*(x-mu)^2 + log(x > 0)
  if(log){
    return(res)
  }else{
    return(exp(res))
  }
}



#' Random generator for the "monomial perturbation of normal" distribution
#'
#' A function to generate from the mpon distribution. This is an alternative approach based on Sun et al (2021). Assumes alpha > 1.
#' (depreciated? See modified half normal distribution)
#'
#' @param n currently, n = 1. Changing this argument has no effect on the code.
#' @param alpha parameter
#' @param gamma parameter
#' @param mu parameter
#' @param max_steps maximum number of steps for rejection sampler. Returns error if met.
#' @return a random draw from the mpon distribution
#' @details Generates a RV with (unnormalized) density function f(x) = x^(alpha-1)exp(-gamma*(x-mu)^2)
#' @export
#' @examples
#' n <- 10000
#' alpha <- 5
#' gamma <- 1
#' mu <- 2
#' y <- rep(NA, n)
#' for(i in 1:n){
#'   y[i] <- rmpon(1, alpha, gamma, mu)
#' }
#' hist(y, breaks=30, freq=F)
#' c <- integrate(dmpon, lower=0, upper=20, alpha=alpha, gamma=gamma, mu=mu)$value
#' curve(dmpon(x, alpha, gamma, mu)/c, add=TRUE, lwd=2, col='blue')
rmpon_sun <- function(n=1, alpha, gamma, mu, max_steps=1e6){
  # Reparameterize in terms of Sun et al. (2021)
  beta <- gamma
  gamma <- 2*mu*beta

  # Calculate optimal choices of enevelope parameters (Thm 1)
  mu_opt    <- (gamma + sqrt(gamma^2 + 8*(alpha - 1)*beta))/(4*beta)
  delta_opt <- max(beta + (gamma^2 - gamma*sqrt(gamma^2 + 8*alpha*beta))/(4*alpha), 0)

  # Compute K1 and K2 (Thm 1)
  ka1 <- log(2*sqrt(pi)) + (alpha-1)*log(sqrt(beta)*(alpha-1)/(2*beta*mu_opt-gamma))
  kb1 <- -(alpha-1)+beta*mu_opt^2
  K1  <- ka1 + kb1
  ka2 <- (alpha/2)*log(beta) + lgamma(alpha/2) - alpha/2*log(delta_opt)
  kb2 <- (gamma^2/(4*(beta-delta_opt)))
  K2  <- ka2 + kb2

  # Sample with algorithm 1
  # Case one (K1 < K2)
  if(K1 < K2){
    cnt <- 1
    while(TRUE){
      cnt <- cnt + 1
      #xx <- rnorm(1, mu_opt, sqrt(1/(2*beta))) # INCONSISTENCY in the Sun et al paper
      xx <- rnorm(1, mu_opt, 1/(2*beta))
      uu <- log(runif(1))
      #aa <- (alpha-1)*log(xx) - log(mu_opt) + (2*beta*mu_opt-gamma)*(mu_opt-xx)
      aa <- (alpha-1)*log(xx*(2*beta*mu_opt - gamma)/(alpha-1)) +
            xx*(gamma - 2*mu_opt*beta) + alpha - 1 # The "corrected" acceptance probability
      if(xx > 0 & uu < aa){
        return(xx)
      }
      if(cnt > max_steps){
        stop("More than 1e9 iterations needed
             to obtain a sample from modified half normal\n
             \n\talpha: ", alpha,
             "\n\tgamma: ", beta,
             "\n\tmu: ", mu)
      }
    }
  }
  # Case two (K1 >= K2)
  else{
    cnt <- 1
    while(TRUE){
      cnt <- cnt + 1
      xx <- sqrt(rgamma(1, alpha/2, delta_opt))
      uu <- log(runif(1))
      aa <- -(beta-delta_opt)*xx^2 + gamma*xx - gamma^2/(4*(beta-delta_opt))
      if(xx > 0 & uu < aa){
        return(xx)
      }
      if(cnt > max_steps){
        stop("More than 1e9 iterations needed
             to obtain a sample from modified half normal\n
             \n\talpha: ", alpha,
             "\n\tgamma: ", beta,
             "\n\tmu: ", mu)
      }
    }
  }
}

