#' @title TBASS - Bayesian MARS with a Student's t likelihood
#'
#' @param X an Nxp matrix of predictor variables.
#' @param y an Nx1 matrix (or vector) of response values.
#' @param df degrees of freedom, default is 5.
#'
#' @details Inverse-gamma prior for local variance components. If df > 2, then
#' the scale is set to (df-2)/df so that the variance of y is equal to w.
#'
#' @export
tbass <- function(X, y, df=5, ...){
  v_prior <- list(type="GIG", p=-df/2, a=0, b=df)
  scale <- ifelse(df > 2, (df-2)/df, 1)
  md=gbass(X, y, v_prior=v_prior, scale=scale, ...)
  md$df <- df
  class(md) <- c("tbass", "gbass")
  return(md)
}


#' @title HBASS - Bayesian MARS with a Horseshoe (or Strawderman-Berger) likelihood
#'
#' @aliases bass_shoe
#'
#' @param X an Nxp matrix of predictor variables.
#' @param y an Nx1 matrix (or vector) of response values.
#' @param likelihood default "h". Also accepts "sb" for Strawderman-Berger
#'
#' @details The Horseshoe and Strawderman-Berger likelihoods can be unified under the class of
#' Generalized Beta Prime distributions. This function provides a wrapper for these likelilhoods. For additional flexibility, use gbass().
#' These distributions have most of their mass near zero but with heavy tails, making this an interesting choice for datasets with potential outliers.
#'
#' @export
hbass <- function(X, y, likelihood="h", ...){
  w_prior <- list(type="GBP", p=1, a=1/2, b=1/2, prop_sigma=0.5)
  if(likelihood=="h"){
    v_prior <- list(type="GBP", p=1, a=1/2, b=1/2, prop_sigma=0.5)
  }else if(likelihood == "sb"){
    v_prior <- list(type="GBP", p=1/2, a=1, b=1/2, prop_sigma=0.5)
  }else{
    stop("likelihood not recognized. See documentation")
  }
  scale <- var(y)
  md=gbass(X, y, w_prior=w_prior, v_prior=v_prior, scale=scale, ...)
  class(md) <- c("hbass", "gbass")
  return(md)
}
bass_shoe <- hbass


#' @title QBASS - Bayesian MARS with an asymmetric Laplace likelhood
#'
#' @param X an Nxp matrix of predictor variables.
#' @param y an Nx1 matrix (or vector) of response values.
#' @param q quantile of interest. Default is 0.5 for median regression.
#' @param prop_sigma_v proposal SD for the local variance factors
#' @param w_prior prior for global variance factor. Default is Jeffreys prior
#'
#' @details Performs quantile regression for quantile q. For many quantiles, consider running in parallel. Can be used to view Sobol decomposition as a function of quantile.
#'
#'
#' @export
qbass <- function(X, y, q=0.5, prop_sigma_v=0.25,
                  w_prior = list(type="GIG", p=0, a=0, b=0, prop_sigma=0.2), ...){
  v_prior <- list(type="GIG", p=1, a=2, b=0, prop_sigma=prop_sigma_v)
  scale <- 2/(q*(1-q))
  m_beta <- 1/q - 1/(1-q)
  s_beta <- 0
  md=gbass(X, y, v_prior=v_prior, w_prior=w_prior, scale=scale, m_beta=m_beta, s_beta=s_beta, ...)
  md$q <- q
  class(md) <-  c("qbass", "gbass")
  return(md)
}







