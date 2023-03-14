#' @name funcs
#' @rdname funcs
#'
#' @title Bivariate functions for Denison, Mallick and Smith
#'
#' @param x a vector of inputs. length(x) must be at least 2, but only the first two entries are used.
#' @details The five test functions from Denison, Mallick and Smith (1998) are referred to as
#' \enumerate{
#'    \item Simple
#'    \item Radial
#'    \item Harmonic
#'    \item Additive
#'    \item Complex interaction
#' }
NULL

#' @rdname funcs
#' @examples
#' n <- 100 #Number of observations
#' p <- 4   #Number of variables (beyond p = 2, variables are inert)
#' X <- matrix(runif(n*p), nrow=n)
#' y <- apply(X, 1, ff1)
#' @export
ff1 <- function(x) 10.391*((x[1]-0.4)*(x[2]-0.6) + 0.36)

#' @rdname funcs
#' @export
ff2 <- function(x){ r <- (x[1]-0.5)^2 + (x[2] - 0.5)^2; 24.234*(r*(0.75 - r))}

#' @rdname funcs
#' @export
ff3 <- function(x){ xx1 <- x[1] - 0.5; xx2 <- x[2] - 0.5; 42.659*(0.1 + xx1*(0.05 + xx1^4 - 10*xx1^2*xx2^2 + 5*xx2^4))}

#' @rdname funcs
#' @export
ff4 <- function(x) 1.3356*(1.5*(1-x[1]) + exp(2*x[1] - 1)*sin(3*pi*(x[1] - 0.6)^2) + exp(3*(x[2]-0.5))*sin(4*pi*(x[2] - 0.9)^2))

#' @rdname funcs
#' @export
ff5 <- function(x) 1.9*(1.35 + exp(x[1])*sin(13*(x[1]-0.6)^2)*exp(-x[2])*sin(7*x[2]))





