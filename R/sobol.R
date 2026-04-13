#' Sobol decomposition for GBASS models
#'
#' Computes Sobol sensitivity indices for a fitted \code{gbass} model by first
#' converting it to a compatible \code{BASS} object with \code{gbass2bass()},
#' and then calling \code{BASS::sobol()}.
#'
#' @param object A fitted object of class \code{"gbass"}.
#' @param ... Additional arguments passed to \code{BASS::sobol()}.
#'
#' @return An object of class \code{"bassSob"} returned by \code{BASS::sobol()}.
#'
#' @details
#' This is a thin wrapper around \code{BASS::sobol()} for GBASS models.
#' Users who want direct access to the converted \code{BASS} object can call
#' \code{gbass2bass()} explicitly and then use \code{BASS::sobol()} themselves.
#'
#' @seealso \code{\link{gbass2bass}}, \code{\link[BASS]{sobol}}
#'
#' @examples
#' ff1 <- function(x) 10.391*((x[1]-0.4)*(x[2]-0.6) + 0.36)
#'
#' n <- 100
#' p <- 4
#' X <- matrix(runif(n * p), nrow = n)
#' y <- apply(X, 1, ff1)
#'
#' mod <- gbass(X, y, nmcmc = 1000, nburn = 900)
#'
#' # Direct wrapper
#' sob <- gsobol(mod)
#'
#' # Equivalent manual conversion
#' bm <- gbass2bass(mod)
#' sob2 <- BASS::sobol(bm)
#'
#' @export
gsobol <- function(object, ...) {
  BASS::sobol(gbass2bass(object), ...)
}
gsobol <- function(object, ...) {
  BASS::sobol(gbass2bass(object), ...)
}
