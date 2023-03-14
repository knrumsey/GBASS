#' Build Prior
#'
#' This function provides a short cut for specifying GIG or GBP priors. See main help page of gbass() for more details.
#'

#' @export
#' @examples
#' #The following are equivalent
#' w_prior <- list(type="GBP", p=1, a=1/2, b=1/2, prop_sigma = 0.2)
#' w_prior <- build_prior("GBP", c(1, 1/2, 1/2), 0.2)

#' @name build_prior
#' @rdname build_prior
#'
#' @title Build GIG and GBP priors
#'
#' @param type either "GIG" or "GBP"
#' @param pars a vector (p, a, b) of parameters
#' @param prop_sigma the proposal standard deviation (only used for GBP or when symmetry parameter (beta) is non-zero)
#' @details Shortcut for specifying priors for the gbass() function. prop_sigma may not be necessary for GIG prior.
NULL

#' @rdname build_prior
#' @examples
#' # The following are equivalent
#' w_prior <- list(type="GBP", p=1, a=1/2, b=1/2, prop_sigma = 0.2)
#' w_prior <- build_prior("GBP", 1, 1/2, 1/2, 0.2)
#' w_prior <- build_GBP(1, 1/2, 1/2, 0.2)
#' @export
build_prior <- function(type, p, a, b, prop_sigma=NULL){
  pr <- list(type=type, p=p, a=a, b=b, prop_sigma=prop_sigma)
  return(pr)
}

#' @rdname build_prior
#' @export
build_GIG <- function(p, a, b, prop_sigma=NULL){
  if(is.null(prop_sigma)){
    pr <- list(type="GIG", p=p, a=a, b=b)
  }else{
    pr <- list(type="GIG", p=p, a=a, b=b, prop_sigma=prop_sigma)
  }

  return(pr)
}

#' @rdname build_prior
#' @export
build_GBP <- function(p, a, b, prop_sigma=NULL){
  pr <- list(type="GBP", p=p, a=a, b=b, prop_sigma=prop_sigma)
  return(pr)
}


