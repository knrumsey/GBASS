#' Construct a prior specification for GBASS
#'
#' @param type "GIG" or "GBP"
#' @param p,a,b prior hyperparameters
#' @param lower_bound optional lower bound
#' @param prop_sigma optional proposal sd (log scale)
#' @param adapt logical; use adaptive MH?
#' @param adapt_delay delay before adapting
#' @param adapt_thin thinning for adaptation updates
#'
#' @export
build_prior <- function(type = c("GIG", "GBP"),
                        p,
                        a,
                        b,
                        lower_bound = NULL,
                        prop_sigma = NULL,
                        adapt = NULL,
                        adapt_delay = NULL,
                        adapt_thin = NULL) {

  type <- match.arg(type)

  if (!is.numeric(p) || length(p) != 1L || is.na(p)) {
    stop("p must be a numeric scalar.")
  }
  if (!is.numeric(a) || length(a) != 1L || is.na(a)) {
    stop("a must be a numeric scalar.")
  }
  if (!is.numeric(b) || length(b) != 1L || is.na(b)) {
    stop("b must be a numeric scalar.")
  }

  out <- list(
    type = type,
    p = p,
    a = a,
    b = b
  )

  if (!is.null(lower_bound)) {
    if (!is.numeric(lower_bound) || length(lower_bound) != 1L || is.na(lower_bound)) {
      stop("lower_bound must be a numeric scalar.")
    }
    out$lower_bound <- lower_bound
  }

  if (!is.null(prop_sigma)) {
    if (!is.numeric(prop_sigma) || length(prop_sigma) != 1L ||
        is.na(prop_sigma) || prop_sigma <= 0) {
      stop("prop_sigma must be a positive numeric scalar.")
    }
    out$prop_sigma <- prop_sigma
  }

  if (!is.null(adapt)) {
    if (!is.logical(adapt) || length(adapt) != 1L || is.na(adapt)) {
      stop("adapt must be a logical scalar.")
    }
    out$adapt <- adapt
  }

  if (!is.null(adapt_delay)) {
    if (!is.numeric(adapt_delay) || length(adapt_delay) != 1L ||
        is.na(adapt_delay) || adapt_delay < 0 ||
        adapt_delay != as.integer(adapt_delay)) {
      stop("adapt_delay must be a nonnegative integer.")
    }
    out$adapt_delay <- as.integer(adapt_delay)
  }

  if (!is.null(adapt_thin)) {
    if (!is.numeric(adapt_thin) || length(adapt_thin) != 1L ||
        is.na(adapt_thin) || adapt_thin < 1 ||
        adapt_thin != as.integer(adapt_thin)) {
      stop("adapt_thin must be a positive integer.")
    }
    out$adapt_thin <- as.integer(adapt_thin)
  }

  out
}

#' @rdname build_prior
#' @export
build_GIG <- function(p, a, b,
                      lower_bound = NULL,
                      prop_sigma = NULL,
                      adapt = NULL,
                      adapt_delay = NULL,
                      adapt_thin = NULL) {

  build_prior(
    type = "GIG",
    p = p,
    a = a,
    b = b,
    lower_bound = lower_bound,
    prop_sigma = prop_sigma,
    adapt = adapt,
    adapt_delay = adapt_delay,
    adapt_thin = adapt_thin
  )
}

#' @rdname build_prior
#' @export
build_GBP <- function(p, a, b,
                      lower_bound = NULL,
                      prop_sigma = NULL,
                      adapt = NULL,
                      adapt_delay = NULL,
                      adapt_thin = NULL) {

  build_prior(
    type = "GBP",
    p = p,
    a = a,
    b = b,
    lower_bound = lower_bound,
    prop_sigma = prop_sigma,
    adapt = adapt,
    adapt_delay = adapt_delay,
    adapt_thin = adapt_thin
  )
}
