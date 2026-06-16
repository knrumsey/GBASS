#' Convert GBASS object to BASS-like object
#'
#' Converts an object of class \code{gbass} to an object with class \code{bass},
#' so that downstream BASS utilities such as Sobol decomposition can be used.
#'
#' @param gm An object of class \code{gbass}.
#' @param knot_strategy Character string specifying how GBASS knot locations are
#'   represented in the converted BASS object. If \code{"exact"}, artificial rows
#'   are appended to \code{xx.des} so that the original GBASS knots are represented
#'   exactly. If \code{"nearest"}, each knot is approximated by the nearest
#'   observed marginal input value in \code{gm$X}.
#'
#' @return An object with class \code{c("bass", "gbass")}. The returned object is
#'   intended mainly for downstream BASS utilities such as \code{BASS::sobol()}.
#' @export
gbass2bass <- function(gm, knot_strategy = c("exact", "nearest")) {
  knot_strategy <- match.arg(knot_strategy)

  if (!inherits(gm, "gbass")) {
    stop("gm must be an object of class 'gbass'")
  }

  needed <- c("a", "M", "X", "basis", "lookup")
  missing_fields <- needed[!needed %in% names(gm)]
  if (length(missing_fields) > 0) {
    stop("gbass object is missing required fields: ",
         paste(missing_fields, collapse = ", "))
  }

  out <- list()
  out$degree <- 1

  nmcmc <- length(gm$a)
  maxb <- max(gm$M)

  beta <- matrix(NA_real_, nrow = nmcmc, ncol = maxb + 1)
  for (i in seq_len(nmcmc)) {
    beta[i, seq_along(gm$a[[i]])] <- gm$a[[i]]
  }

  out$beta <- beta
  out$nbasis <- gm$M
  out$p <- ncol(gm$X)

  model.lookup <- integer(nmcmc)
  model.lookup[1] <- 1L

  if (nmcmc > 1L) {
    for (i in 2:nmcmc) {
      same_model <- gm$M[i] == gm$M[i - 1] &&
        identical(gm$basis[[i]], gm$basis[[i - 1]])

      model.lookup[i] <- if (same_model) model.lookup[i - 1] else model.lookup[i - 1] + 1L
    }
  }

  out$model.lookup <- model.lookup
  out$n.models <- max(model.lookup)

  out$type <- "_des"
  out$des <- TRUE
  out$func <- FALSE
  out$cat <- FALSE

  max.int <- max(vapply(gm$lookup, function(zz) length(zz$s), integer(1)))

  n.int.des <- matrix(0, nrow = out$n.models, ncol = maxb)
  signs.des <- vars.des <- knotInd.des <- array(
    NA_integer_,
    dim = c(out$n.models, maxb, max.int)
  )
  knots.des <- array(NA_real_, dim = c(out$n.models, maxb, max.int))

  xx.des <- gm$X
  template_row <- colMeans(gm$X)
  knot_error <- numeric(0)

  add_knot_row <- function(var, knot) {
    new_row <- template_row
    new_row[var] <- knot
    xx.des <<- rbind(xx.des, new_row)
    nrow(xx.des)
  }

  nearest_knot_index <- function(var, knot) {
    idx <- which.min(abs(xx.des[, var] - knot))
    knot_error <<- c(knot_error, abs(xx.des[idx, var] - knot))
    idx
  }

  get_knot_index <- function(var, knot) {
    if (knot_strategy == "exact") {
      add_knot_row(var = var, knot = knot)
    } else {
      nearest_knot_index(var = var, knot = knot)
    }
  }

  for (i in seq_len(out$n.models)) {
    ind <- which(model.lookup == i)
    draw_idx <- ind[1]

    if (gm$M[draw_idx] > 0) {
      basis_ids <- unlist(gm$basis[[draw_idx]])

      for (j in seq_len(gm$M[draw_idx])) {
        basis_id <- basis_ids[j]
        basis_obj <- gm$lookup[[basis_id]]

        J <- basis_obj$J
        ind_j <- seq_len(J)

        n.int.des[i, j] <- J
        signs.des[i, j, ind_j] <- basis_obj$s
        vars.des[i, j, ind_j] <- basis_obj$u
        knots.des[i, j, ind_j] <- basis_obj$t

        for (ell in ind_j) {
          knotInd.des[i, j, ell] <- get_knot_index(
            var = basis_obj$u[ell],
            knot = basis_obj$t[ell]
          )
        }
      }
    }
  }

  out$xx.des <- xx.des
  out$n.int.des <- n.int.des
  out$signs.des <- signs.des
  out$vars.des <- vars.des
  out$knotInd.des <- knotInd.des
  out$knots.des <- knots.des
  out$cx <- rep("numeric", out$p)
  out$range.des <- rbind(rep(0, out$p), rep(1, out$p))

  if (knot_strategy == "nearest") {
    out$knot_error <- knot_error
    out$max_knot_error <- if (length(knot_error) > 0) max(knot_error) else 0
    out$mean_knot_error <- if (length(knot_error) > 0) mean(knot_error) else 0
  }

  out$nburn <- 0
  out$thin <- 1
  out$nmcmc <- nmcmc
  out$pdes <- out$p
  out$pcat <- 0
  out$pfunc <- 0
  out$maxInt.des <- max.int
  out$maxInt.cat <- 0
  out$maxInt.func <- 0

  class(out) <- c("bass", class(gm))
  out
}
