#' Fay-Herriot BMARS (fhbass2)
#'
#' A version of BMARS which implements a Fay-Herriot style model:
#'   y_i = B(x_i) + sqrt(w * v_i) * eps_i,   eps_i ~ N(0,1)
#' with a GIG prior on v_i centered around supplied delta_i = d_i, and a GIG prior on w.
#' Includes optional Metropolis updates for gamma_v and gamma_w with half-Cauchy priors.
#'
#' IMPORTANT DESIGN CHOICES (vs your draft):
#'  - Removes the "beta * v * sqrt(w)" mean term entirely (beta fixed to 0).
#'  - Fixes scoping bugs in MH log-densities (closures inside fhbass2).
#'  - Fixes predict method dispatch + indexing, and adds posterior predictive option.
#'  - Separates "prior scale" from "initial value" for gamma parameters.
#'
#' @param X an Nxp matrix of predictors scaled to (0,1) recommended.
#' @param y length-N response vector.
#' @param d length-N vector of sampling variances (delta_i).
#' @param maxInt maximum interaction degree for BMARS bases.
#' @param maxBasis maximum number of basis functions.
#' @param npart minimum number of nonzero entries in a candidate basis.
#' @param nmcmc number of MCMC iterations.
#' @param nburn burn-in iterations (<= nmcmc).
#' @param thin thinning interval.
#' @param moveProbs probabilities for birth/death/mutate moves.
#' @param gamma_v_prior_scale half-Cauchy prior scale for gamma_v (>0).
#' @param gamma_w_prior_scale half-Cauchy prior scale for gamma_w (>0). If NULL, computed from fh_sample_var.
#' @param fh_sample_var prior guess of "between-area" variance (FH sampling error variance proxy).
#' @param gamma_v_init initial value for gamma_v (defaults to gamma_v_prior_scale).
#' @param gamma_w_init initial value for gamma_w (defaults to gamma_w_prior_scale).
#' @param do_mh_gamma logical; if FALSE, gamma_v and gamma_w are held fixed at initial values.
#' @param a_tau,b_tau IG prior parameters for tau (penalty/coeff scale parameter).
#' @param a_lambda,b_lambda Gamma prior for lambda (basis count penalty).
#' @param scale fixed variance scaling (default 1).
#' @param Iw0 nominal weights for interaction degree sampling.
#' @param Zw0 nominal weights for variable sampling.
#' @param verbose print progress.
#'
#' @return an object of class "fhbass2" containing posterior draws and basis lookup.
#' @export
fhbass2 <- function(X, y, d,
                    maxInt=3, maxBasis=1000, npart=NULL,
                    nmcmc=10000, nburn=9001, thin=1,
                    moveProbs=rep(1/3,3),
                    gamma_v_prior_scale=100,
                    gamma_w_prior_scale=NULL,
                    fh_sample_var=NULL,
                    gamma_v_init=NULL,
                    gamma_w_init=NULL,
                    do_mh_gamma=TRUE,
                    a_tau=1/2, b_tau=NULL,
                    a_lambda=1, b_lambda=1,
                    scale=1,
                    Iw0 = rep(1, maxInt), Zw0 = rep(1, ncol(X)),
                    verbose=TRUE){

  # -----------------------
  # Basic checks/standardize
  # -----------------------
  if(is.null(dim(X))) X <- matrix(X, ncol=1)
  if(nrow(X) != length(y)) stop("nrow(X) and length(y) should match")
  N <- nrow(X)
  p <- ncol(X)
  if(length(d) != N) stop("d must be a vector of length nrow(X)")
  if(any(d <= 0)) stop("All d must be positive")

  if(min(X) < 0 || max(X) > 1) warning("Found values of X outside (0,1). Basis generation assumes (0,1) scaling.")
  maxInt <- min(p, maxInt)
  Iw0 <- Iw0[seq_len(maxInt)]
  if(length(Zw0) != p) stop("Zw0 must have length ncol(X)")
  if(is.null(npart)) npart <- min(20, 0.1*N)
  if(is.null(b_tau)) b_tau <- N/2

  if(nburn > nmcmc) stop("nburn must be <= nmcmc")
  nkeep <- length(seq(nburn, nmcmc, by=thin))

  # -----------------------
  # Fay-Herriot prior calibration for gamma_w prior scale
  # -----------------------
  delta <- as.numeric(d)

  if(is.null(gamma_w_prior_scale)){
    if(is.null(fh_sample_var)){
      fh_sample_var <- stats::var(y)
    }
    # Crude heuristic: total var ~ mean(delta) + fh_sample_var (or similar)
    # We want a positive scale. Use fallback quantiles if needed.
    v_bar <- mean(delta) * (1 + 1/max(1e-8, gamma_v_prior_scale))
    gamma_w_prior_scale <- v_bar / (fh_sample_var - v_bar)
    if(gamma_w_prior_scale <= 0){
      v_bar <- mean(delta)
      gamma_w_prior_scale <- v_bar / (fh_sample_var - v_bar)
    }
    if(gamma_w_prior_scale <= 0){
      v_bar <- as.numeric(stats::quantile(delta, 0.25))
      gamma_w_prior_scale <- v_bar / (fh_sample_var - v_bar)
    }
    if(!is.finite(gamma_w_prior_scale) || gamma_w_prior_scale <= 0){
      stop("Cannot find a good prior scale for gamma_w. Please specify gamma_w_prior_scale directly.")
    }
  }

  # Initial values for gammas
  if(is.null(gamma_v_init)) gamma_v_init <- gamma_v_prior_scale
  if(is.null(gamma_w_init)) gamma_w_init <- gamma_w_prior_scale
  gamma_v <- as.numeric(gamma_v_init)
  gamma_w <- as.numeric(gamma_w_init)
  if(gamma_v <= 0 || gamma_w <= 0) stop("gamma_v_init and gamma_w_init must be > 0")

  # Priors (as lists) - updated each iteration when gammas change
  v_prior_fh <- list(type="GIG", p=1/2, a=gamma_v / delta, b=gamma_v * delta, lb=0)
  w_prior    <- list(type="GIG", p=1/2, a=gamma_w*(1+1/gamma_v), b=gamma_w/(1+1/gamma_v), lb=1/N^2)

  # -----------------------
  # MH tuning/accept tracking for gammas
  # -----------------------
  prop_sigma_gamma_v <- mean(delta) / 10
  prop_sigma_gamma_w <- stats::var(y) / 10
  accept_gamma_v <- 0
  accept_gamma_w <- 0

  # -----------------------
  # Allocate storage
  # -----------------------
  a_mc <- vector("list", nkeep)
  basis_mc <- vector("list", nkeep)

  M_mc <- lam_mc <- tau_mc <- w_mc <- gamma_v_mc <- gamma_w_mc <- ss <- rep(NA_real_, nkeep)
  v_mc <- matrix(NA_real_, nrow=nkeep, ncol=N)

  lookup <- list()
  basis_index <- integer(0)
  nlook <- 1L
  kk <- 1L

  # -----------------------
  # Initialize BMARS state
  # -----------------------
  M   <- 0L
  lam <- rgamma(1, a_lambda, b_lambda)
  tau <- (b_tau+.1)/(a_tau+.1)

  # Start v at delta; start w at 1 (or draw from prior if you prefer)
  v <- delta
  w <- 1.0

  B <- matrix(1, nrow=N, ncol=1)
  tX <- t(X)

  # Helper: compute U, U2 for current B, v, tau, w
  recompute_linear_cache <- function(B, v, tau, w, y, scale){
    Vinv <- Matrix::Diagonal(x = 1/as.numeric(v))
    U  <- solve(symchol(crossprod(B, Vinv) %*% B + scale/tau * Matrix::Diagonal(ncol(B))))
    # Here z = y (beta removed)
    U2 <- Matrix::crossprod(y/v, B) %*% U
    list(Vinv=Vinv, U=U, U2=U2)
  }

  cache <- recompute_linear_cache(B, v, tau, w, y, scale)
  U <- cache$U
  U2 <- cache$U2
  Vinv <- cache$Vinv

  cnt1 <- cnt2 <- rep(0, 3)

  if(verbose){
    cat(c('MCMC iteration',0,myTimestamp(),'nbasis:',M), '\n')
  }

  # -----------------------
  # Closures for MH log-densities (fix scoping)
  # -----------------------
  log_halfcauchy <- function(x, s){
    if(x <= 0) return(-Inf)
    # constant log(2) cancels in MH, but keep for clarity
    log(2) + dcauchy(x, location=0, scale=s, log=TRUE)
  }

  # log p(gamma_v | v, delta) up to additive constant
  logdens_gamma_v <- function(gv, v, delta, gamma_v_prior_scale){
    if(gv <= 0) return(-Inf)
    a_v <- gv / delta
    b_v <- gv * delta
    # v | gv density (product over i) * prior(gv)
    sum(dgig(v, p = 1/2, a = a_v, b = b_v, log = TRUE)) + log_halfcauchy(gv, gamma_v_prior_scale)
  }

  # log p(gamma_w | w, gamma_v) up to constant
  logdens_gamma_w <- function(gw, w, gamma_v, N_eff, gamma_w_prior_scale){
    if(gw <= 0) return(-Inf)
    a_w <- gw * (1 + 1/gamma_v)
    b_w <- gw / (1 + 1/gamma_v)
    # w | gw density * prior(gw)
    dgig(w, p = 1/2, a = a_w, b = b_w, log = TRUE) + log_halfcauchy(gw, gamma_w_prior_scale)
  }

  logdens_kappa_w <- function(kw, w, gamma_v, kappa_w_prior_scale){
    if(!is.finite(kw) || kw <= 0) return(-Inf)
    if(!is.finite(w)  || w  <= 0) return(-Inf)
    if(!is.finite(gamma_v) || gamma_v <= 0) return(-Inf)

    gamma_w <- 1/kw

    a_w <- gamma_w * (1 + 1/gamma_v)
    b_w <- gamma_w / (1 + 1/gamma_v)

    # Likelihood term: w ~ GIG(p=1/2, a=a_w, b=b_w)
    dgig(w, p = 1/2, a = a_w, b = b_w, log = TRUE) +
      log_halfcauchy(kw, kappa_w_prior_scale)
  }

  # -----------------------
  # Main MCMC
  # -----------------------
  for(k in 1:nmcmc){

    move <- move_type(M, maxBasis, moveProbs)

    # ======================================================
    # BIRTH
    # ======================================================
    if(move == "B"){
      J_cand <- sample(maxInt, 1, prob=Iw0)
      u_cand <- sample(p, J_cand, replace=FALSE, prob=Zw0)
      s_cand <- 1 - 2*rbinom(J_cand, 1, 0.5)
      t_cand <- runif(J_cand)

      B_new <- makeBasis(s_cand, u_cand, t_cand, tX, 1)

      if(sum(B_new > 0) >= npart){
        B_cand <- cbind(B, B_new)
        U_cand <- tryCatch(
          solve(symchol(crossprod(B_cand, Vinv) %*% B_cand + scale/tau*Matrix::Diagonal(ncol(B_cand)))),
          error=function(e) FALSE
        )

        if(isFALSE(U_cand)){
          log_accept_prob_B <- -Inf
        }else{
          U2_cand <- Matrix::crossprod(y/v, B_cand) %*% U_cand

          log_accept_prob_B <- sum(log(abs(diag(U_cand)))) - sum(log(abs(diag(U)))) +
            1/(2*w*scale) * (TCP(U2_cand) - TCP(U2))

          log_accept_prob_B <- log_accept_prob_B +
            log(scale)/2 - log(tau)/2 + log(lam) +
            log(moveProbs[2]) - log(moveProbs[1]) -
            log(M+1) - log(maxInt) - lchoose(p, J_cand) -
            log(Iw0[J_cand]/sum(Iw0)) - log(dmwnchBass(Zw0, u_cand))
        }
      }else{
        log_accept_prob_B <- -Inf
      }

      cnt1[1] <- cnt1[1] + 1
      if(log(runif(1)) < as.numeric(log_accept_prob_B)){
        cnt2[1] <- cnt2[1] + 1
        B <- B_cand
        U <- U_cand
        U2 <- U2_cand
        Iw0[J_cand] <- Iw0[J_cand] + 1
        Zw0[u_cand] <- Zw0[u_cand] + 1
        M <- M + 1L
        lookup[[nlook]] <- list(J=J_cand, s=s_cand, t=t_cand, u=u_cand)
        basis_index <- c(basis_index, nlook)
        nlook <- nlook + 1L
      }
    }

    # ======================================================
    # DEATH
    # ======================================================
    if(move == "D"){
      m0 <- sample(M, 1)
      basis_cand <- lookup[[basis_index[m0]]]
      J_cand <- basis_cand$J
      u_cand <- basis_cand$u

      Zw_cand <- Zw0
      Zw_cand[u_cand] <- Zw_cand[u_cand] - 1

      B_cand <- B[,-(m0+1), drop=FALSE]
      U_cand <- tryCatch(
        solve(symchol(crossprod(B_cand, Vinv) %*% B_cand + scale/tau*Matrix::Diagonal(ncol(B_cand)))),
        error=function(e) FALSE
      )

      if(isFALSE(U_cand)){
        log_accept_prob_D <- -Inf
      }else{
        U2_cand <- Matrix::crossprod(y/v, B_cand) %*% U_cand

        log_accept_prob_D <- sum(log(abs(diag(U_cand)))) - sum(log(abs(diag(U)))) +
          1/(2*w*scale) * (TCP(U2_cand) - TCP(U2))

        log_accept_prob_D <- log_accept_prob_D -
          log(scale)/2 + log(tau)/2 - log(lam) -
          log(moveProbs[2]) + log(moveProbs[1]) +
          log(M) + log(maxInt) + lchoose(p, J_cand) +
          log((Iw0[J_cand]-1)/(sum(Iw0)-1)) + log(dmwnchBass(Zw_cand, u_cand))
      }

      cnt1[2] <- cnt1[2] + 1
      if(log(runif(1)) < as.numeric(log_accept_prob_D)){
        cnt2[2] <- cnt2[2] + 1
        B <- B_cand
        U <- U_cand
        U2 <- U2_cand
        Iw0[J_cand] <- Iw0[J_cand] - 1
        Zw0 <- Zw_cand
        M <- M - 1L
        basis_index <- basis_index[-m0]
      }
    }

    # ======================================================
    # MUTATE
    # ======================================================
    if(move == "M"){
      m0 <- sample(M, 1)
      basis_cand <- lookup[[basis_index[m0]]]
      J_cand <- basis_cand$J
      u_cand <- basis_cand$u
      s_cand <- 1 - 2*rbinom(J_cand, 1, 0.5)
      t_cand <- runif(J_cand)

      B_new <- makeBasis(s_cand, u_cand, t_cand, tX, 1)

      if(sum(B_new > 0) >= npart){
        B_cand <- B
        B_cand[,(m0+1)] <- B_new

        U_cand <- tryCatch(
          solve(symchol(crossprod(B_cand, Vinv) %*% B_cand + scale/tau*Matrix::Diagonal(ncol(B_cand)))),
          error=function(e) FALSE
        )

        if(isFALSE(U_cand)){
          log_accept_prob_M <- -Inf
        }else{
          U2_cand <- Matrix::crossprod(y/v, B_cand) %*% U_cand
          log_accept_prob_M <- sum(log(abs(diag(U_cand)))) - sum(log(abs(diag(U)))) +
            1/(2*w*scale) * (TCP(U2_cand) - TCP(U2))
        }
      }else{
        log_accept_prob_M <- -Inf
      }

      cnt1[3] <- cnt1[3] + 1
      if(log(runif(1)) < as.numeric(log_accept_prob_M)){
        cnt2[3] <- cnt2[3] + 1
        B <- B_cand
        U <- U_cand
        U2 <- U2_cand
        lookup[[nlook]] <- list(J=J_cand, s=s_cand, t=t_cand, u=u_cand)
        basis_index[m0] <- nlook
        nlook <- nlook + 1L
      }
    }

    # ======================================================
    # GIBBS STEPS (beta removed; z = y)
    # ======================================================
    # Coefficients
    a <- TCP(U, U2) + sqrt(scale*w) * (U %*% rnorm(M+1))
    yhat <- as.numeric(B %*% a)
    r <- y - yhat
    pen <- sum(a^2)

    # Lambda, tau
    lam <- rgamma(1, a_lambda + M, b_lambda + 1)
    tau <- 1/rgamma(1, a_tau + (M+1)/2, b_tau + pen/(2*w))

    # Update w (GIG full conditional)
    # w | ... ~ GIG(p = 1/2 - (N+M+1)/2, a = w_prior$a, b = w_prior$b + sum(r^2/v)/scale + pen/tau)
    w_cand <- rgig2(
      p = 1/2 - (N + M + 1)/2,
      a = w_prior$a,
      b = w_prior$b + sum((r^2)/v)/scale + pen/tau
    )
    if(is.finite(w_cand) && w_cand > w_prior$lb) w <- w_cand

    # Update v_i (vectorized FH GIG update)
    # v_i | ... ~ GIG(p = 1/2 - 1/2 = 0, a = gamma_v/delta, b = gamma_v*delta + r_i^2/(w*scale))
    v <- rgig2.vec_fh(
      p = 1/2 - 1/2,
      a = (gamma_v / delta),
      b = as.numeric(gamma_v * delta + (r^2)/(w*scale))
    )

    # ======================================================
    # Update gammas via MH (optional)
    # ======================================================
    if(do_mh_gamma){
      # gamma_v
      res_v <- mh_update_gamma_scalar(
        gamma = gamma_v,
        logdens = function(gv) logdens_gamma_v(gv, v=v, delta=delta, gamma_v_prior_scale=gamma_v_prior_scale),
        prop_sigma = prop_sigma_gamma_v
      )
      gamma_v <- res_v$gamma
      accept_gamma_v <- accept_gamma_v + as.integer(res_v$accept)

      # gamma_w (depends on gamma_v)
      if(k == 791){
        #browser()
      }
      res_w <- mh_update_gamma_scalar(
        gamma = gamma_w,
        logdens = function(gw) logdens_gamma_w(gw, w=w, gamma_v=gamma_v, N_eff=(N+M+1),
                                               gamma_w_prior_scale=gamma_w_prior_scale),
        prop_sigma = prop_sigma_gamma_w
      )
      gamma_w <- min(max(0.01, res_w$gamma), 100)

      accept_gamma_w <- accept_gamma_w + as.integer(res_w$accept)

      # adapt every 50
      if (k %% 50 == 0) {
        prop_sigma_gamma_v <- prop_sigma_gamma_v * exp(0.3 * ((accept_gamma_v / 50) - 0.44))
        prop_sigma_gamma_w <- prop_sigma_gamma_w * exp(0.3 * ((accept_gamma_w / 50) - 0.44))
        accept_gamma_v <- 0
        accept_gamma_w <- 0
      }
    }

    # Refresh priors with updated gammas
    v_prior_fh <- list(type="GIG", p=1/2, a=gamma_v / delta, b=gamma_v * delta, lb=0)
    w_prior    <- list(type="GIG", p=1/2, a=gamma_w*(1+1/gamma_v), b=gamma_w/(1+1/gamma_v), lb=1/N^2)

    # Recompute caches
    Vinv <- Matrix::Diagonal(x=1/as.numeric(v))
    U <- solve(symchol(crossprod(B, Vinv) %*% B + scale/tau*Matrix::Diagonal(ncol(B))))
    U2 <- Matrix::crossprod(y/v, B) %*% U

    # Store draws
    if(k >= nburn && ((k-nburn) %% thin) == 0){
      ss[kk] <- mean((y-yhat)^2)
      M_mc[kk] <- M
      lam_mc[kk] <- lam
      tau_mc[kk] <- tau
      w_mc[kk] <- w
      v_mc[kk,] <- v
      gamma_v_mc[kk] <- gamma_v
      gamma_w_mc[kk] <- gamma_w

      basis_mc[[kk]] <- basis_index
      a_mc[[kk]] <- a
      kk <- kk + 1L
    }

    if(verbose && k %% 1000 == 0){
      cat(c('MCMC iteration',k,myTimestamp(),'nbasis:',M), '\n')
    }
  }

  obj <- list(
    nbasis=M_mc, w=w_mc, v=v_mc,
    tau=tau_mc, lam=lam_mc,
    a=a_mc, basis=basis_mc, lookup=lookup,
    cnt1=cnt1, cnt2=cnt2, ss=ss,
    delta=delta,
    gamma_v=gamma_v_mc, gamma_w=gamma_w_mc,
    gamma_v_prior_scale=gamma_v_prior_scale,
    gamma_w_prior_scale=gamma_w_prior_scale,
    X=X, y=y,
    call=match.call()
  )
  class(obj) <- "fhbass2"
  return(obj)
}


# ============================================================
# Predict method for fhbass2
# ============================================================
#' @export
predict.fhbass2 <- function(object, newdata=NULL, d_new=NULL,
                            mcmc.use=NULL,
                            type=c("theta","yrep"),
                            fallback=c("error","mean"),
                            seed=NULL){

  type <- match.arg(type)
  fallback <- match.arg(fallback)

  if(is.null(newdata)) newdata <- object$X
  if(is.null(dim(newdata))) newdata <- matrix(newdata, ncol=ncol(object$X))
  Nnew <- nrow(newdata)

  if(is.null(mcmc.use)) mcmc.use <- seq_along(object$nbasis)

  # Resolve d_new if needed
  if(type == "yrep"){
    if(is.null(d_new)){
      if(fallback == "error"){
        stop("type='yrep' requires d_new, or set fallback='mean' to use mean training delta.")
      } else {
        d_new <- rep(mean(object$delta), Nnew)
      }
    }
    if(length(d_new) != Nnew) stop("d_new must have length nrow(newdata)")
    if(any(d_new <= 0)) stop("d_new must be positive")
  }

  if(!is.null(seed)) set.seed(seed)

  tX <- t(newdata)
  nd <- length(mcmc.use)
  out <- matrix(NA_real_, nrow=nd, ncol=Nnew)

  jj <- 1L
  for(i in mcmc.use){
    basis_curr <- object$lookup[unlist(object$basis[[i]])]
    B_curr <- matrix(1, nrow=Nnew, ncol=length(basis_curr)+1)
    if(length(basis_curr) > 0){
      for(m in seq_along(basis_curr)){
        thetaB <- basis_curr[[m]]
        B_curr[,m+1] <- makeBasis(thetaB$s, thetaB$u, thetaB$t, tX, 1)
      }
    }
    theta_draw <- as.numeric(B_curr %*% object$a[[i]])

    if(type == "theta"){
      out[jj,] <- theta_draw
    } else {
      # posterior predictive draw:
      # yrep = theta + sqrt(w * v_pred) * eps
      # For new points we only know d_new; we use v_pred ~ GIG(p=1/2, a=gamma_v/d_new, b=gamma_v*d_new)
      gv <- object$gamma_v[i]
      gw <- object$gamma_w[i]
      w  <- object$w[i]

      v_pred <- rgig2.vec_fh(
        p = 1/2,               # prior predictive for v at new point
        a = gv / d_new,
        b = gv * d_new
      )
      out[jj,] <- theta_draw + sqrt(w * v_pred) * rnorm(Nnew)
    }
    jj <- jj + 1L
  }

  out
}


# ============================================================
# Helper: MH update for a positive scalar on log-scale
# ============================================================
mh_update_gamma_scalar <- function(gamma, logdens, prop_sigma){
  log_gamma_prop <- log(gamma) + rnorm(1, 0, prop_sigma)
  gamma_prop <- exp(log_gamma_prop)

  log_alpha <- logdens(gamma_prop) - logdens(gamma)
  if(log(runif(1)) < log_alpha){
    list(gamma=gamma_prop, accept=TRUE)
  } else {
    list(gamma=gamma, accept=FALSE)
  }
}


# ============================================================
# GIG density (used only for MH; use stable implementation if available)
# ============================================================
dgig <- function(x, p, a, b, log = FALSE) {
  if (any(x <= 0)) stop("x must be positive")
  if (any(a <= 0) || any(b <= 0)) stop("a and b must be positive")

  # vectorized a,b to match x
  if(length(a) == 1) a <- rep(a, length(x))
  if(length(b) == 1) b <- rep(b, length(x))

  sqrt_ab <- sqrt(a * b)
  log_kernel <- (p - 1) * log(x) - 0.5 * (a * x + b / x)

  # Normalizing constant:
  # f(x) = (a/b)^(p/2) / (2 K_p(sqrt(ab))) x^(p-1) exp(-0.5(ax + b/x))
  # logC = (p/2)log(a/b) - log(2) - log(K_p(sqrt(ab)))
  logC <- (p/2) * log(a / b) - log(2) - log(besselK(sqrt_ab, p))

  out <- logC + log_kernel
  if(log) out else exp(out)
}
