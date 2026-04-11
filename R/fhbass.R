#' Fay-Herriot BMARS
#'
#' A version of BMARS which attempts to replicate the Fay-Herriot model.
#'
#' @param X an Nxp matrix of predictor variables.
#' @param y an Nx1 matrix (or vector) of response values.
#' @param d A vector of delta values (measured variance at each location). See details.
#' @param maxInt integer for maximum degree of interaction in spline basis functions. Defaults to the number of predictors, which could result in overfitting.
#' @param maxBasis maximum number of basis functions. This should probably only be altered if you run out of memory.
#' @param npart of non-zero points in a basis function. If the response is functional, this refers only to the portion of the basis function coming from the non-functional predictors. Defaults to 20 or 0.1 times the number of observations, whichever is smaller.
#' @param nmcmc number of mcmc iterations
#' @param nburn number of burn-in samples
#' @param thin thinning for mcmc
#' @param moveProbs a vector defining the probabilities for (i) birth (ii) death and (iii) mutation. Default is rep(1/3,3).
#' @param gamma_v_scale scale of the half-Cauchy prior for \code{gamma_v}.
#' @param gamma_w_scale scale of the half-Cauchy prior for \code{gamma_w}.
#' @param fh_sample_var a prior guess for the classical Fay-Herriot "sampling error variance". Used to generate prior for \code{gamma_w}.
#' @param a_tau prior for tau
#' @param b_tau prior for tyau
#' @param a_lambda prior for lambda
#' @param b_lambda prior for lambda
#' @param m_beta prior for beta
#' @param s_beta prior for beta
#' @param scale fixed variance parameter. default is one.
#' @param Iw0 vector of nominal weights for degree of interaction, used in generating candidate basis functions. Should have length equal to Jmax and have positive entries.
#' @param Zw0 vector of nominal weights for variable selection, used in generating candidate basis functions. Should have length equal to ncol(X) and have positive entries.
#' @param verbose Logical. Should gbass print completion status? Default TRUE
#' @details The prior \code{v_prior_fh} should be a named list with named fields \code{df} (a scalar) and \code{m, vhat} (n-vectors).
#' @import Matrix
#' @export
#' @examples
#' fm <- function (x){
#'   x[1]^2 + x[1] * x[2] +  x[2]^3 / 9
#' }
#'
#' fv <- function (x){
#'   1.3356 * (1.5 * (1 - x[1]) + exp(2 * x[1] - 1) * sin(3 * pi *
#'     (x[1] - 0.6)^2) + exp(3 * (x[2] - 0.5)) * sin(4 * pi * (x[2] -
#'      0.9)^2))
#' }
#'
#' # GENERATE DATA
#' n <- 50
#' sigma_v <- 0.1
#' X <- lhs::randomLHS(n, 2)
#' mu <- apply(X, 1, fm)
#' mu <- mu / sd(mu)
#' N_pop <- 1 + round(apply(X, 1, fv))
#' d <- 1/N_pop
#' y <- rnorm(n, mu, sqrt(d)) + rnorm(n, 0, sigma_v)
#'
#' # Fit model
#' fit <- fhbass(X, y, d)
#'
#' @export
fhbass <- function(X, y, d,
                  maxInt=3, maxBasis=1000, npart=NULL, nmcmc=10000, nburn=9001, thin=1,
                  moveProbs=rep(1/3,3),
                  gamma_v_scale=100, gamma_w_scale=NULL, fh_sample_var=NULL,
                  a_tau=1/2, b_tau=NULL,
                  a_lambda=1, b_lambda=1,
                  m_beta=0, s_beta=0,
                  scale=1,
                  Iw0 = rep(1, maxInt), Zw0 = rep(1, ncol(X)),
                  verbose=TRUE){
  v_prior_fh <- d
  if(is.null(dim(X))) X <- matrix(X, ncol=1)
  if(max(X) > 1 | min(X) < 0) warning("Found values of X outside of (0, 1).")
  if(nrow(X) != length(y)) stop("nrow(X) and length(y) should match")
  N <- nrow(X)
  p <- ncol(X)
  maxInt <- min(ncol(X), maxInt)
  Iw0 <- Iw0[seq_len(maxInt)]
  nkeep <- length(seq(nburn, nmcmc, by=thin))
  if(is.null(npart)) npart <- min(20, 0.1*N)
  if(is.null(b_tau)) b_tau <- N/2
  if((m_beta != 0 | s_beta != 0)){
    w_prior_prop_sigma <- var(y)/scale/2
  }
  if(min(X) < 0 | max(X) > 1) warning("Data out of range (0,1). Are you sure you want to do this?")

  # Fay-Herriot stuff
  if(length(v_prior_fh) != N) stop("v_prior_fh must be a vector of length = nrow(X)")
  if(!is.null(gamma_w_scale) && !is.null(fh_sample_var)) stop("Only one of gamma_w_scale or fh_sample_var can be specified.")
  gamma_v <- gamma_v_scale
  if(is.null(gamma_w_scale)){
    if(is.null(fh_sample_var)){
      fh_sample_var <- var(y)
    }
    v_bar <- mean(v_prior_fh)*(1 + 1/gamma_v)
    gamma_w_scale <- v_bar / (fh_sample_var - v_bar)
    if(gamma_w_scale <= 0){
      v_bar <- mean(v_prior_fh)
      gamma_w_scale <- v_bar / (fh_sample_var - v_bar)
    }
    if(gamma_w_scale <= 0){
      v_bar <- quantile(v_prior_fh, 0.25)
      gamma_w_scale <- v_bar / (fh_sample_var - v_bar)
    }
    if(gamma_w_scale <= 0){
      stop("Cannot find a good prior for gamma_w. Please specify gamma_w_scale directly.")
    }
  }
  gamma_w <- gamma_w_scale

  delta <- v_prior_fh
  v_prior_fh <- list(type="GIG", p=1/2, a=gamma_v / delta, b=gamma_v * delta, lb=0)
  w_prior <- list(type="GIG", p=1/2, a=gamma_w*(1+1/gamma_v), b=gamma_w/(1+1/gamma_v), lb=1/N^2)

  # Initial values
  prop_sigma_gamma_v <- mean(delta) / 10
  prop_sigma_gamma_w <- var(y) / 10
  cnt_gamma_v <- 0
  cnt_gamma_w <- 0
  accept_gamma_v <- 0
  accept_gamma_w <- 0

  #ALLOCATE STORAGE SPACE
  a_mc <- list()
  M_mc <- lam_mc <- tau_mc <- w_mc <- bet_mc <- ss <- bias_mc <- s2_mc <- gamma_v_mc <- gamma_w_mc <- rep(NA, nkeep)
  v_mc   <- matrix(NA, nrow=nkeep, ncol=N)
  lookup <- basis_mc <- list()
  basis_index <- integer(0)
  nlook <- 1
  kk    <- 1

  #INITIALIZE PARAMETERS
  M   <- 0
  tau <- (b_tau+.1)/(a_tau+.1)
  v   <- delta
  w   <- gamma_w_scale
  lam <- rgamma(1, a_lambda, b_lambda)
  bet <- rnorm(1, m_beta, s_beta)
  a   <- mean(y)
  z    <- y - bet*v*sqrt(w)
  B    <- matrix(1, nrow=N)
  Vinv <- Matrix::Diagonal(x=1/v)
  # 3/14/23 using crossprod for slight speedup
  #U    <- solve(symchol(t(B)%*%Vinv%*%B + scale/tau*Diagonal(M+1)))
  #U2   <- t(z/v)%*%B%*%U
  U    <- solve(symchol(crossprod(B, Vinv)%*%B + scale/tau*Diagonal(M+1)))
  U2   <- crossprod(z/v, B)%*%U
  cnt1 <- cnt2 <- rep(0, 3)
  if(w_prior$type == "GBP" || abs(bet) > 1e-9){
    cntw <- 0
  }
  tX <- t(X)
  if(verbose){
    pr<-c('MCMC iteration',0,myTimestamp(),'nbasis:',M)
    cat(pr,'\n')
  }


  for(k in 1:nmcmc){
    move <- move_type(M, maxBasis, moveProbs)
    # ======================================================
    #        BIRTH STEP
    # ======================================================
    if(move == "B"){
      #Propose a new basis function
      J_cand <- sample(maxInt, 1, prob=Iw0)
      u_cand <- sample(p, J_cand, replace=FALSE, prob=Zw0)
      s_cand <- 1 - 2*rbinom(J_cand, 1, 0.5)
      t_cand <- runif(J_cand)
      #B_new <- rep(1, N)
      #for(j in 1:J_cand) B_new <- B_new * pmax(0, s_cand[j]*(X[,u_cand[j]]-t_cand[j]))
      B_new <- makeBasis(s_cand,u_cand,t_cand,tX,1)
      if(sum(B_new > 0) >= npart){
        B_cand  <- cbind(B, B_new)
        #U_cand  <- tryCatch(solve(symchol(t(B_cand)%*%Vinv%*%B_cand + scale/tau*Diagonal(M+2))), error=function(e) FALSE)
        U_cand  <- tryCatch(solve(symchol(crossprod(B_cand, Vinv)%*%B_cand + scale/tau*Diagonal(M+2))),
                            error=function(e) FALSE)
        if(isFALSE(U_cand)){
          log_accept_prob_B <- -Inf
        }else{
          #U2_cand <- t(z/v)%*%B_cand%*%U_cand
          U2_cand <- crossprod(z/v, B_cand)%*%U_cand
          log_accept_prob_B <- sum(log(abs(diag(U_cand))))-sum(log(abs(diag(U)))) +
            1/(2*w*scale)*(TCP(U2_cand)-TCP(U2))
          #Add line here if Sigma != Identity
          log_accept_prob_B <- log_accept_prob_B+log(scale)/2-log(tau)/2+log(lam)+log(moveProbs[2])-log(moveProbs[1])-log(M+1)-log(maxInt)-lchoose(p,J_cand)
          log_accept_prob_B <- log_accept_prob_B-log(Iw0[J_cand]/sum(Iw0))-log(dmwnchBass(Zw0, u_cand))
        }
      }else{
        log_accept_prob_B <- -Inf
      }

      if(is.nan(as.numeric(log_accept_prob_B))) browser()
      cnt1[1] <- cnt1[1] + 1
      if(log(runif(1)) < as.numeric(log_accept_prob_B)){
        cnt2[1] <- cnt2[1] + 1
        B <- B_cand
        U <- U_cand
        U2 <- U2_cand
        Iw0[J_cand] <- Iw0[J_cand] + 1
        Zw0[u_cand] <- Zw0[u_cand] + 1
        M <- M + 1
        lookup[[nlook]] <- list(J=J_cand, s=s_cand, t=t_cand, u=u_cand)
        basis_index <- c(basis_index, nlook)
        nlook <- nlook + 1
      }
    }
    # ======================================================
    #        DEATH STEP
    # ======================================================
    if(move == "D"){
      m0 <- sample(M, 1)
      basis_cand <- lookup[[basis_index[m0]]]
      J_cand <- basis_cand$J
      u_cand <- basis_cand$u
      Zw_cand <- Zw0
      Zw_cand[u_cand] <- Zw_cand[u_cand] - 1
      B_cand <- B[,-(m0+1)]
      #U_cand  <- tryCatch(solve(symchol(t(B_cand)%*%Vinv%*%B_cand + scale/tau*Diagonal(M))), error=function(e) FALSE)
      U_cand  <- tryCatch(solve(symchol(crossprod(B_cand, Vinv)%*%B_cand + scale/tau*Diagonal(M))),
                          error=function(e) FALSE)
      if(isFALSE(U_cand)){
        log_accept_prob_D <- -Inf
      }else{
        #U2_cand <- t(z/v)%*%B_cand%*%U_cand
        U2_cand <- crossprod(z/v, B_cand)%*%U_cand
        log_accept_prob_D<- sum(log(abs(diag(U_cand))))-sum(log(abs(diag(U)))) +
          1/(2*w*scale)*(TCP(U2_cand)-TCP(U2))
        #Add line here if Sigma != Identity
        log_accept_prob_D <- log_accept_prob_D-log(scale)/2+log(tau)/2-log(lam)-log(moveProbs[2])+log(moveProbs[1])+log(M)+log(maxInt)+lchoose(p,J_cand)
        log_accept_prob_D <- log_accept_prob_D+log((Iw0[J_cand]-1)/(sum(Iw0)-1))+log(dmwnchBass(Zw_cand, u_cand))
      }

      cnt1[2] <- cnt1[2] + 1
      if(is.nan(as.numeric(log_accept_prob_D))) browser()
      if(log(runif(1)) < as.numeric(log_accept_prob_D)){
        cnt2[2] <- cnt2[2] + 1
        B <- B_cand
        U <- U_cand
        U2 <- U2_cand
        Iw0[J_cand] <- Iw0[J_cand] - 1
        Zw0 <- Zw_cand
        M <- M - 1
        basis_index <- basis_index[-m0]
      }
    }
    # ======================================================
    #        MUTATE STEP
    # ======================================================
    if(move == "M"){
      m0 <- sample(M, 1)
      basis_cand <- lookup[[basis_index[m0]]]
      J_cand <- basis_cand$J
      u_cand <- basis_cand$u
      s_cand <- 1 - 2*rbinom(J_cand, 1, 0.5)
      t_cand <- runif(J_cand)
      #B_new <- rep(1, N)
      #for(j in 1:J_cand) B_new <- B_new*pmax(0, s_cand[j]*(X[,u_cand[j]] - t_cand[j]))
      B_new <- makeBasis(s_cand,u_cand,t_cand,tX,1)
      if(sum(B_new > 0) > npart){
        B_cand <- B
        B_cand[,(m0+1)] <- B_new
        #U_cand  <- tryCatch(solve(symchol(t(B_cand)%*%Vinv%*%B_cand + scale/tau*Diagonal(M+1))), error=function(e) FALSE)
        U_cand  <- tryCatch(solve(symchol(crossprod(B_cand, Vinv)%*%B_cand + scale/tau*Diagonal(M+1))),
                            error=function(e) FALSE)
        if(isFALSE(U_cand)){
          log_accept_prob_M <- -Inf
        }else{
          #U2_cand <- t(z/v)%*%B_cand%*%U_cand
          U2_cand <- crossprod(z/v, B_cand)%*%U_cand
          log_accept_prob_M<- sum(log(abs(diag(U_cand))))-sum(log(abs(diag(U)))) +
            1/(2*w*scale)*(sum(U2_cand^2)-sum(U2^2))
          #Add line here if Sigma != Identity
        }
      }else{
        log_accept_prob_M <- -Inf
      }
      if(is.nan(as.numeric(log_accept_prob_M))) browser()
      cnt1[3] <- cnt1[3] + 1
      if(log(runif(1)) < as.numeric(log_accept_prob_M)){
        cnt2[3] <- cnt2[3] + 1
        B <- B_cand
        U <- U_cand
        U2 <- U2_cand
        lookup[[nlook]] <- list(J=J_cand, s=s_cand, t=t_cand, u=u_cand)
        basis_index[m0] <- nlook
        nlook <- nlook + 1
      }
    }

    # ======================================================
    #        GIBBS STEPS
    # ======================================================
    a <- TCP(U, U2) + sqrt(scale*w)*(U%*%rnorm(M+1))
    yhat <- B%*%a
    r <- y - yhat
    pen <- sum(a^2)
    lam <- rgamma(1, a_lambda[1]+M, b_lambda+1)
    tau <- 1/rgamma(1, a_tau+(M+1)/2, b_tau+pen/(2*w))
    bet <- rnorm(1, (s_beta^2*sum(r)/sqrt(w) + scale*m_beta)/(s_beta^2*sum(v) + scale),
                 sqrt(scale*s_beta^2/(s_beta^2*sum(v)+scale)))
    if(w_prior$type == "GIG"){
      if(abs(bet) < 1e-9){
        w_cand <- rgig2(p=w_prior$p - (N+M+1)/2,
                        a=w_prior$a,
                        b=w_prior$b + sum(r^2/v)/scale + pen/tau)
        w <- ifelse(w_cand < w_prior$lb, w, w_cand)
      }else{
        w_cand <- exp(log(w) + rnorm(1, 0, w_prior_prop_sigma))
        log_alpha_w <- (w_prior$p - (N+M+1)/2)*(log(w_cand)-log(w)) -
          0.5*(w_prior$a*(w_cand-w) +
                 (w_prior$b + sum(r^2/v)/scale + pen/tau)*(1/w_cand - 1/w) -
                 2*bet/scale*sum(r)*(1/sqrt(w_cand)-1/sqrt(w))) +
          log(w_cand > w_prior$lb)
        if(log_alpha_w > log(runif(1))){
          cntw <- cntw + 1
          w <- w_cand
        }
      }
    }else{
      stop("Must use GIG prior for Fay-Herriot BASS")
    }
    if(v_prior_fh$type == "GIG"){
      v <- rgig2.vec_fh(p=v_prior_fh$p-1/2,
                     a=v_prior_fh$a+bet^2/scale,
                     b=as.numeric(v_prior_fh$b + r^2/(w*scale)))
    }else{
      stop("Must use GIG prior for Fay-Herriot BASS")
    }

    # Update gamma's
    res_v <- mh_update_gamma(gamma_v, logdens_gamma_v, prop_sigma_gamma_v, gamma_v_scale, cnt_gamma_v)
    gamma_v <- res_v$gamma
    cnt_gamma_v <- res_v$cnt
    accept_gamma_v <- res_v$accept

    res_w <- mh_update_gamma(gamma_w, logdens_gamma_w, prop_sigma_gamma_w, gamma_w_scale, cnt_gamma_w)
    gamma_w <- res_w$gamma
    cnt_gamma_w <- res_w$cnt
    accept_gamma_w <- res_w$accept

    if (k %% 50 == 0) {
      accept_rate_gamma_v <- accept_gamma_v / 50
      prop_sigma_gamma_v <- prop_sigma_gamma_v * exp(0.3 * (accept_rate_gamma_v - 0.44))
      accept_gamma_v <- 0

      accept_rate_gamma_w <- accept_gamma_w / 50
      prop_sigma_gamma_w <- prop_sigma_gamma_w * exp(0.3 * (accept_rate_gamma_w - 0.44))
      accept_gamma_w <- 0
    }

    v_prior_fh <- list(type="GIG", p=1/2, a=gamma_v / delta, b=gamma_v * delta, lb=0)
    w_prior <- list(type="GIG", p=1/2, a=gamma_w*(1+1/gamma_v), b=gamma_w/(1+1/gamma_v), lb=1/N^2)


    z    <- y - bet*v*sqrt(w)
    Vinv <- Matrix::Diagonal(x=1/v)
    #U    <- solve(symchol(t(B)%*%Vinv%*%B + scale/tau*Diagonal(M+1)))
    #U2   <- t(z/v)%*%B%*%U
    U     <- solve(symchol(crossprod(B, Vinv)%*%B + scale/tau*Diagonal(M+1)))
    U2    <- crossprod(z/v, B)%*%U
    if(v_prior_fh$type == "GIG"){
      bias <- sqrt(w)*bet*mean(delta)
      s2   <- scale*w*mean(delta) + w*bet^2*var(delta)
    }

    if(k >= nburn & ((k-nburn) %% thin) == 0){
      ss[kk]     <- mean((y-yhat)^2)
      M_mc[kk]   <- M
      lam_mc[kk] <- lam
      tau_mc[kk] <- tau
      w_mc[kk]   <- w
      bet_mc[kk] <- bet
      v_mc[kk,]  <- v
      basis_mc[[kk]] <- basis_index
      a_mc[[kk]] <- a

      if(v_prior_fh$type == "GIG"){
        bias_mc[kk] <- bias
        s2_mc[kk]  <- s2
      }

      # Fay-herriot stuff
      gamma_v_mc[kk] <- gamma_v
      gamma_w_mc[kk] <- gamma_w

      kk <- kk + 1
    }

    if(verbose & k%%1000 == 0){
      pr<-c('MCMC iteration',k,myTimestamp(),'nbasis:',M)
      cat(pr,'\n')
    }

  } #END MCMC ITERATIONS
  #browser()
  obj <- list(nbasis=M_mc, w=w_mc, v=v_mc,
              tau=tau_mc, lam=lam_mc, beta=bet_mc,
              a=a_mc, basis=basis_mc, lookup=lookup,
              cnt1=cnt1, cnt2=cnt2, ss=ss,
              v_prior_fh=v_prior_fh, M=M_mc, scale=scale, s2=s2_mc, bias=bias_mc,
              gamma_v=gamma_v_mc, gamma_w=gamma_w_mc,
              X=X, y=y)
  class(obj) <- "gbass"
  return(obj)
} #END FUNCTION


## NEED TO UPDATE ALL OF THESE PREDICT FUNCTIONS
#' @export
predict.fhbass <- function(object, newdata=NULL, mcmc.use=NULL){
  if(is.null(newdata)){
    newdata <- object$X
  }
  N <- nrow(newdata)
  if(is.null(mcmc.use)) mcmc.use <- seq_along(object$M)
  res <- matrix(NA, ncol=length(mcmc.use), nrow=N)
  tX <- t(newdata)
  for(i in mcmc.use){
    basis_curr <- object$lookup[unlist(object$basis[i])]
    B_curr <- matrix(1, nrow=N, ncol=length(basis_curr)+1)
    for(m in seq_along(basis_curr)){
      thetaB <- basis_curr[[m]]
      B_curr[,m+1] <- makeBasis(thetaB$s,thetaB$u,thetaB$t,tX,1)
      #for(j in 1:thetaB$J){
      #  B_curr[,m+1] <- B_curr[,m+1]*pmax(0, thetaB$s[j]*(newdata[,thetaB$u[j]]-thetaB$t[j]))
      #}
    }
    res[,i] <- B_curr%*%object$a[[i]]
  }
  return(t(res))
}

mh_update_gamma <- function(gamma, logdens, prop_sigma, scale, cnt) {
  log_gamma_prop <- log(gamma) + rnorm(1, 0, prop_sigma)
  gamma_prop <- exp(log_gamma_prop)

  log_alpha <- logdens(gamma_prop) - logdens(gamma)

  if (log(runif(1)) < log_alpha) {
    cnt <- cnt + 1
    list(gamma = gamma_prop, cnt = cnt, accept = TRUE)
  } else {
    list(gamma = gamma, cnt = cnt, accept = FALSE)
  }
}

logdens_gamma_v <- function(gv) {
  a_v <- gv / delta
  b_v <- gv * delta
  sum(dgig(v, p = v_prior_fh$p - 1/2, a = a_v, b = b_v, log = TRUE)) +
    dcauchy(gv, scale = gamma_v_scale, log = TRUE)
}

logdens_gamma_w <- function(gw) {
  a_w <- gw * (1 + 1 / gamma_v)
  b_w <- gw / (1 + 1 / gamma_v)
  dgig(w, p = w_prior$p - (N + M + 1) / 2, a = a_w, b = b_w, log = TRUE) +
    dcauchy(gw, scale = gamma_w_scale, log = TRUE)
}

dgig <- function(x, p, a, b, log = FALSE) {
  if (any(x <= 0)) stop("x must be positive")

  sqrt_ab <- sqrt(a * b)
  log_c <- (p - 1) * log(x) - 0.5 * (a * x + b / x)
  log_norm <- -log(2) - log(besselK(sqrt_ab, p)) + (p * log(sqrt(b / a)))

  out <- log_c - log_norm
  if (log) return(out)
  else return(exp(out))
}
