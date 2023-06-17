#' Generalized Bayesian MARS
#'
#' A function for Bayesian non-linear regression under various likelihood functions.
#'
#' @param X an Nxp matrix of predictor variables.
#' @param y an Nx1 matrix (or vector) of response values.
#' @param w_prior a named list specifying the prior for the global variance component. See details.
#' @param v_prior a named list specifying the (shared) prior for the local variance components. See details.
#' @param maxInt integer for maximum degree of interaction in spline basis functions. Defaults to the number of predictors, which could result in overfitting.
#' @param maxBasis maximum number of basis functions. This should probably only be altered if you run out of memory.
#' @param npart of non-zero points in a basis function. If the response is functional, this refers only to the portion of the basis function coming from the non-functional predictors. Defaults to 20 or 0.1 times the number of observations, whichever is smaller.
#' @param nmcmc number of mcmc iterations
#' @param nburn number of burn-in samples
#' @param thin thinning for mcmc
#' @param moveProbs a vector defining the probabilities for (i) birth (ii) death and (iii) mutation. Default is rep(1/3,3).
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
#' @details Currently, the prior for w and v_i must belong to the class of Generalized inverse Gaussian (GIG) or Generalized Beta Prime (GBP) priors. The list should have the following named fields
#' \enumerate{
#'    \item type. either "GIG" or "GBP".
#'    \item p, a, b. Hyperparameters for the prior. p,a,b > 0 for GBP. See ?rgig2 for details on GIG parameters.
#'    \item prop_sigma. The proposal standard deviation for Metropolis-Hastings. Only needed if type="GBP" or if type="GIG" and beta is not fixed at zero.
#'    \item lb. An optional lower bound which truncates the prior for w. This argument is ignored when specified for v_prior.
#' }
#' The build_prior function can be used to construct these priors.
#'
#' @return The returned value is a named list with components for each of the MCMC parameters. The acceptance rates for each move type is returned. If applicable, we also return acceptance rates for w and the v_i.
#' @note Some comments about current deficiencies in the code.
#' \enumerate{
#'    \item basis function parameters are stored as lists.
#'    \item burn-in and thinning is not implemented intelligently.
#'    \item continuous uniform prior for knot locations.
#'    \item assumes a ridge prior for basis coefficients.
#' }
#' @import Matrix
#' @export
#' @examples
#' n <- 100 #Number of observations
#' p <- 4   #Number of variables (beyond p = 2, variables are inert)
#' X <- matrix(runif(n*p), nrow=n)
#' y <- apply(X, 1, ff1)
#' mod <- gbass(X, y, nmcmc=1000, nburn=901, thin=2)
#'
#' @export
#'
gbass <- function(X, y,
                  w_prior=list(type="GIG", p=0, a=0, b=0),
                  v_prior=list(type="GIG", p=-15, a=0, b=30),
                  maxInt=3, maxBasis=1000, npart=NULL, nmcmc=10000, nburn=9001, thin=1,
                  moveProbs=rep(1/3,3),
                  a_tau=1/2, b_tau=NULL,
                  a_lambda=1, b_lambda=1,
                  m_beta=0, s_beta=0,
                  scale=1,
                  Iw0 = rep(1, maxInt), Zw0 = rep(1, ncol(X)),
                  verbose=TRUE){

  if(is.null(dim(X))) X <- matrix(X, ncol=1)
  if(max(X) > 1 | min(X) < 0) warning("Found values of X outside of (0, 1).")
  if(nrow(X) != length(y)) stop("nrow(X) and length(y) should match")
  N <- nrow(X)
  p <- ncol(X)
  maxInt <- min(ncol(X), maxInt)
  nkeep <- length(seq(nburn, nmcmc, by=thin))
  if(is.null(npart)) npart <- min(20, 0.1*N)
  if(is.null(b_tau)) b_tau <- N/2
  if(is.null(w_prior$lb)) w_prior$lb <- 1/N
  if(is.null(v_prior$lb)) v_prior$lb <- 0
  if((m_beta != 0 | s_beta != 0) & is.null(w_prior$prop_sigma)){
    warning("Setting w_prior$prop_sigma = scale/var(y)/2. Consider specifying w_prior.")
    w_prior$prop_sigma <- var(y)/scale/2
  }
  if(min(X) < 0 | max(X) > 1) warning("Data out of range (0,1). Are you sure you want to do this?")

  #ALLOCATE STORAGE SPACE
  a_mc <- list()
  M_mc <- lam_mc <- tau_mc <- w_mc <- bet_mc <- ss <- bias_mc <- s2_mc <- rep(NA, nkeep)
  v_mc   <- matrix(NA, nrow=nkeep, ncol=N)
  lookup <- basis_mc <- list()
  basis_index <- integer(0)
  nlook <- 1
  kk    <- 1

  #INITIALIZE PARAMETERS
  M   <- 0
  tau <- (b_tau+.1)/(a_tau+.1)
  v   <- rep(1, N)
  w   <- 2*var(y)/scale #Multiply by 2 tries to avoid collapsing to degenerate solution
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
  if(v_prior$type == "GIG"){
    mu_v <- mu_gig(v_prior$p, v_prior$a, v_prior$b)
    s2_v <- var_gig(v_prior$p, v_prior$a, v_prior$b)
    bias <- sqrt(w)*bet*mu_v
    s2   <- scale*w*mu_v + w*bet^2*s2_v
  }

  cnt1 <- cnt2 <- rep(0, 3)
  if(w_prior$type == "GBP" || abs(bet) > 1e-9){
    cntw <- 0

  }
  if(v_prior$type == "GBP"){
    cntv <- 0
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
        w_cand <- exp(log(w) + rnorm(1, 0, w_prior$prop_sigma))
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
      w_cand <- exp(log(w) + rnorm(1, 0, w_prior$prop_sigma))
      log_alpha_w <- (w_prior$p*w_prior$a - (N+M+1)/2)*(log(w_cand)-log(w)) -
        (sum(r^2/v)/(2*scale) + pen/(2*tau))*(1/w_cand - 1/w) -
        (w_prior$a + w_prior$b)*(log(1+w_cand^w_prior$p)-log(1+w^w_prior$p)) +
        bet/scale*sum(r)*(1/sqrt(w_cand) - 1/sqrt(w)) +
        log(w_cand > w_prior$lb)
      if(log_alpha_w > log(runif(1))){
        cntw <- cntw + 1
        w <- w_cand
      }
    }
    if(v_prior$type == "GIG"){
      v <- rgig2.vec(p=v_prior$p-1/2,
                     a=v_prior$a+bet^2/scale,
                     b=as.numeric(v_prior$b + r^2/(w*scale)))
    }else{
      v <- rgbp.vec(as.numeric(v), v_prior, w, scale, as.numeric(r), bet)
    }

    z    <- y - bet*v*sqrt(w)
    Vinv <- Matrix::Diagonal(x=1/v)
    #U    <- solve(symchol(t(B)%*%Vinv%*%B + scale/tau*Diagonal(M+1)))
    #U2   <- t(z/v)%*%B%*%U
    U     <- solve(symchol(crossprod(B, Vinv)%*%B + scale/tau*Diagonal(M+1)))
    U2    <- crossprod(z/v, B)%*%U
    if(v_prior$type == "GIG"){
      bias <- sqrt(w)*bet*mu_v
      s2   <- scale*w*mu_v + w*bet^2*s2_v
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
      if(v_prior$type == "GIG"){
        bias_mc[kk] <- bias
        s2_mc[kk]  <- s2
      }
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
              v_prior=v_prior, M=M_mc, scale=scale, s2=s2_mc, bias=bias_mc,
              X=X, y=y)
  class(obj) <- "gbass"
  return(obj)
} #END FUNCTION



## makes negative values 0
pos<-function(vec){
  #replace(vec,vec<0,0)
  (abs(vec)+vec)/2
}

## largest value of basis function, assuming x's in [0,1], used for scaling
const<-function(signs,knots,degree){
  cc<-prod(((signs+1)/2 - signs*knots))^degree
  if(cc==0)
    return(1)
  return(cc)
} # since a product, can find for functional & categorical pieces separately, take product

## make basis function (from continuous variables)
makeBasis<-function(signs,vars,knots,datat,degree){
  cc<-const(signs,knots,degree)
  temp1<-pos(signs*(datat[vars,,drop=F]-knots))^degree # this only works for t(data)...
  if(length(vars)==1){
    return(c(temp1)/cc)
  } else{
    temp2<-1
    for(pp in 1:length(vars)){ # faster than apply
      temp2<-temp2*temp1[pp,]
    }
    return(temp2/cc)
  }
}
