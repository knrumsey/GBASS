#' Generalized Bayesian MARS with a Normal Wald Likelihood
#'
#' Fits a generalized BMARS model with Normal-Wald likelihood. General purpose regression with flexible error distribution (unimodal).
#'
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
#' @param lag_beta number of steps for which beta is fixed to m_beta (usually zero)
#' @param m_gamma prior for gamma
#' @param s_gamma prior for gamma
#' @param vsig proposal sd for v_i
#' @param scale fixed variance parameter. default is one.
#' @param Iw0 vector of nominal weights for degree of interaction, used in generating candidate basis functions. Should have length equal to Jmax and have positive entries.
#' @param Zw0 vector of nominal weights for variable selection, used in generating candidate basis functions. Should have length equal to ncol(X) and have positive entries.
#' @param verbose Logical. Should gbass print completion status? Default TRUE
#' @details See details for gbass(). Additionally, a normal prior must be specified for asymmetry parameter (beta) and tail heaviness parameter (gamma)
#' See also nw_gamma_prior() and nw_triangle()
#'
#' @return The returned value is a named list with components for each of the MCMC parameters. The acceptance rates for each move type is returned. If applicable, we also return acceptance rates for w and the v_i.
#' @import Matrix
#' @export
#' @examples
#' n <- 100 #Number of observations
#' p <- 4   #Number of variables (beyond p = 2, variables are inert)
#' X <- matrix(runif(n*p), nrow=n)
#' y <- apply(X, 1, ff1)
#' y <- y + rgamma(n, 1, 0.1)/4
#' gamma_prior <- nw_gamma_prior(p1=0.5, p2=1e-4) #Get reccomended prior for gamma
#' mod <- nwbass(X, y, m_beta=0, s_beta=30, m_gamma=gamma_prior[1], s_gamma=gamma_prior[2],
#'               nmcmc=1000, nburn=901)
#' nw_triangle(mod)
#'
nwbass <- function(X, y,
                  w_prior=list(type="GIG", p=0, a=0, b=0, prop_sigma=0.2),
                  maxInt=3, maxBasis=1000, npart=NULL, nmcmc=10000, nburn=9001, thin=1,
                  moveProbs=rep(1/3,3),
                  a_tau=1/2, b_tau=NULL,
                  a_lambda=1, b_lambda=1,
                  m_beta=0, s_beta=100, lag_beta=round(nmcmc*0.5),
                  m_gamma=1, s_gamma=0, vsig=0.25,
                  scale=1,
                  Iw0 = rep(1, maxInt), Zw0 = rep(1, ncol(X)),
                  verbose=TRUE){

  v_prior <- list(type="GIG", p=-1/2, a=1, b=1, prop_sigma=vsig)
  if(nrow(X) != length(y)) stop("nrow(X) and length(y) should match")
  if(w_prior$a != 0) warning("Gibbs sampler assumes w_prior$a = 0")
  if(lag_beta > nburn){
    lag_beta <- nburn
    warning("lag_beta cannot exceed nburn. Setting lag_beta = nburn")
  }
  s_beta_hold <- s_beta
  s_beta <- 0
  N <- nrow(X)
  p <- ncol(X)
  nkeep <- length(seq(nburn, nmcmc, by=thin))
  if(is.null(npart)) npart <- min(20, 0.1*N)
  if(is.null(b_tau)) b_tau <- N/2
  if(is.null(w_prior$lb)) w_prior$lb <- 1/N
  if(is.null(v_prior$lb)) v_prior$lb <- 0

  #ALLOCATE STORAGE SPACE
  a_mc <- list()
  M_mc <- lam_mc <- tau_mc <- w_mc <- bet_mc <- gam_mc <- ss <- rep(NA, nkeep)
  v_mc   <- matrix(NA, nrow=nkeep, ncol=N)
  lookup <- basis_mc <- list()
  basis_index <- integer(0)
  nlook <- 1
  kk    <- 1

  #INITIALIZE PARAMETERS
  M   <- 0
  tau <- (b_tau+.1)/(a_tau+.1)
  gam <- 1
  v   <- rep(1, N)
  w   <- 2*var(y)/scale #Multiply by 2 tries to avoid collapsing to degenerate solution
  lam <- rgamma(1, a_lambda, b_lambda)
  bet <- rnorm(1, m_beta, s_beta)
  a   <- mean(y)
  z    <- y - bet*v*sqrt(w)
  B    <- matrix(1, nrow=N)
  Vinv <- Matrix::Diagonal(x=1/v)
  U    <- solve(symchol(t(B)%*%Vinv%*%B + scale/tau*Diagonal(M+1)))
  U2   <- t(z/v)%*%B%*%U

  cnt1 <- cnt2 <- rep(0, 3)
  if(w_prior$type == "GBP" || abs(bet) > 1e-9){
    cntw <- 0

  }
  if(v_prior$type == "GBP"){
    cntv <- 0
  }
  tX <- Matrix::t(X)
  if(verbose) cat("UPDATING\n==============================\n")
  for(k in 1:nmcmc){
    if(k == lag_beta) s_beta <- s_beta_hold
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
        U_cand  <- tryCatch(solve(symchol(t(B_cand)%*%Vinv%*%B_cand + scale/tau*Diagonal(M+2))), error=function(e) FALSE)
        if(isFALSE(U_cand)){
          log_accept_prob_B <- -Inf
        }else{
          U2_cand <- t(z/v)%*%B_cand%*%U_cand
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
      U_cand  <- tryCatch(solve(symchol(t(B_cand)%*%Vinv%*%B_cand + scale/tau*Diagonal(M))), error=function(e) FALSE)
      if(isFALSE(U_cand)){
        log_accept_prob_D <- -Inf
      }else{
        U2_cand <- t(z/v)%*%B_cand%*%U_cand
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
        U_cand  <- tryCatch(solve(symchol(t(B_cand)%*%Vinv%*%B_cand + scale/tau*Diagonal(M+1))), error=function(e) FALSE)
        if(isFALSE(U_cand)){
          log_accept_prob_M <- -Inf
        }else{
          U2_cand <- t(z/v)%*%B_cand%*%U_cand
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
    gam <- rnorm(1, (s_gamma^2*N+m_gamma)/(s_gamma^2*sum(v)+1), s_gamma/sqrt(s_gamma^2*sum(v)+1))
    if(w_prior$type == "GIG"){
      if(abs(bet) < 1e-9){
        w_cand <- rgig2(p=w_prior$p - (N+M+1)/2,
                        a=w_prior$a,
                        b=w_prior$b + sum(r^2/v)/scale + pen/tau)
        w <- ifelse(w_cand < w_prior$lb, w, w_cand)
      }else{
        #Metropolis Hastings (depreciated)
        #w_cand <- exp(log(w) + rnorm(1, 0, w_prior$prop_sigma))
        #log_alpha_w <- (w_prior$p - (N+M+1)/2)*(log(w_cand)-log(w)) -
        #  0.5*(w_prior$a*(w_cand-w) +
        #         (w_prior$b + sum(r^2/v)/scale + pen/tau)*(1/w_cand - 1/w) -
        #         2*bet/scale*sum(r)*(1/sqrt(w_cand)-1/sqrt(w))) +
        #  log(w_cand > w_prior$lb)
        #if(log_alpha_w > log(runif(1))){
        #  cntw <- cntw + 1
        #  w <- w_cand
        #}
        w_tform <- rmpon(1,
                         N+M+w_prior$p-1,
                         0.5*(w_prior$b + sum(r^2/v)/scale + pen/tau),
                         bet*sum(r)/scale/(sum(r^2/v)/scale + pen/tau))
        w <- w_tform^-2
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
    U    <- solve(symchol(t(B)%*%Vinv%*%B + scale/tau*Diagonal(M+1)))
    U2   <- t(z/v)%*%B%*%U
    v_prior$a <- gam^2       # THIS SEEMS LIKE A MAJOR PROBLEM. WHY DID I DO THIS AND WHY DIDNT I DOCUMENT IT ))))))::::

    if(k >= nburn & ((k-nburn) %% thin) == 0){
      ss[kk]     <- mean((y-yhat)^2)
      M_mc[kk]   <- M
      lam_mc[kk] <- lam
      tau_mc[kk] <- tau
      w_mc[kk]   <- w
      bet_mc[kk] <- bet
      gam_mc[kk] <- gam
      v_mc[kk,]  <- v
      basis_mc[[kk]] <- basis_index
      a_mc[[kk]] <- a
      kk <- kk + 1
    }

    if(verbose){
      if((k%%(nmcmc/10)) == 0){
        cat("***")
        #browser()
      }
    }

  } #END MCMC ITERATIONS
  if(verbose) cat("\n")
  #browser()
  obj <- list(M=M_mc, w=w_mc, v=v_mc, tau=tau_mc, lam=lam_mc, a=a_mc, bet=bet_mc, gam=gam_mc, basis=basis_mc, lookup=lookup,
              cnt1=cnt1, cnt2=cnt2, ss=ss, v_prior=v_prior)
  class(obj) <- "gbass"
  return(obj)
} #END FUNCTION



#' The NW Triangle Plot
#'
#' A function for plotting the NW Triangle for an object created with nwbass()
#'
#' @param obj an object returned by nwbass function
#' @return Creates the NW Triangle plot
#' @details Plots the posterior draws of the steepness and asymmetry parameter, defined as (1+gamma)^(-1/2) and beta/sqrt(gamma^2+beta^2)*steepness respectively.
#' These parameters are location and scale invariant. Normal (0,0) and Cauchy (0,1) occur as limiting cases.
#' @import Matrix
#' @export
#' @examples
#' #not yet
#'
nw_triangle <- function(obj, add=FALSE, ...){
  steepness <- (obj$gam+1)^(-1/2)
  asymmetry <- obj$bet/sqrt(obj$bet^2+obj$gam^2)*(steepness)
  if(!add){
    plot(NULL, xlim=c(-1,1), ylim=c(0, 1), xlab='Asymmetry', ylab='Steepness')
    segments(x0=c(-1, -1, 0), x1=c(1, 0, 1), y0=c(1,1,0), y1=c(1, 0, 1), lwd=3, col=adjustcolor('gray', alpha.f=0.5))
    #points(c(0,0), c(0, 1), cex=1.5)
  }
  points(asymmetry, steepness, ...)
}


#' Function to select prior for gamma
#'
#' A  function to select hyperparameters for gamma in terms of steepness parameter xi = (1+gamma)^(-1/2)
#'
#' @param q1 default 0.1
#' @param q2 default 0.9
#' @param p1 default 0.5
#' @param p2 default 0.05
#' @param par0 optional starting values
#' @param par0 optional  ridge penalty for optimization
#' @return hyperparameter values for prior of gamma
#' @details Plots the posterior draws of the steepness and asymmetry parameter, defined as (1+gamma)^(-1/2) and beta/sqrt(gamma^2+beta^2)*steepness respectively.
#' These parameters are location and scale invariant. Normal (0,0) and Cauchy (0,1) occur as limiting cases.
#' @import Matrix
#' @export
#' @examples
#' #not yet
nw_gamma_prior <- function(q1=0.1, q2=0.9, p1=0.5, p2=0.05, par0=NULL, lambda=0){
  logit <- function(z) log(z/(1-z))
  fff <- function(xx, q1,q2,p1,p2,lambda){
    lhs1 <- 1 + pnorm((1-q1^-2-xx[1])/xx[2]) - pnorm((q1^-2-1-xx[1])/xx[2])
    lhs2 <- pnorm((q2^-2-1-xx[1])/xx[2]) - pnorm((1-q2^-2-xx[1])/xx[2])
    (logit(lhs1)-logit(p1))^2 + (logit(lhs2)-logit(p2))^2 + lambda*sum(xx^2)
  }
  if(is.null(par0)) par0 <- c(100, 100)
  opt <- optim(par0, fff, method="Nelder-Mead", q1=q1, q2=q2, p1=p1, p2=p2, lambda=lambda, control=list(maxit=c(50000)))
  par <- opt$par
  #return(opt)
  #if(opt$value - lambda*sum(par^2) > 0.5 ) warning("Attained a minimum of ", opt$value - lambda*sum(par^2), ", optimization may not have converged")
  names(par) <- c("m_gamma", "s_gamma")
  return(par)
}




