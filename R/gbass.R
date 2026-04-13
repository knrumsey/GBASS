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
  if(nrow(X) != length(y)) stop("nrow(X) and length(y) should match")
  N <- nrow(X)
  p <- ncol(X)
  if(N < 4){
    warning("Trying to fit a model with less than 4 data points may fail.")
  }
  maxInt <- min(ncol(X), maxInt)
  keep_idx <- seq.int(from = nburn, to = nmcmc, by = thin)
  nkeep <- length(keep_idx)
  keep_iter <- rep(FALSE, nmcmc)
  keep_iter[keep_idx] <- TRUE
  if(is.null(npart)) npart <- min(20, 0.1*N)
  if(is.null(b_tau)) b_tau <- N/2
  if(is.null(w_prior$lower_bound)) w_prior$lower_bound <- var(y)/N
  if(is.null(v_prior$lower_bound)) v_prior$lower_bound <- 0
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
  U    <- solve(symchol(crossprod(B, Vinv)%*%B + scale/tau*Diagonal(M+1)))
  U2   <- crossprod(z/v, B)%*%U
  if(v_prior$type == "GIG"){
    mu_v <- mu_gig(v_prior$p, v_prior$a, v_prior$b)
    s2_v <- var_gig(v_prior$p, v_prior$a, v_prior$b)
    bias <- sqrt(w)*bet*mu_v
    s2   <- scale*w*mu_v + w*bet^2*s2_v
  }

  cnt1 <- cnt2 <- rep(0, 3)
  if(w_prior$type == "GBP"){
    if(is.null(w_prior$prop_sigma)){
      w_prior$prop_sigma <- min(max(var(y) / scale / 2, 1e-3), 5)
    }
    if(is.null(w_prior$adapt)){
      w_prior$adapt <- TRUE
    }

    # Initialize parameters for adaptive metropolis
    if(w_prior$adapt){
      if(is.null(w_prior$adapt_delay)){
        w_prior$adapt_delay <- floor(nburn / 10)
      }
      if(is.null(w_prior$adapt_thin)){
        w_prior$adapt_thin <- 1
      }

      w_prior$count <- 1
      w_prior$sum1 <- log(w)
      w_prior$welford <- 0
      adapt_iter <- seq_len(nmcmc) %in% seq(1, nburn, by=w_prior$adapt_thin)
    }

    # Track acceptance rate
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

    # Sample w
    rss_over_v <- sum(r^2 / v) / scale
    sum_r <- sum(r)
    quad_term <- w_prior$b + rss_over_v + pen / tau
    shape_shift <- (N + M + 1) / 2

    if (w_prior$type == "GIG") {
      if (abs(bet) < 1e-9) {
        w_cand <- rgig2(
          p = w_prior$p - shape_shift,
          a = w_prior$a,
          b = quad_term
        )

        if (w_cand >= w_prior$lower_bound) {
          w <- w_cand
        }
      } else {
        w_tform <- rMHN(
          n = 1,
          alpha = N + M - 2 * w_prior$p + 1,
          beta = 0.5 * quad_term,
          gamma = bet * sum_r / scale
        )

        w <- w_tform^(-2)
      }
    } else {
      # Adaptive Metropolis Hastings
      w_cand <- exp(log(w) + rnorm(1, 0, w_prior$prop_sigma))

      log_alpha_w <-
        (w_prior$p * w_prior$a - shape_shift) * (log(w_cand) - log(w)) -
        0.5 * (rss_over_v + pen / tau) * (1 / w_cand - 1 / w) -
        (w_prior$a + w_prior$b) *
        (log(1 + w_cand^w_prior$p) - log(1 + w^w_prior$p)) +
        (bet / scale) * sum_r * (1 / sqrt(w_cand) - 1 / sqrt(w)) +
        log(w_cand > w_prior$lower_bound)

      if (log_alpha_w > log(runif(1))) {
        cntw <- cntw + 1
        w <- w_cand
      }

      if(w_prior$adapt){
        if(adapt_iter[k]){
          theta <- log(w)

          n_adapt  <- w_prior$count + 1
          S1_adapt <- w_prior$sum1 + theta

          w_prior$count <- n_adapt
          w_prior$sum1  <- S1_adapt

          xbar_n <- S1_adapt / n_adapt
          xbar_n_1 <- (S1_adapt - theta) / (n_adapt - 1)
          w_prior$welford <- w_prior$welford + (theta - xbar_n_1) * (theta - xbar_n)

          if(k > w_prior$adapt_delay){
            var_adapt <- w_prior$welford / (n_adapt - 1)
            ps_adapt <- sqrt(2.38^2 * var_adapt + 1e-8)
            w_prior$prop_sigma <- min(max(ps_adapt, 1e-3), 5)
          }
        }
      }
    }

    # Sample v_i's
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
      #browser()
      bias <- sqrt(w)*bet*mu_v
      s2   <- scale*w*mu_v + w*bet^2*s2_v
    }


    if(keep_iter[k]){
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



#' Predict method for GBASS objects
#'
#' Returns posterior draws of either:
#'   (i) the linear predictor / mean surface, or
#'   (ii) the full posterior predictive distribution.
#'
#' @param object an object of class "gbass" (including subclasses like
#'   "tbass", "qbass", "nwbass")
#' @param newdata a matrix of predictor variables. Defaults to training inputs.
#' @param mcmc.use optional vector indicating which posterior draws to use.
#' @param predictive logical. If TRUE, return posterior predictive draws.
#'   If FALSE, return draws of the linear predictor.
#' @param bias_correct logical. Ignored unless predictive = FALSE.
#'   If TRUE, return the posterior mean response rather than just the linear predictor.
#'
#' @details
#' If predictive = FALSE and bias_correct = FALSE, this returns draws of B(x)a.
#'
#' If predictive = FALSE and bias_correct = TRUE, this returns draws of
#'   B(x)a + E(error | posterior draw),
#' i.e. the mean response under the fitted GBASS error model.
#'
#' If predictive = TRUE, this returns posterior predictive draws by simulating
#' a fresh latent local variance v_new and Gaussian draw for each posterior sample.
#'
#' For qbass objects, bias_correct = TRUE is usually not what you want, because
#' qbass is typically being used for quantile regression rather than mean regression.
#'
#' Currently, posterior predictive draws are implemented for GIG-based models.
#' If the fitted object uses a GBP prior for v, predictive = TRUE will stop.
#'
#' @return a matrix with rows corresponding to posterior draws and columns
#'   corresponding to rows of newdata.
#'
#' @export
predict.gbass <- function(object, newdata = NULL, mcmc.use = NULL,
                          predictive = TRUE, bias_correct = FALSE) {

  if (is.null(newdata)) {
    newdata <- object$X
  }
  if (is.null(dim(newdata))) {
    newdata <- matrix(newdata, ncol = ncol(object$X))
  }

  N <- nrow(newdata)

  if (is.null(mcmc.use)) {
    if (!is.null(object$M)) {
      mcmc.use <- seq_along(object$M)
    } else if (!is.null(object$nbasis)) {
      mcmc.use <- seq_along(object$nbasis)
    } else {
      stop("Could not determine available MCMC draws from object.")
    }
  }

  if (predictive && isTRUE(bias_correct)) {
    warning("bias_correct is ignored when predictive = TRUE.")
  }

  if (inherits(object, "qbass") && !predictive && isTRUE(bias_correct)) {
    warning("bias_correct=TRUE is usually not what you want for qbass/quantile regression.")
  }

  # small helper to build basis matrix at newdata for one posterior draw
  build_basis_gbass <- function(object, draw_index, newdata) {
    tX <- t(newdata)
    basis_curr <- object$lookup[unlist(object$basis[[draw_index]])]
    B_curr <- matrix(1, nrow = nrow(newdata), ncol = length(basis_curr) + 1)

    if (length(basis_curr) > 0) {
      for (m in seq_along(basis_curr)) {
        thetaB <- basis_curr[[m]]
        B_curr[, m + 1] <- makeBasis(thetaB$s, thetaB$u, thetaB$t, tX, 1)
      }
    }
    B_curr
  }

  # helper: draw from the prior for a new local variance factor v_new
  draw_v_new <- function(object, draw_index, n) {

    # nwbass: assume intended NW model from paper:
    # v ~ GIG(-1/2, gamma^2, 1)
    if (inherits(object, "nwbass")) {
      if (is.null(object$gamma)) {
        stop("nwbass object does not contain gamma draws.")
      }
      gam <- object$gamma[draw_index]
      return(replicate(n, rgig2(-1/2, gam^2, 1)))
    }

    # generic gbass / tbass / qbass path
    vp <- object$v_prior

    if (is.null(vp$type)) {
      stop("object$v_prior$type is missing.")
    }

    if (vp$type == "GIG") {
      return(replicate(n, rgig2(vp$p, vp$a, vp$b)))
    }

    stop("predictive=TRUE is currently implemented only for GIG-based v_prior objects.")
  }

  # helper: mean of v under the prior, for bias correction
  mean_v_prior <- function(object, draw_index) {

    if (inherits(object, "nwbass")) {
      if (is.null(object$gamma)) {
        stop("nwbass object does not contain gamma draws.")
      }
      gam <- object$gamma[draw_index]
      # For GIG(-1/2, gamma^2, 1), E(v) = 1/gamma
      return(1 / gam)
    }

    vp <- object$v_prior

    if (vp$type == "GIG") {
      return(mu_gig(vp$p, vp$a, vp$b))
    }

    stop("bias_correct=TRUE currently requires a GIG-based v_prior (or nwbass).")
  }

  ndraw <- length(mcmc.use)
  res <- matrix(NA_real_, nrow = ndraw, ncol = N)

  jj <- 1L
  for (i in mcmc.use) {

    B_curr <- build_basis_gbass(object, i, newdata)
    eta <- as.numeric(B_curr %*% object$a[[i]])

    # draw-specific parameters
    w_i <- object$w[i]
    beta_i <- if (!is.null(object$beta)) object$beta[i] else 0
    scale_i <- if (!is.null(object$scale)) object$scale else 1

    if (!predictive) {
      if (isTRUE(bias_correct)) {
        mu_v_i <- mean_v_prior(object, i)
        res[jj, ] <- eta + sqrt(w_i) * beta_i * mu_v_i
      } else {
        res[jj, ] <- eta
      }
    } else {
      v_new <- draw_v_new(object, i, N)
      z_new <- rnorm(N)
      eps_new <- sqrt(w_i) * (beta_i * v_new + sqrt(scale_i * v_new) * z_new)
      res[jj, ] <- eta + eps_new
    }

    jj <- jj + 1L
  }

  return(res)
}


#' Plot function for class "gbass"
#'
#' This function plots an object of class "gbass"
#'
#' @param object an object of class "gbass" created by the gbass, tbass or hbass functions
#' @param newdata a matrix of predictor variables with ncol(newdata) == ncol(X)
#' @param mcmc.use a vector subsetting which posterior samples to use for prediction. Default is to retain all samples.
#' @details Returns a matrix of posterior predictions.
#'
#' @export
plot.gbass <- function(x, ...){
  if(!("gbass" %in% class(x)))
    stop("x must be an object of class gbass")
  par(mfrow=c(2,2))
  plot(x$M, type='l',
       xlab="MCMC iteration (post-burn)",
       ylab="number of basis functions")
  plot(x$w, type='l',
       xlab="MCMC iteration (post-burn)",
       ylab="global error variance (w)")
  yhat <- apply(predict(x), 2, mean)
  plot(x$y, yhat,
       xlab="observed",
       ylab="posterior mean",
       main="Training Fit")
  resid <- x$y - yhat
  hist(resid,
       main="Posterior mean residuals",
       xlab="residuals",
       ylab="Density",
       breaks="FD",
       freq=FALSE)
  curve(dnorm(x, mean(resid), sd(resid)),
        add=TRUE,
        col='red')
  par(mfrow=c(1,1))
}
