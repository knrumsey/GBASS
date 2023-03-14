#' Predict function for class "gbass"
#'
#' This function generates posterior predictions using predictors Xnew and a model of class "gbmars"
#'
#' @param object an object of class "gbass" created by the gbass, tbass or hbass functions
#' @param newdata a matrix of predictor variables with ncol(newdata) == ncol(X)
#' @param mcmc.use a vector subsetting which posterior samples to use for prediction. Default is to retain all samples.
#' @details Returns a matrix of posterior predictions.
#'
#' @export
predict.gbass <- function(object, newdata=NULL, mcmc.use=NULL){
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


#' #' Predict function for class "gbmars" (Depreciated)
#' #'
#' #' This function generates posterior predictions using predictors Xnew and a model of class "gbmars"
#' #'
#' #' @param model an object of class "gbmars" created by the gbmars() function.
#' #' @param Xnew a matrix of predictor variables with ncol(Xnew) == ncol(X)
#' #' @param burn burn-in
#' #' @param thin thinning
#' #' @details add stuff here
#' #' @examples
#' #' foo <- 1 + 1
#' #'
#' predict.gbmars <- function(model, Xnew, burn=1, thin=1){
#'   N <- nrow(Xnew)
#'   p <- ncol(Xnew)
#'   K <- model$nsamp
#'   index <- seq((burn+1), model$iter, by=thin)
#'   y <- matrix(NA, nrow=length(index), ncol=N)
#'   cnt <- 1
#'   for(k in index){
#'     a_temp <- model$a[[k]]
#'     B_temp <- buildB(model, Xnew, k)
#'     y[cnt,] <- as.numeric(B_temp%*%a_temp)
#'     cnt <- cnt + 1
#'   }
#'   return(y)
#' }


#Function to build basis column
buildB <- function(model, Xnew, k){
  N_new <- nrow(Xnew)
  B <- matrix(1, nrow=N_new)
  k_curr <- model$k[[k]]
  M <- length(k_curr)
  if(M < 1) return(B)
  for(m in 1:M){
    J_cand <- model$J[[k_curr[m]]]
    u_cand <- model$u[[k_curr[m]]]
    s_cand <- model$s[[k_curr[m]]]
    t_cand <- model$t[[k_curr[m]]]
    B_new <- rep(1, N_new)
    for(j in 1:J_cand){
      B_new <- B_new * pmax(0, s_cand[j]*(Xnew[,u_cand[j]] - t_cand[j]))
    }
    B <- cbind(B, B_new)
  }
  return(B)
}

