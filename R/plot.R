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
