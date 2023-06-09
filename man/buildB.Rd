% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict.R
\name{buildB}
\alias{buildB}
\title{#' Predict function for class "gbmars" (Depreciated)
#'
#' This function generates posterior predictions using predictors Xnew and a model of class "gbmars"
#'
#' @param model an object of class "gbmars" created by the gbmars() function.
#' @param Xnew a matrix of predictor variables with ncol(Xnew) == ncol(X)
#' @param burn burn-in
#' @param thin thinning
#' @details add stuff here
#' @examples
#' foo <- 1 + 1
#'
predict.gbmars <- function(model, Xnew, burn=1, thin=1){
N <- nrow(Xnew)
p <- ncol(Xnew)
K <- model$nsamp
index <- seq((burn+1), model$iter, by=thin)
y <- matrix(NA, nrow=length(index), ncol=N)
cnt <- 1
for(k in index){
a_temp <- model$a[\link{k}]
B_temp <- buildB(model, Xnew, k)
y\link{cnt,} <- as.numeric(B_temp\%*\%a_temp)
cnt <- cnt + 1
}
return(y)
}}
\usage{
buildB(model, Xnew, k)
}
\description{
#' Predict function for class "gbmars" (Depreciated)
#'
#' This function generates posterior predictions using predictors Xnew and a model of class "gbmars"
#'
#' @param model an object of class "gbmars" created by the gbmars() function.
#' @param Xnew a matrix of predictor variables with ncol(Xnew) == ncol(X)
#' @param burn burn-in
#' @param thin thinning
#' @details add stuff here
#' @examples
#' foo <- 1 + 1
#'
predict.gbmars <- function(model, Xnew, burn=1, thin=1){
N <- nrow(Xnew)
p <- ncol(Xnew)
K <- model$nsamp
index <- seq((burn+1), model$iter, by=thin)
y <- matrix(NA, nrow=length(index), ncol=N)
cnt <- 1
for(k in index){
a_temp <- model$a[\link{k}]
B_temp <- buildB(model, Xnew, k)
y\link{cnt,} <- as.numeric(B_temp\%*\%a_temp)
cnt <- cnt + 1
}
return(y)
}
}
