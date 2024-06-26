% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nwbass.R
\name{nw_est_mom}
\alias{nw_est_mom}
\title{Method of Moments for Normal-Wald}
\usage{
nw_est_mom(
  data = NULL,
  stats = NULL,
  mu = NA,
  delta = NA,
  beta = NA,
  alpha = NA,
  triangle = FALSE,
  ...
)
}
\arguments{
\item{data}{the data vector}

\item{stats}{a vector of length 4 containing mean, variance, skew, kurtosis. Ignored when data is provided.}

\item{mu}{location parameter, if fixed}

\item{delta}{scale parameter, if fixed}

\item{beta}{skewness parameter, if fixed}

\item{alpha}{tail heaviness parameter, if fixed}

\item{triangle}{logical. When TRUE, only the steepness and asymmetry values are returned.}

\item{...}{additional parameters passed to optim.}
}
\value{
estimated parameters
}
\description{
A function to estimate parameters of the Normal-Wald distribution
}
\details{
Method of moments estimators for NW parameters. The stats vector can contain NA values when parameters are fixed. If the mean is to be estimated, then \code{stats[1]} must be provided.
}
\examples{
n <- 500
y <- rgamma(n, 3, 1.5) + rlnorm(n, 1, 0.5)
skew <- mean(((y-mean(y))/sd(y))^3) # Sample skewness
kurt <- mean(((y-mean(y))/sd(y))^4) # Sample kurtosis

nw_est_mom(stats=c(NA, NA, skew, kurt), mu=0, delta=1, triangle=TRUE)
nw_est_mom(stats=c(NA, var(y), skew, kurt), mu=0, triangle=TRUE)
}
