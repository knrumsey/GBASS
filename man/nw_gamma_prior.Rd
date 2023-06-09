% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nwbass.R
\name{nw_gamma_prior}
\alias{nw_gamma_prior}
\title{Function to select prior for gamma}
\usage{
nw_gamma_prior(
  q1 = 0.1,
  q2 = 0.9,
  p1 = 0.5,
  p2 = 0.05,
  par0 = NULL,
  lambda = 0
)
}
\arguments{
\item{q1}{default 0.1}

\item{q2}{default 0.9}

\item{p1}{default 0.5}

\item{p2}{default 0.05}

\item{par0}{optional  ridge penalty for optimization}
}
\value{
hyperparameter values for prior of gamma
}
\description{
A  function to select hyperparameters for gamma in terms of steepness parameter xi = (1+gamma)^(-1/2)
}
\details{
Plots the posterior draws of the steepness and asymmetry parameter, defined as (1+gamma)^(-1/2) and beta/sqrt(gamma^2+beta^2)*steepness respectively.
These parameters are location and scale invariant. Normal (0,0) and Cauchy (0,1) occur as limiting cases.
}
\examples{
#not yet
}
