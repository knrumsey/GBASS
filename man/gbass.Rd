% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gbass.R
\name{gbass}
\alias{gbass}
\title{Generalized Bayesian MARS}
\usage{
gbass(
  X,
  y,
  w_prior = list(type = "GIG", p = 0, a = 0, b = 0),
  v_prior = list(type = "GIG", p = -15, a = 0, b = 30),
  maxInt = 3,
  maxBasis = 1000,
  npart = NULL,
  nmcmc = 10000,
  nburn = 9001,
  thin = 1,
  moveProbs = rep(1/3, 3),
  a_tau = 1/2,
  b_tau = NULL,
  a_lambda = 1,
  b_lambda = 1,
  m_beta = 0,
  s_beta = 0,
  scale = 1,
  Iw0 = rep(1, maxInt),
  Zw0 = rep(1, ncol(X)),
  verbose = TRUE
)
}
\arguments{
\item{X}{an Nxp matrix of predictor variables.}

\item{y}{an Nx1 matrix (or vector) of response values.}

\item{w_prior}{a named list specifying the prior for the global variance component. See details.}

\item{v_prior}{a named list specifying the (shared) prior for the local variance components. See details.}

\item{maxInt}{integer for maximum degree of interaction in spline basis functions. Defaults to the number of predictors, which could result in overfitting.}

\item{maxBasis}{maximum number of basis functions. This should probably only be altered if you run out of memory.}

\item{npart}{of non-zero points in a basis function. If the response is functional, this refers only to the portion of the basis function coming from the non-functional predictors. Defaults to 20 or 0.1 times the number of observations, whichever is smaller.}

\item{nmcmc}{number of mcmc iterations}

\item{nburn}{number of burn-in samples}

\item{thin}{thinning for mcmc}

\item{moveProbs}{a vector defining the probabilities for (i) birth (ii) death and (iii) mutation. Default is rep(1/3,3).}

\item{a_tau}{prior for tau}

\item{b_tau}{prior for tyau}

\item{a_lambda}{prior for lambda}

\item{b_lambda}{prior for lambda}

\item{m_beta}{prior for beta}

\item{s_beta}{prior for beta}

\item{scale}{fixed variance parameter. default is one.}

\item{Iw0}{vector of nominal weights for degree of interaction, used in generating candidate basis functions. Should have length equal to Jmax and have positive entries.}

\item{Zw0}{vector of nominal weights for variable selection, used in generating candidate basis functions. Should have length equal to ncol(X) and have positive entries.}

\item{verbose}{Logical. Should gbass print completion status? Default TRUE}
}
\value{
The returned value is a named list with components for each of the MCMC parameters. The acceptance rates for each move type is returned. If applicable, we also return acceptance rates for w and the v_i.
}
\description{
A function for Bayesian non-linear regression under various likelihood functions.
}
\details{
Currently, the prior for w and v_i must belong to the class of Generalized inverse Gaussian (GIG) or Generalized Beta Prime (GBP) priors. The list should have the following named fields
\enumerate{
\item type. either "GIG" or "GBP".
\item p, a, b. Hyperparameters for the prior. p,a,b > 0 for GBP. See ?rgig2 for details on GIG parameters.
\item prop_sigma. The proposal standard deviation for Metropolis-Hastings. Only needed if type="GBP" or if type="GIG" and beta is not fixed at zero.
\item lb. An optional lower bound which truncates the prior for w. This argument is ignored when specified for v_prior.
}
The build_prior function can be used to construct these priors.
}
\note{
Some comments about current deficiencies in the code.
\enumerate{
\item basis function parameters are stored as lists.
\item burn-in and thinning is not implemented intelligently.
\item continuous uniform prior for knot locations.
\item assumes a ridge prior for basis coefficients.
}
}
\examples{
n <- 100 #Number of observations
p <- 4   #Number of variables (beyond p = 2, variables are inert)
X <- matrix(runif(n*p), nrow=n)
y <- apply(X, 1, ff1)
mod <- gbass(X, y, nmcmc=1000, nburn=901, thin=2)

}
