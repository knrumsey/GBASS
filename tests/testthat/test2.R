test_that("simple polynomial example with Gaussian measure", {
  cat('simple polynomial example test with Gaussian measure')

  f <- function(x){
    x[1]^2 + x[1]*x[2] + x[2]^3/9
  }

  # Sim Study Parameters
  N <- 1000
  p <- 3

  # Get true value of C matrix (with monte carlo)
  measure <- function() rnorm(p, 0.5, 0.1)
  Cmc <- C_mc(f, measure, nmc=1e5)

  X <- matrix(runif((N-2) * p), nrow=N-2, ncol=p)
  Yf <- apply(X, 1, f)

  mod_bass <- BASS::bass(X, Yf, verbose=FALSE)
  pr <- list()
  pr[[1]]  <- list(dist="normal", trunc=c(-Inf, Inf), mean=0.5, sd=0.1, weights=1)
  pr[[2]]  <- list(dist="normal", trunc=c(-Inf, Inf), mean=0.5, sd=0.1, weights=1)
  pr[[3]]  <- list(dist="normal", trunc=c(-Inf, Inf), mean=0.5, sd=0.1, weights=1)

  Cbass <- C_bass(mod_bass, prior=pr)

  d1 <- sum(abs(Cbass - Cmc))
  expect_that(d1, is_less_than(0.04))
})
