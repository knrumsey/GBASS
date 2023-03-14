test_that("simple polynomial example for bivariate concordance", {
  cat('simple polynomial example test: bivariate concordance check')

  f <- function(x){
    x[1]^2 + x[1]*x[2]
  }
  g <- function(x){
    x[1]^2 + x[1]*x[2] + x[2]^3/9
  }
  Ctrue <- matrix(c(8/3, 10/9, 11/12, 7/18), nrow=2, byrow=2) #Worked this out on paper

  #Monte Carlo
  #C0 <- conc_analysis(f, g, 2)$C$Cfg
  #C0 <- Cfg_mc(f, g, measure=2)

  #BASS
  X <- matrix(runif(1500), ncol=3)
  X <- rbind(X, rep(0, 3))
  X <- rbind(X, rep(1, 3))
  Yf <- apply(X, 1, f)
  Yg <- apply(X, 1, g)

  mod1  <- BASS::bass(X, Yf)
  mod2 <- BASS::bass(X, Yg)

  C1_list <- Cfg_bass(mod1, mod2, mcmc.use=seq(1, 1000, by=10))
  C1 <- matrix(0, nrow=3, ncol=3)
  for(i in 1:100){
    C1 <- C1 + C1_list[[i]]/100
  }

  dd <- sum((C1[1:2, 1:2] - Ctrue)^2)
  expect_that(dd, is_less_than(0.01))
})
