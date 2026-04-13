test_that("NW Regression with GBASS", {
  cat('simple NW regression example')

  set.seed(1)
  ff1 <- function(x) 10.391*((x[1]-0.4)*(x[2]-0.6) + 0.36)
  X <- lhs::maximinLHS(300, 2)
  y <- 6*apply(X, 1, ff1) + (rgamma(300, 1.5, 1) - 1.5)
  mod <- nwbass(X, y,
               s_beta=10, m_gamma=100, s_gamma=30,
               nmcmc=2000, nburn=1001)

  X2 <- lhs::maximinLHS(100, 2)
  y2 <- 6*apply(X2, 1, ff1) #+ (rgamma(100, 1.5, 1) - 1.5)
  yhat <- apply(predict(mod, X2)
                ,2, mean)
  d1 <- sqrt(mean((y2-yhat)^2))
  expect_that(d1, is_less_than(sd(y2)/3))
})


