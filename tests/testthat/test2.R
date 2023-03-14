test_that("Quantile regression with GBASS", {
  cat('simple quantile regression example')

  X <- lhs::maximinLHS(300, 2)
  y <- 6*apply(X, 1, ff1) + (rgamma(300, 1.5, 1) - 1.5)
  mod <- qbass(X, y, q=0.75,
               nmcmc=2000, nburn=1001)

  X2 <- lhs::maximinLHS(100, 2)
  y2 <- 6*apply(X2, 1, ff1) + (qgamma(0.75, 1.5, 1) - 1.5)
  yhat <- apply(predict(mod, X2)
                , 2, mean)
  d1 <- sqrt(mean((y2-yhat)^2))
  expect_that(d1, is_less_than(0.5))
})
