test_that("GBASS with t likelihood", {
  cat('simple test function')

  X <- lhs::maximinLHS(300, 2)
  y <- apply(X, 1, ff1)
  mod <- gbass(X, y, nmcmc=2000, nburn=1001)
  yhat <- apply(predict(mod), 2, mean)
  d1 <- sqrt(mean((y-yhat)^2))
  expect_that(d1, is_less_than(0.4))
})
