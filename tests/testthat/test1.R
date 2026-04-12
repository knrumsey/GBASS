test_that("GBASS with t likelihood", {
  cat('simple test function')

  ff1 <- function (x){
    1.3356 * (1.5 * (1 - x[1]) + exp(2 * x[1] - 1) * sin(3 * pi * (x[1] - 0.6)^2) + exp(3 * (x[2] - 0.5)) * sin(4 * pi * (x[2] - 0.9)^2))
  }
  X <- lhs::maximinLHS(300, 2)
  y <- apply(X, 1, ff1)
  mod <- gbass(X, y, nmcmc=2000, nburn=1001)
  yhat <- apply(predict(mod), 2, mean)
  d1 <- sqrt(mean((y-yhat)^2))
  expect_that(d1, is_less_than(0.4))
})
