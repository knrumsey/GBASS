test_that("GBASS with t likelihood", {
  cat('simple test function')

  set.seed(1)
  ff1 <- function(x) 10.391*((x[1]-0.4)*(x[2]-0.6) + 0.36)
  X <- lhs::maximinLHS(300, 2)
  y <- 6*apply(X, 1, ff1)
  mod <- gbass(X, y, nmcmc=2000, nburn=1001)

  X2 <- lhs::maximinLHS(100, 2)
  y2 <- 6*apply(X2, 1, ff1)
  yhat <- apply(predict(mod, X2)
                , 2, mean)
  d1 <- sqrt(mean((y2-yhat)^2))
  expect_that(d1, is_less_than(sd(y)/3))
})
