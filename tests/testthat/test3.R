test_that("NW Regression with GBASS", {
  cat('simple NW regression example')

  ff1 <- function (x){
    1.3356 * (1.5 * (1 - x[1]) + exp(2 * x[1] - 1) * sin(3 * pi * (x[1] - 0.6)^2) + exp(3 * (x[2] - 0.5)) * sin(4 * pi * (x[2] - 0.9)^2))
  }
  X <- lhs::maximinLHS(300, 2)
  y <- 6*apply(X, 1, ff1) + (rgamma(300, 1.5, 1) - 1.5)
  mod <- nwbass(X, y,
               nmcmc=2000, nburn=1001)

  X2 <- lhs::maximinLHS(100, 2)
  y2 <- 6*apply(X2, 1, ff1) + (rgamma(100, 1.5, 1) - 1.5)
  yhat <- apply(predict(mod, X2)
                ,2, mean)
  d1 <- sqrt(mean((y2-yhat)^2))
  expect_that(d1, is_less_than(0.6))
})


xx <- c(0.5, .5)
y <- 6*ff1(xx) + (rgamma(3000, 1.5, 1) - 1.5)
yhat <- predict(mod, matrix(xx, nrow=1))
