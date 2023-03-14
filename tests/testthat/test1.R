test_that("GBASS with t-distribution", {
  cat('simple test function')

  X <- lhs::maximinLHS(300, 2)
  y <- apply(X, 1, ff1)
  mod <- gbass(X, y, nmcmc=2000, nburn=1001)



  d1 <- sum(abs(Cbass - Ctrue))
  expect_that(d1, is_less_than(0.04))
})
