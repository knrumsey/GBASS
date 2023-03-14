library(GBASS)
library(BASS)
library(duqling)
library(lhs)

X <- maximinLHS(1000, 2)
y <- 10*apply(X, 1, dms_simple) + (rgamma(1000, 2, 0.05) - 40)/sqrt(2/0.05^2)

mod <- nwbass(X, y, s_beta=10, m_gamma=100, s_gamma=30, lag_beta=1, scale=1,
               w_prior=list(type="GIG", p=0, a=0, b=0, prop_sigma=0.2), maxInt = 2)
yhat <- apply(predict(mod, X, mcmc.use=seq(1, 500, by=1)), 1, mean)
yhat <- yhat + mean(sqrt(mod$w)*mod$bet/mod$gam)

par(mfrow=c(1,2))
plot(y, yhat)
abline(0, 1, lwd=2, lty=3, col='orange')
hist(y-yhat, freq=F, breaks=30)
curve(dgamma(40 + x*sqrt(2/.05^2), 2, 0.05)*sqrt(2/.05^2), add=T, col='orange', lwd=2)


mod2 <- bass(X, y)
yhat2 <- apply(predict(mod2, X), 2, mean)
par(mfrow=c(1,2))
plot(y, yhat2)
abline(0, 1, lwd=2, lty=3, col='orange')
hist(y-yhat2, freq=F, breaks=30)
curve(dgamma(40 + x*sqrt(2/.05^2), 2, 0.05)*sqrt(2/.05^2), from=-3, to=10, add=TRUE, col='orange', lwd=2)




#Pick a reference point
xx <- c(0.5, 0.5)


# Plot the true response distribution
curve(28.28*dgamma(28.28*x - 28.28*10*dms_simple(xx) + 40, 2, 0.04),
      from=30, to=45,
      xlab="output", ylab="density",
      lwd=2)


# Plot the predicted response distribution based on BASS
ffbass <- function(xx, mod, x0){
  res <- length(xx)
  tmp <- predict(mod, matrix(x0, nrow=1))
  for(i in seq_along(xx)){
    res[i] <- mean(dnorm(xx[i], tmp, sqrt(mod$s2)))
  }
  return(res)
}
curve(ffbass(x, mod2, x0=xx), add=TRUE, col='dodgerblue', lwd=2)

# Plot the predicted response distribution based on NWBASS
tmp <- predict(mod, matrix(xx, nrow=1))
MCMC <- 1000
K    <- 100  #get K samples for each of MCMC iteration
ddd <- rep(NA, K*MCMC)
for(m in 1:MCMC){
  xi <- rnorm(K)

}






