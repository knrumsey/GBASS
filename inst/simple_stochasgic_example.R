library(GBASS)
library(BASS)
library(duqling)
library(lhs)

X <- maximinLHS(1000, 2)
y <- 10*apply(X, 1, dms_simple) + (rgamma(1000, 2, 0.05) - 40)/sqrt(2/0.05^2)

mod <- nwbass(X, y, s_beta=10, m_gamma=100, s_gamma=30, lag_beta=1, scale=1,
               w_prior=list(type="GIG", p=0, a=0, b=0, prop_sigma=0.2), maxInt = 2,
              nmcmc = 50000, thin=5)
yhat <- apply(predict(mod, X), 2, mean)
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
curve(28.28*dgamma(28.28*x - 28.28*10*dms_simple(xx) + 40, 2, 0.05),
      from=30, to=45,
      xlab="output", ylab="density",
      lwd=2,
      ylim=c(0, 0.8))


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
K    <- 100  # get K samples for each MCMC iteration
ddd <- rep(NA, K*MCMC)

for(m in 1:MCMC){
  xi <- rnorm(K)
  vv <- replicate(K, rgig2(-1/2, mod$gam[m]^2, 1))
  ddd[((m-1)*K + 1):(m*K)] <- tmp[m] +
    sqrt(mod$w[m]) * (mod$bet[m] * vv + sqrt(vv) * xi)
}

lines(density(ddd, bw=0.4), col='firebrick', lwd=2)
legend("topright",
       legend=c("Truth", "BASS", "NWBASS"),
       col=c("black", "dodgerblue", "firebrick"),
       lwd=2,
       bty="n")




preds1 <- predict(mod, X, predictive = TRUE)
preds2 <- predict(mod, X, predictive = FALSE, bias_correct = FALSE)
preds3 <- predict(mod, X, predictive = FALSE, bias_correct = TRUE)

# This looks perfect.
plot(y, apply(preds1, 2, mean))
abline(0, 1, col='red')

# This looks bad (as expected)
plot(y, apply(preds2, 2, mean))
abline(0, 1, col='red')

# This seems to "over-adjust" worse than expected
# I don't understand why it doesn't match the preds1 case better
plot(y, apply(preds3, 2, mean))
abline(0, 1, col='red')


