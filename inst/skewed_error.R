# Semi-random seed
if(!exists("foobar")) foobar <- 101
set.seed(foobar)
foobar <- foobar+1
# Seed done

X <- lhs::maximinLHS(200, 2)
X <- rbind(X, X, X, X, X)
n <- nrow(X)
y <- apply(X, 1, ff4)*6 - (rgamma(n, 1.5, 1) - 1.5)

#mod2 <- nwbass2(X, y, maxInt=2, m_gamma = 1, s_gamma = 0.1, scale=1000, m_beta=0, s_beta=1)
mod1 <- bass(X, y)
tic()
mod2 <- nwbass2(X, y, maxInt=2, m_gamma=90, s_gamma=25, scale=1, m_beta=0, s_beta=1, lag_beta=100)
toc()
plot(mod2)


nn <- 1000
ww <- mean(mod2$w)
cc <- 1
beta <- mean(mod2$beta)
xi <- rnorm(nn, 0, sqrt(cc))
gamma <- mean(mod2$gamma)
vv <- rep(NA, nn)
for(i in 1:nn){
  vv[i] <- rgig2(-1/2, gamma^2, 1)
}
epsilon <- sqrt(ww)*(beta*vv + sqrt(vv)*xi)

par(mfrow=c(1,1))
plot(density(-(rgamma(1000000, 1.5, 1) - 1.5)), lwd=3, xlim=c(-3, 6), ylim=c(0, 0.52))
lines(density(epsilon - mean(epsilon)), col='orange', lwd=3, lty=2)
curve(dnorm(x, 0, mean(sqrt(mod1$s2))), add=TRUE, col='dodgerblue', lwd=3, lty=3)





#
# w_prior <- build_GIG(0, 0, 0)
# w_prior$lb <- 0
# mod <- gbass(X, y, maxInt=2, w_prior=w_prior)
# plot(mod)
#
#
# w_prior <- build_GIG(0, 0, 0)
# w_prior$lb <- 0
# w_prior$prop_sigma <- 0.1
# qmod <- qbass(X, y, q=0.75, w_prior=w_prior, a_lambda=1, b_lambda=0.01)
#
#
# Xt <- lhs::maximinLHS(500, 2)
# yt <- apply(Xt, 1, ff4)*6
# yhat <- apply(predict(qmod, Xt), 2, mean)
# plot(yt, yhat)
# abline(0, 1, col='red')
#
# plot(yt + qgamma(0.75, 1.5, 1) - 1.5, yhat)
# abline(0, 1, col='red')
