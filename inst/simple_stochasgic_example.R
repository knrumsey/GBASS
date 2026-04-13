library(GBASS)
library(BASS)
library(duqling)
library(lhs)

X <- maximinLHS(1000, 2)
y <- 10*apply(X, 1, dms_simple) + (rgamma(1000, 2, 0.05) - 40)/sqrt(2/0.05^2)

# Fit nwbass model
mod1 <- nwbass(X, y, s_beta=10, m_gamma=100, s_gamma=30, lag_beta=1, scale=1,
               w_prior=list(type="GIG", p=0, a=0, b=0, prop_sigma=0.2), maxInt = 2)

# Fit regular bass model
mod2 <- bass(X, y)

#Pick a reference point
xx <- c(0.5, 0.5)

# Plot the true response distribution
curve(28.28*dgamma(28.28*x - 28.28*10*dms_simple(xx) + 40, 2, 0.05),
      from=32, to=43,
      xlab="output", ylab="density",
      lwd=2,
      ylim=c(0, 0.65))

# Plot the predicted response distribution based on NWBASS
preds1 <- predict(mod1, matrix(xx, nrow=1))
d1 <- density(preds1)
polygon(d1, border="orange", lwd=2)

# Plot the predicted response distribution based on BASS
preds2 <- predict(mod2, matrix(xx, nrow=1), nugget=TRUE)
d2 <- density(preds2)
polygon(d2, border="dodgerblue", lwd=2)

# Add leged
legend("topright",
       legend=c("Truth", "BASS", "NWBASS"),
       col=c("black", "dodgerblue", "firebrick"),
       lwd=2,
       bty="n")

