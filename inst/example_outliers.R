library(GBASS)
library(BASS)
library(lhs)

# Define function
ishigami <- function (x){
  ab = c(7, 0.1)
  RR <- cbind(rep(-base::pi, 3), rep(base::pi, 3))
  x[1:3] <- x[1:3] * (RR[, 2] - RR[, 1]) + RR[, 1]
  res <- sin(x[1]) + ab[1] * sin(x[2])^2 + ab[2] * x[3]^4 * sin(x[1])
  return(res)
}

# Get data
X <- lhs::randomLHS(500, 3)
y <- apply(X, 1, ishigami)
y <- y + rnorm(500, 0, sd(y) / 100)

# Add outliers to first 100 data points
y[1:100] <- y[1:100] + rt(100, 5)

# Fit regular bass model
mod1 <- bass(X, y)

# Fit t-bass model
mod2 <- tbass(X, y, df=5)

# Fit horseshoe model
mod3 <- hbass(X, y)

# Get predictions
pred1 <- predict(mod1, X, nugget=TRUE)
pred2 <- predict(mod2, X, predictive=TRUE)
pred3 <- predict(mod3, X, predictive=TRUE)

# Make some plots
par(mfrow=c(1,3))
colors <- rep("black", 500)
colors[1:100] <- 'red'

plot(colMeans(pred1), y, col=colors)
abline(0, 1, col='dodgerblue', lwd=2)
plot(colMeans(pred2), y, col=colors)
abline(0, 1, col='dodgerblue', lwd=2)
plot(colMeans(pred3), y, col=colors)
abline(0, 1, col='dodgerblue', lwd=2)
par(mfrow=c(1,1))

# Test on test set
X_test <- lhs::randomLHS(1000, 3)
f_test <- apply(X_test, 1, ishigami)

# Small observational noise, but no outliers
y_test <- f_test + rnorm(1000, 0, sd(y) / 100)

# Posterior mean predictions
pred1_test <- predict(mod1, X_test, nugget = FALSE)
pred2_test <- predict(mod2, X_test, predictive = FALSE)
pred3_test <- predict(mod3, X_test, predictive = FALSE)

yhat1 <- colMeans(pred1_test)
yhat2 <- colMeans(pred2_test)
yhat3 <- colMeans(pred3_test)

# Compare against the true signal
rmse <- function(truth, fit) sqrt(mean((truth - fit)^2))
mae  <- function(truth, fit) mean(abs(truth - fit))

results <- data.frame(
  model = c("BASS", "tBASS", "hBASS"),
  RMSE = c(
    rmse(f_test, yhat1),
    rmse(f_test, yhat2),
    rmse(f_test, yhat3)
  ),
  MAE = c(
    mae(f_test, yhat1),
    mae(f_test, yhat2),
    mae(f_test, yhat3)
  )
)

print(results)
cat("\nBest by RMSE:", results$model[which.min(results$RMSE)], "\n")

# Simple plot
par(mfrow = c(1, 3))
plot(f_test, yhat1, pch = 16, cex = 0.5,
     main = paste0("BASS\nRMSE = ", round(results$RMSE[1], 3)),
     xlab = "True signal", ylab = "Predicted mean")
abline(0, 1, col = "dodgerblue", lwd = 2)

plot(f_test, yhat2, pch = 16, cex = 0.5,
     main = paste0("tBASS\nRMSE = ", round(results$RMSE[2], 3)),
     xlab = "True signal", ylab = "Predicted mean")
abline(0, 1, col = "dodgerblue", lwd = 2)

plot(f_test, yhat3, pch = 16, cex = 0.5,
     main = paste0("hBASS\nRMSE = ", round(results$RMSE[3], 3)),
     xlab = "True signal", ylab = "Predicted mean")
abline(0, 1, col = "dodgerblue", lwd = 2)







