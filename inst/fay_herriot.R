library(duqling)
show_function_2d("banana")
show_function_2d("twin_galaxies")

n <- 50
sigma_v <- 0
X <- lhs::randomLHS(n, 2)
points(X)


fm <- function (x){
  x[1]^2 + x[1] * x[2] +  x[2]^3 / 9
}

fv <- function (x){
  1.3356 * (1.5 * (1 - x[1]) + exp(2 * x[1] - 1) * sin(3 * pi *
    (x[1] - 0.6)^2) + exp(3 * (x[2] - 0.5)) * sin(4 * pi * (x[2] -
     0.9)^2))
}

# GENERATE DATA
n <- 50
sigma_v <- 0.1
X <- lhs::randomLHS(n, 2)
mu <- apply(X, 1, fm)
mu <- mu / sd(mu)
N_pop <- 1 + round(apply(X, 1, fv))
d <- 1/N_pop
y <- rnorm(n, mu, sqrt(d)) + rnorm(n, 0, sigma_v)

# Fit model
fit <- fhbass(X, y, d)
