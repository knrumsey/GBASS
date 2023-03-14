library(duqling)
library(lhs)
library(BASS)
library(BART)


# 5 dimensional input
stochastic_piston <- function(x){
  x[2] <- x[2] + 1
  Ta <- rbeta(1, 10, 1.2)
  Pa <- runif(1, 0.5 - 0, 0.5 + 0)
  xx <- c(x[1:4], Pa, Ta, x[5])
  piston(xx, scale01=TRUE)
}

#Simulate data
X <- maximinLHS(500, 5)
#X <- augmentLHS(X, m=4000)
X <- rbind(rbind(rbind(X, X), X), X)
y <- apply(X, 1, stochastic_piston)
y <- y/sd(y)
hist(y)

# BASS Model
mod_bass <- bass(X, y)

# GBASS Model
mod_gbass <- nwbass(X, y,
                    a_lambda=4, b_lambda=400/nrow(X),
                    m_beta = 0, s_beta=3,
                    m_gamma=100, s_gamma=35,
                    w_prior=list(type="GIG", p=-5, a=0, b=5, prop_sigma=1000),
                    nmcmc=5000, nburn=4001,
                    vsig=0.1)

