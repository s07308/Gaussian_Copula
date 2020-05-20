library(mvtnorm)

## generate (Y, delta, X), where X is ternary and Y is lifetime with censoring
## sample size
size <- 500

## correlation of latent Gaussian distribution
rho <- 0.5 

## correlation matrix of latent Gaussian
corr.matrix <- matrix(c(1, rho, rho, 1), nrow = 2) 

## Cutoff points to make ternary variable
cut.off <- c(-1, 1) 

## Draw a random sample from the latent Gaussian distribution
Z <- rmvnorm(n = size, mean = c(0, 0), sigma = corr.matrix)

## Assume that h0(t)=exp(t), and make the lifetime variable
T1 <- exp(Z[, 1])

## Do discretization to make ternary variable
X1 <- ifelse(Z[, 2] < cut.off[1], 0,
             ifelse(Z[, 2] < cut.off[2], 1, 2))

## generate the censoring variable
C <- rexp(n = size, rate = 0.5)

## Make the censored observation for lifetime
Y1 <- ifelse(T1 < C, T1, C)

## Indicator to indicate that the lifetime is censored or not
delta <- ifelse(T1 < C, 1, 0)