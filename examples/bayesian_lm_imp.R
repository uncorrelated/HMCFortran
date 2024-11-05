library(hmcfortran)

gen_data <- function(n, mu = c(2, -3), S = matrix(c(2, -1.5, -1.5, 2), 2, 2)){
    library(mvtnorm)
    X <- rmvnorm(n, mu, S)
    beta <- c(1, -1, 1)
    y <- cbind(1, X) %*% beta + rnorm(n)
    x <- X[,1]
    z <- X[,2]
    p <- 1/(1  + exp(-(4 + y)))
    mx <- x 
    mx[runif(n) < p] <- NA
    data.frame(y, x, mx, z)
}

set.seed(1029)
df01 <- gen_data(50)

np <- 3
r <- hmc.blm.imp(y ~ mx + z, 1, 1, rep(0, np), diag(100, np), rep(0, np - 1), diag(100, np - 1), data = df01, N = 3000)

plot(r, ask = TRUE)

summary(r)

log_ml(r)

log_BF(r, 2, 0)
