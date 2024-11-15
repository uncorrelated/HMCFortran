library(hmcfortran)

gen_data <- function(n, beta = c(0.5, -1, 0.5, 0), mu = c(-0.02,-0.054,0.071,0.025), 
    S = matrix(c(15.206,7.327,26.692,13.811,7.327,14.883,24.389,14.606,26.692,24.389,162.93,37.233,13.811,14.606,37.233,31.147), 4, 4)){
    library(mvtnorm)
    X <- rmvnorm(n, mu, S)
	colnames(X) <- sprintf("x%d", 1:ncol(X))
    y <- X %*% beta
    m <- 3
    ng <- 3 # 応答の種類
    g <- seq(2, 3, length.out = ng - 1)
    a <- matrix(g, n, length(g), byrow = TRUE) - replicate(length(g), c(y))
    p <- 1/(1 + 1/exp(a))
    q <- runif(n)
    r <- rep(1, n)
    for(j in 2:ng){
        r[q >= p[, j - 1]] <- j
    }
    no_missing <- list(ans = ordered(r, levels = 1:ng, labels = letters[1:ng]), X = X, g = g, b = beta)
    missing <- (r == 1) & 0.75 > runif(length(r)) 
    missing_r <- (1:nrow(X))[missing]
    missing_c <- round(runif(length(missing_r), 0.5, 0.5 + length(mu) - 1))
    for(i in 1:length(missing_r)){
        X[missing_r[i], missing_c[i]] <- NA
    }
    missing <- list(ans = ordered(r, levels = 1:ng, labels = letters[1:ng]), X = X, g = g, b = beta)
    list(no_missing = no_missing, missing = missing)
}

set.seed(1234)

data <- gen_data(100)

df01 <- with(data$missing, data.frame(ans, X))

beta.mu <- rep(0, ncol(df01) - 2)
beta.Sigma <- diag(1000, ncol(df01) - 2)

nok_ans <- length(levels(df01$ans))
gamma.mu <- seq(1, 2, length.out = nok_ans - 1)
gamma.Sigma <- diag(1000, nok_ans - 1, nok_ans - 1)

X_mu <- rep(0, 4)
X_mu_Sigma <- diag(100, 4)

r.imp  <- hmc.ologit.imp(ans ~ x1 + x2 + x3, 
    beta.mu, beta.Sigma, 
    gamma.mu, gamma.Sigma, 
    X_mu, X_mu_Sigma, 
    frml.aux = ~ x4, data = df01, N = 3000, nchains = 2)

log_BF(r.imp, 1, 0)

prdct.dist <- predict(r.imp)
hist(prdct.dist, breaks = 100)
