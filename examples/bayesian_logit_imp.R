library(hmcfortran)

gen_data <- function(n, beta = c(1, 1, -1, 0), mu = c(0.936,-0.089,0.51), S = matrix(c(22.481,-8.381,14.516,-8.381,5.315,-5.607,14.516,-5.607,9.769), 3)){
    library(mvtnorm)
    X <- cbind(1, rmvnorm(n, mu, S))
	colnames(X) <- c("Const.", sprintf("x%d", 1:(ncol(X) - 1)))
    y <- X %*% beta
	p <- exp(y)/(1 + exp(y))
	ans <- 0 + 1*(p > runif(length(p)))
    Xm <- X
    missing <- (ans == 0) & 0.75 > runif(length(ans)) 
    missing_r <- (1:nrow(X))[missing]
    missing_c <- 1 + round(runif(length(missing_r), 0.5, 0.5 + length(mu) - 1))
    for(i in 1:length(missing_r)){
        Xm[missing_r[i], missing_c[i]] <- NA
    }
	df01 <- cbind(data.frame(ans), X[, -1])
	df02 <- cbind(data.frame(ans), Xm[, -1])
    list(no_missing = df01, missing = df02)
}

set.seed(1110)
data <- gen_data(100)

r_lm <- lm(ans ~ x1 + x2, data$missing)
beta.mu <- coef(r_lm)
beta.sigma <- diag(100, 3, 3)

r.imp  <- hmc.blogit.imp(ans ~ x1 + x2, beta.mu, beta.sigma, rep(0, 3), diag(100, 3), frml.aux = ~ x3, data = data$missing, N = 3000, nchains = 2)
