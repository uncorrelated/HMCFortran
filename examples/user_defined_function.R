library(hmcfortran)

set.seed(2024)
# data generator for practice
gendata <- function(nr, nc, sd = 1){
	X <- matrix(c(rep(1, nr), runif(nr*(nc-1), -2, 2)), nr)
	colnames(X) <- c("Const.", sprintf("x%d", 1:(nc - 1)))
	beta <- c(0.5, rep(1, nc - 1))*(-1)^(1:nc)
	y <- X %*% beta + rnorm(nr, sd = sd)
	cbind(y, as.data.frame(X))
}
data <- gendata(100, 3)

frml <- y ~ x1 + x2
data <- model.frame(frml, data)
X <- model.matrix(frml, data)
y <- model.response(data)

objf <- function(p, X, y){
    nc <- ncol(X)
    sum(dnorm(y - X %*% p[2:(1 + nc)], sd = abs(p[1]), log = TRUE))
}

init.p <- c(1, runif(ncol(X), -1, 1))
r_hmc <- hmc.usrfunc(init.p, objf, X = X, y = y, N = 3000, BI = 300, nchains = 2)

plot(r_hmc)

gelman.diag(r_hmc)

summary(r_hmc)

log_ml(r_hmc)

objfg <- function(p, X, y){
    theta <- p[1]
    beta <- p[-1]
    res <- y - X %*% beta
    dtheta <- -1*nrow(X)/theta + sum(res^2)/theta^3
    dbeta <- apply(X*rep(res, ncol(X)), 2, sum)/theta^2
    c(dtheta, dbeta)
}

r_hmc_g <- hmc.usrfunc(init.p, objf, objfg, X = X, y = y, N = 5000, BI = 1000, nchains = 2)

sprintf("elapsed time with numerical gradient: %.3f", r_hmc$system.time[3])
sprintf("elapsed time with analytical gradient: %.3f", r_hmc_g$system.time[3])
