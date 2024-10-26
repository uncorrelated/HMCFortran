library(hmcfortran)

# data generator for practice
gendata <- function(nr, nc, sd = 1){
	X <- matrix(c(rep(1, nr), runif(nr*(nc-1), -2, 2)), nr)
	colnames(X) <- c("Const.", sprintf("x%d", 1:(nc - 1)))
	beta <- c(0.5, rep(1, nc - 1))*(-1)^(1:nc)
	y <- X %*% beta + rnorm(nr, sd = sd)
	cbind(y, as.data.frame(X))
}

set.seed(2024)

np <- 4 # number of parameters
data <- gendata(150, np - 1)

# prior/hyper parameters
ig.alpha <- 1
ig.gamma <- 1
beta.mean <- rep(0, np - 1)
beta.Sigma <- diag(1e+10, np - 1, np - 1)

# do Hamiltonian Monte Carlo method to estimate a liner model
r_hmc_lm <- hmc.blm(y ~ x1 + x2, ig.alpha, ig.gamma, beta.mean, beta.Sigma, data = data, N = 3000, nchains = 2)

# get the marginal likelihood
mlf_1 <- log_ml(r_hmc_lm)

# alternative model: third parameter, a, is zero.
r_hmc_lm_0 <- hmc.blm(y ~ x2, 1, 1, rep(0, np - 2), diag(1e+10, np - 2, np - 2), data = data, N = 3000)

# get the marginal likelihood of alternative model
mlf_0 <- log_ml(r_hmc_lm_0)

cat("Laplace Approximated BF_{10}:", exp(mlf_1 - mlf_0), "\n")

BF <- exp(log_BF(r_hmc_lm, 3, 0))
cat("IWAMDE Savage-Dickey ratio BF_{01}:", BF, " BF_{10}:", 1/BF, "\n")

# show the prediction, where explanatory variables are the mean of X.
prdct <- predict(r_hmc_lm, X = apply(r_hmc_lm$input$X, 2, mean))
hist(prdct, breaks = 20)

# call coda::plot
plot(r_hmc_lm)

# calcurate R-hat
gelman.diag(r_hmc_lm)

# call other functions of coda package
# gelman.plot(r_hmc_lm$mcmc.list)
# geweke.plot(r_hmc_lm$mcmc.list)
# crosscorr.plot(r_hmc_lm$mcmc.list)
# autocorr.diag(r_hmc_lm$mcmc.list)
# autocorr.plot(r_hmc_lm$mcmc.list)

# call coda::summary
summary(r_hmc_lm)

# call coda::HPDinterval
HPDinterval(r_hmc_lm)

# statistics of the sample of MCMC
mean(r_hmc_lm)
median(r_hmc_lm)
quantile(r_hmc_lm, 0.5)
