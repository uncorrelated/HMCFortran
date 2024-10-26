library(hmcfortran)

set.seed(1234)

# data generator for practice
gendata <- function(n = 100){
	b <- c(1, 2, -3)
	X <- cbind(1, matrix(runif(n*(length(b)-1)), n))
	colnames(X) <- c("Const.", "x", "z")
	p <- exp(X%*%b)/(1+exp(X%*%b))
	y <- 0 + 1*(p > runif(length(p)))
	cbind(y, as.data.frame(X[,-1]))
}

df01 <- gendata(1000)

r_glm <- glm(y ~ x + z, data = df01, family = binomial())

# prior/hyper parameters
beta.mu <- c(0, 0, 0)
beta.sigma <- diag(100, 3, 3)

# do Hamiltonian Monte Carlo method to estimate a logistic model
r_hmc_logit <- hmc.blogit(y ~ x + z, beta.mu, beta.sigma, data = df01, N = 3000, BI = 1000)

plot(r_hmc_logit)

gelman.diag(r_hmc_logit)

summary(r_hmc_logit)

predicted_distribution <- predict(r_hmc_logit, apply(r_hmc_logit$input$X, 2, mean))

hist(predicted_distribution, breaks = 20)

# get the Laplace Approximated Marginal Log-Likelihood
# log_ml(r_hmc_logit)

# Bayes Factors of the x = 0 constrain model
log_BF(r_hmc_logit, 1, 0)

