library(hmcfortran)

gendata <- function(n = 100, link = "log"){
    x <- runif(n, -1, 1)
    z <- runif(n, -1, 1)
    Xb <- 2*x - z
    if("identity" == link){
        lambda <- Xb + min(0, Xb)
    } else if("sqrt" == link){
        lambda <- sqrt(Xb + min(0, Xb))
    } else {
        lambda <- exp(Xb)
    }
    y <- rpois(n, lambda)
    data.frame(y, x, z)
}

set.seed(123)
df01 <- gendata(100)

r_hmc <- hmc.poisson_exp(y ~ x + z, rep(0, 3), diag(100, 3), data = df01)

plot(r_hmc)

gelman.diag(r_hmc)

summary(r_hmc)

# calculate BF_{01}
exp(log_BF(r_hmc, 2, 0))
