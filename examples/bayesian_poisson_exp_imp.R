library(hmcfortran)

mkdata <- function(n = 100, link = "log"){
    w <- runif(n, -1, 1)
    x <- runif(n, -1, 1)/2 + w/2
    z <- runif(n, -1, 1)/3 + 2*w/3
    Xb <- 2*x - z
    if("identity" == link){
        lambda <- Xb + min(0, Xb)
    } else if("sqrt" == link){
        lambda <- sqrt(Xb + min(0, Xb))
    } else {
        lambda <- exp(Xb)
    }
    y <- rpois(n, lambda)

    no_missing <- data.frame(y, x, z, w)

    x[(2 - y)/3 > runif(length(y))] <- NA
    z[(3 - y)/4 > runif(length(y))] <- NA

    list(missing = data.frame(y, x, z, w), no_missing = no_missing)
}

# set.seed(123)
data <- mkdata(100)

r_hmc <- hmc.poisson_exp.imp(y ~ x + z, rep(0, 3), diag(100, 3), rep(0, 3), diag(100, 3), frml.aux = ~ w, data = data$missing, N = 1500, nchains = 2)
