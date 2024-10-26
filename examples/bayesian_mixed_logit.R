library(hmcfortran)

gendata <- function(noc = 3, nrow = 100, ncol = 3){

    nb <- (ncol + 1)*noc
    beta <- matrix((-1)^(1:nb + 1) * 1:nb, ncol + 1)
    X <- matrix(c(rep(1, nrow), runif(ncol*nrow, -1, 1)), nrow)
    colnames(X) <- sprintf("x%d", 1:(ncol + 1) - 1)

    names <- letters[1:noc]
    cnd <- matrix(runif(nrow*noc, 1, 5), nrow, noc, dimnames = list(NULL, paste("cnd", names[1:noc], sep = ".")))
    alpha <- -3

    y <- alpha * cnd + X %*% beta
    colnames(y) <- sprintf("y%d", 1:noc)

    p_n <- exp(y)
    p_d <- apply(p_n, 1, sum)
    p <- p_n/p_d
    colnames(p) <- sprintf("pr%d", 1:noc) 

    q <- p
    for(j in 2:noc){
        q[,j] <- q[,j-1] + q[,j]
    }
    r <- runif(nrow)
    r <- replicate(noc, r)
    a <- apply(q < r, 1, sum) + 1

    ans <- names[a]

    cbind(ans = as.factor(ans), cnd, as.data.frame(X[, -1]))
}

set.seed(1234)
df01 <- gendata(3, 500, 3)

frml_cnd <- ans ~ cnd.a + cnd.b + cnd.c
frml_mnl <- ans ~ x1 + x2 + x3

mu <- rep(0, 1 + 4*2)
Sigma <- diag(100, length(mu))

r_hmc <- hmc.mlogit(frml_mnl, frml_cnd, mu, Sigma, data = df01, N = 6000)

plot(r_hmc)

gelman.diag(r_hmc)

summary(r_hmc)
