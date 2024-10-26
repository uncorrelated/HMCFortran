library(hmcfortran)

set.seed(1234)

ds <- with(new.env(), {
    set.seed(111)
    ng <- 3 # 応答の種類
    g <- seq(3, 4, length.out = ng - 1)
    n <- 300 # サンプルサイズ
    m <- 3 # 説明変数の数
    X <- matrix(runif(m*n), n, m)
#    b <- matrix(runif(m, min = 1, max = 2), m, 1)
	b <- matrix((-1)^(1:m + 1)*(1:m), m, 1)
    y <- X %*% b
    a <- matrix(g, n, length(g), byrow = TRUE) - replicate(length(g), c(y))
    p <- 1/(1 + 1/exp(a))
    q <- runif(n)
    r <- rep(1, n)
    for(j in 2:ng){
        r[q >= p[, j - 1]] <- j
    }
    list(ans = ordered(r, levels = 1:ng, labels = letters[1:ng]), X = X, g = g, b = b)
})

beta.mu <- rep(0, ncol(ds$X))
beta.sigma <- diag(1000, ncol(ds$X))

nok_ans <- length(levels(ds$ans))
gamma.mu <- seq(1, 2, length.out = nok_ans - 1)
gamma.sigma <- diag(1000, nok_ans - 1, nok_ans - 1)

r_hmc <- hmc.ologit(ans ~ X, beta.mu, beta.Sigma, gamma.mu, gamma.Sigma, data = ds)

plot(r_hmc)
summary(r_hmc)
