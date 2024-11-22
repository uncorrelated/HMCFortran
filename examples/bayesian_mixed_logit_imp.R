gen_data <- function(noc = 3, nrow = 100, ncol = 3, mu = c(0, 0, 0), Sigma = matrix(c(0.147,0.135,0.143,0.135,0.299,-0.027,0.143,-0.027,0.316), 3)){
    library(mvtnorm)

    # 線形部分の値をつくる
    nb <- (ncol + 1)*noc
    beta <- matrix((-1)^(1:nb + 1) * 1:nb, ncol + 1)
    beta[ncol + 1, ] <- 0
    X <- cbind(1, rmvnorm(nrow, mu, Sigma))
    colnames(X) <- sprintf("x%d", 1:(ncol + 1) - 1)

    # 線形部分の係数の確認
    rownames(beta) <- colnames(X)
    # print(beta[, 2:3] - beta[,1])

	# conditionalな変数の値をつくる
    names <- letters[1:noc]
    cnd <- matrix(runif(nrow*noc, 1, 5), nrow, noc, dimnames = list(NULL, paste("cnd", names[1:noc], sep = ".")))
    alpha <- -3 # conditionalな変数の係数

    y <- alpha * cnd + X %*% beta
    colnames(y) <- sprintf("y%d", 1:noc)

    # 選択肢ごとの確率をつくる
    p_n <- exp(y)
    p_d <- apply(p_n, 1, sum)
    p <- p_n/p_d
    colnames(p) <- sprintf("pr%d", 1:noc)

    # 乱数から選択肢を選ぶ
    q <- p
    for(j in 2:noc){
        q[,j] <- q[,j-1] + q[,j]
    }
    r <- runif(nrow)
    r <- replicate(noc, r)
    a <- apply(q < r, 1, sum) + 1

    ans <- names[a]

    X_missing <- X
    # x1にだけ欠損値が入る
    # for(j in 2:ncol(X)){
    for(j in 2:2){
        # missing completely at random
        X_missing[runif(nrow(X)) < 0.2, j] <- NA
    }
    # データフレームに選択肢と説明変数をまとめる
    list(no_missing = cbind(ans = as.factor(ans), cnd, as.data.frame(X[, -1])), 
        missing = cbind(ans = as.factor(ans), cnd, as.data.frame(X_missing[, -1])))
}

set.seed(1234)
data <- gen_data(nrow = 150)

library(hmcfortran)

frml_mnl <- ans ~ x1 + x2
frml_cnd <- ~ cnd.a + cnd.b + cnd.c
frml_aux <- ~ x3

noc <- 1
nok <- max(as.integer(data$missing$ans))
beta.mu <- rep(1, noc + (nok - 1)*3)
beta.Sigma <- diag(1000, length(beta.mu))

X_mu <- rep(0, 6)
X_mu_Sigma <- diag(100, length(X_mu))

r.fbim <- hmc.mlogit.imp(frml_mnl, frml_cnd, beta.mu, beta.Sigma, X_mu, X_mu_Sigma, frml.aux = frml_aux, data = data$missing, nchains = 2, N = 1500)
