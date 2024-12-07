### 線形回帰で使う関数群 ###

#
# frml: モデル式
# data: データフレーム
# ig.alpha: 分散の事前分布になる逆ガンマ関数のハイパーパラメーター
# ig.beta: 〃
# beta.mu: 線形予測部分のパラメーターの事前分布になる多変量正規分布の平均を表すハイパーパラメーター
# beta.Sigma: 〃 の分散共分散行列を 〃
# N: サンプリングする数
# BI: Burin-in（捨てる初期部分）
# adjustSpan: epsilonsを調整する期間/-1を入れると初期値の調整のみ行う
# nchains: 生成する乱数列の本数
# L, epsilons, seed: Hamiltonian Monte Carlo法のパラメーター
#
hmc.blm <- function(frml, ig.alpha, ig.beta, beta.mu, beta.Sigma, data = NULL, 
	N = 3000, BI = as.integer(N*0.1), adjustSpan = 0, L = 20, epsilons = NULL, nchains = 2, seed = NULL){

	seed <- set_seeds(seed, nchains)

	data <- model.frame(frml, data)
	X <- model.matrix(frml, data)
	y <- model.response(data)

	init.p <- c(0.5, runif(ncol(X)))
	np <- length(init.p)

	if(length(beta.mu) != np - 1) stop("The length of beta.mu is wrong.")
	if(!is.matrix(beta.Sigma)) stop("beta.Sigma must be a matrix.")
	if(np - 1 != ncol(beta.Sigma) || np - 1 != nrow(beta.Sigma)) stop("The length of beta.Sigma is wrong.")

	hp <- c(ig.alpha, ig.alpha, beta.mu, chol2inv(chol(beta.Sigma)))
	nhp <- length(hp)

	r_optim <- .Fortran("optim_bayesian_lm",
		nrow(X), ncol(X), as.double(X), as.double(y),
		np, as.double(init.p),
		nhp, as.double(hp),
		as.double(0), as.double(matrix(0, np, np)), integer(1))
	r_optim <- list(p = r_optim[[6]], f = r_optim[[9]], h = matrix(r_optim[[10]], np), info = r_optim[[11]])

	init.p <- r_optim$p
	epsilons <- adjust_epsilons(epsilons, L, r_optim$h)

	st <- system.time({
		# 並行処理の準備
		library(parallel)
		library(doParallel)
		library(foreach)
		cl <- makeCluster(max(1, detectCores()))
		registerDoParallel(cl)

		r_sub <- foreach(i = 1:nchains) %dopar% {

			r <- .Fortran("hmc_fit_bayesian_lm",
				nrow(X), ncol(X), X, y,
				as.integer(N + BI), matrix(as.double(0), N + BI, np),
				np, as.double(init.p),
				nhp, as.double(hp),
				epsilons, as.integer(adjustSpan),
				as.integer(L),
				FALSE,
				as.double(diag(np)),
				FALSE,
				as.integer(seed[i]),
				as.double(0))

			list(s = r[[6]], ar = r[[18]])
		}

		# 並行処理の終了
		stopCluster(cl)
	})

	r <- makeResultObj(
		"blm", 
		r_sub, 
		BI, 
		colnames = c("(SD)", colnames(X)), 
		alpha = 0.05,
		r_optim,
		ar = mean(sapply(r_sub, "[[", 2)), 
		st = st)

	# 推定に使った変数など
	r$input <- list(
		frml = frml,
		X = X,
		y = y,
		hp = hp,
		init.p = init.p
	)

	# MCMCで使ったパラメーター
	r$mcmc.param <- list(
		BI = BI,
		L = L,
		epsilons = epsilons,
		nchains = nchains,
		seed = seed)

	# 事後分布の要約統計量
	r$posterior <- list(
		beta = list(mu = mean(r)[-1], Sigma = vcov(r)[-1, -1]),
		igamma = list(mean = mean(r)[[1]], var = max(1, vcov(r)[1, 1]))
	)
	r$posterior$igamma$alpha <- with(r$posterior$igamma, mean^2 / var + 2)
	r$posterior$igamma$beta <- with(r$posterior$igamma, mean*(alpha - 1))

	r
}

# 周辺尤度を求める
# log_ml.blm <- function(r){
# 	with(r$input, .Fortran("log_ml_bayesian_lm",
# 		nrow(X), ncol(X), as.double(X), as.double(y),
# 		length(init.p), as.double(init.p),
# 		length(hp), as.double(hp),
# 		as.double(0))[[9]])
# }

# 制約モデルのベイズファクター
log_BF.blm <- function(r, index, value){
	sample <- as.matrix(r$mcmc.list)
	with(r$input, .Fortran("bayes_factor_lm",
		nrow(X), ncol(X), as.double(X), as.double(y),
		nrow(sample), ncol(sample), sample,
		as.integer(length(hp)), as.double(hp),
		as.integer(length(index)), as.integer(index), as.double(value),
		as.double(0))[[13]])
}

# 予測値
predict.blm <- function(r, X = NULL, P = NULL){
	if(is.null(X)) X <- r$input$X
	if(is.null(P)) P <- as.matrix(r$mcmc.list)
	if(!is.matrix(X)) X <- matrix(X, 1, length(X))
	y <- .Fortran("hmc_predict_lm",
		nrow(X), ncol(X), as.double(X),
		nrow(P), ncol(P), as.double(P),
		double(nrow(X) * nrow(P)))[[7]]
}

### 二項ロジットモデル ###
hmc.blogit <- function(frml, beta.mu, beta.Sigma, data = NULL, 
	N = 3000, BI = as.integer(N*0.1), adjustSpan = 0, L = 20, epsilons = NULL, nchains = 2, seed = NULL){

	seed = set_seeds(seed, nchains)

	data <- model.frame(frml, data)
	X <- model.matrix(frml, data)
	y <- model.response(data)

	init.p <- runif(ncol(X))
	np <- length(init.p)

	hp <- c(beta.mu, chol2inv(chol(beta.Sigma)))
	nhp <- length(hp)

	r_optim <- .Fortran("optim_bayesian_logit",
		nrow(X), ncol(X), as.double(X), as.double(y),
		np, as.double(init.p),
		nhp, as.double(hp),
		as.double(0), as.double(matrix(0, np, np)), integer(1))
	r_optim <- list(p = r_optim[[6]], f = r_optim[[9]], h = matrix(r_optim[[10]], np), info = r_optim[[11]])

	init.p <- r_optim$p
	epsilons <- adjust_epsilons(epsilons, L, r_optim$h)

	st <- system.time({
		# 並行処理の準備
		library(parallel)
		library(doParallel)
		library(foreach)
		cl <- makeCluster(max(1, detectCores()))
		registerDoParallel(cl)

		r_sub <- foreach(i = 1:nchains) %dopar% {

			r <- .Fortran("hmc_fit_bayesian_logit",
				nrow(X), ncol(X), X, y,
				as.integer(N + BI), matrix(as.double(0), N + BI, np),
				length(init.p), as.double(init.p),
				as.integer(nhp), as.double(hp),
				epsilons, as.integer(adjustSpan),
				as.integer(L),
				FALSE,
				as.double(diag(np)),
				FALSE,
				as.integer(seed[i]),
				as.double(0))

			list(s = r[[6]], ar = r[[18]])
		}

		# 並行処理の終了
		stopCluster(cl)
	})

	# 予約統計量などをつくっておく
	r <- makeResultObj(
		"blogit", 
		r_sub, 
		BI, 
		colnames = colnames(X), 
		alpha = 0.05,
		r_optim,
		ar = mean(sapply(r_sub, "[[", 2)), 
		st = st)

	# 推定に使った変数
	r$input <- list(
		frml = frml,
		X = X,
		y = y,
		hp = hp,
		init.p = init.p
	)

	# MCMCで使ったパラメーター
	r$mcmc.param <- list(
		BI = BI,
		L = L,
		epsilons = epsilons,
		nchains = nchains,
		seed = seed)

	# 事後分布の要約統計量
	r$posterior <- list(beta = list(mu = mean(r), Sigma = vcov(r)))

	r
}

# 周辺尤度を求める
# log_ml.blogit <- function(r){
# 	with(r$input, .Fortran("log_ml_bayesian_logit",
# 		nrow(X), ncol(X), as.double(X), as.double(y),
# 		length(init.p), as.double(init.p),
# 		length(hp), as.double(hp),
# 		as.double(0))[[9]])
# }

# 制約モデルのベイズファクター
log_BF.blogit <- function(r, index, value){
	sample <- as.matrix(r$mcmc.list)
	with(r$input, .Fortran("bayes_factor_logit",
		nrow(X), ncol(X), as.double(X), as.double(y),
		nrow(sample), ncol(sample), sample,
		as.integer(length(hp)), as.double(hp),
		as.integer(length(index)), as.integer(index), as.double(value),
		as.double(0))[[13]])
}

predict.blogit <- function(r, X = NULL, P = NULL){
	if(is.null(X)) X <- r$input$X
	if(is.null(P)) P <- as.matrix(r$mcmc.list)
	if(!is.matrix(X)) X <- matrix(X, 1, length(X))
	.Fortran("hmc_predict_logit",
		nrow(X), ncol(X), as.double(X),
		nrow(P), ncol(P), as.double(P),
		double(nrow(X) * nrow(P)))[[7]]
}

### 順序ロジットで使う関数群 ###

hmc.ologit <- function(frml, beta.mu, beta.Sigma, gamma.mu, gamma.Sigma, data = NULL, 
	N = 3000, BI = as.integer(N*0.1), adjustSpan = 0, L = 20, epsilons = NULL, nchains = 2, seed = NULL){

	seed = set_seeds(seed, nchains)

	frml <- update(frml, ~ . -1)
	data <- model.frame(frml, data)
	ans <- model.response(data)
	nok_ans <- length(levels(ans))
	X <- model.matrix(frml, data)

	init.p <- as.double(c(runif(ncol(X), -1, 1), seq(1, 2, length.out = nok_ans - 1)))
	np <- length(init.p)

	hp <- c(beta.mu, chol2inv(chol(beta.sigma)), gamma.mu, chol2inv(chol(gamma.sigma)))
	nhp <- length(hp)

	r_optim <- .Fortran("optim_bayesian_ologit",
		nrow(X), ncol(X), as.double(X), as.double(ans),
		np, as.double(init.p),
		nhp, as.double(hp),
		as.double(0), as.double(matrix(0, np, np)), integer(1))
	r_optim <- list(p = r_optim[[6]], f = r_optim[[9]], h = matrix(r_optim[[10]], np), info = r_optim[[11]])

	init.p <- r_optim$p
	epsilons <- adjust_epsilons(epsilons, L, r_optim$h)

	st <- system.time({
		# 並行処理の準備
		library(parallel)
		library(doParallel)
		library(foreach)
		cl <- makeCluster(max(1, detectCores()))
		registerDoParallel(cl)

		r_sub <- foreach(i = 1:nchains) %dopar% {
			# dll <- paste("r_if_bayesian_logit", .Platform$dynlib.ext, sep="")
			# fs <- dyn.load(dll)

			r <- .Fortran("hmc_fit_bayesian_ologit",
				nrow(X), ncol(X), X, as.double(ans),
				as.integer(N + BI), matrix(as.double(0), N + BI, np),
				length(init.p), as.double(init.p),
				as.integer(nhp), as.double(hp),
				epsilons, as.integer(adjustSpan),
				as.integer(L),
				FALSE,
				as.double(diag(np)),
				FALSE,
				as.integer(seed[i]),
				as.double(0))

			# dyn.unload(dll)

			list(s = r[[6]], ar = r[[18]])
		}

		# 並行処理の終了
		stopCluster(cl)
	})

	colnames <- c(colnames(X), sprintf("gamma%d", 1:(nok_ans - 1)))
	r <- makeResultObj(
		"ologit", 
		r_sub, 
		BI, 
		colnames = colnames, 
		alpha = 0.05,
		r_optim,
		ar = mean(sapply(r_sub, "[[", 2)), 
		st = st)

	# 推定に使った変数
	r$input <- list(
		frml = frml,
		X = X,
		ans = ans,
		hp = hp,
		init.p = init.p
	)

	# MCMCで使ったパラメーター
	r$mcmc.param <- list(
		BI = BI,
		L = L,
		epsilons = epsilons,
		nchains = nchains,
		seed = seed)

	# 事後分布の要約統計量
	r$posterior <- list(
		beta = list(mu = mean(r)[1:ncol(X)], Sigma = vcov(r)[1:ncol(X), 1:ncol(X)]),
		gamma = list(mu = mean(r)[(1 + ncol(X)):np], Sigma = vcov(r)[(1 + ncol(X)):np, (1 + ncol(X)):np])
	)

	r
}

# 制約モデルのベイズファクター
log_BF.ologit <- function(r, index, value){
	sample <- as.matrix(r$mcmc.list)
	with(r$input, .Fortran("bayes_factor_ologit",
		nrow(X), ncol(X), as.double(X), as.double(ans),
		nrow(sample), ncol(sample), sample,
		as.integer(length(hp)), as.double(hp),
		as.integer(length(index)), as.integer(index), as.double(value),
		as.double(0))[[13]])
}

predict.ologit <- function(r, X = NULL, P = NULL){
	if(is.null(X)) X <- r$input$X
	if(is.null(P)) P <- as.matrix(r$mcmc.list)
	if(!is.matrix(X)) X <- matrix(X, 1, length(X))
	.Fortran("hmc_predict_ologit",
		nrow(X), ncol(X), as.double(X),
		nrow(P), ncol(P), as.double(P),
		double(nrow(X) * nrow(P)))[[7]]
}

### 混合ロジットで使う関数群 ###

#
# frml.mnl: モデル式（多項ロジット部分）
# frml.cnd: モデル式（条件付ロジット部分）
# data: データフレーム
# beta.mu: 線形予測部分のパラメーターの事前分布になる多変量正規分布の平均を表すハイパーパラメーター
# beta.Sigma: 〃 の分散共分散行列を 〃
# N: サンプリングする数
# BI: Burin-in（捨てる初期部分）
# adjustSpan: epsilonsを調整する期間/-1を入れると初期値の調整のみ行う
# nchains: 生成する乱数列の本数
# L, epsilons, seed: Hamiltonian Monte Carlo法のパラメーター
#
hmc.mlogit <- function(frml.mnl = NULL, frml.cnd = NULL, beta.mu, beta.Sigma, data = NULL, 
	N = 3000, BI = as.integer(N*0.2), adjustSpan = 0, L = 20, epsilons = NULL, nchains = 2, seed = NULL){

	seed = set_seeds(seed, nchains)

	if(is.null(frml.mnl)) stop("Specify a multinominal ligit model. If there is no variable, write as like 'dependent ~ 1'.")

	# 応答変数を抜き出す
	df02 <- model.frame(frml.mnl, data)
	y <- as.double(model.response(df02))
	nok <- as.integer(max(y)) # 選択肢の種類

	# 説明変数を整理
	X <- model.matrix(terms(df02), df02)
	if(!is.null(frml.cnd)){
		frml.cnd <- update(frml.cnd, ~ . + 0)
		Z <- model.matrix(frml.cnd, data)
		if(0 != ncol(Z) %% nok) stop("The number of kinds of responses, ", nok,", doesn't match to the formula of conditions: ", ncol(Z))
		noc <- as.integer(ncol(Z)/nok) # conditionalな変数の数
	} else {
		noc <- 0
		Z <- NULL
	}
	M <- cbind(X, Z)

	noev <- (nok - 1)*ncol(X) + noc
	if(noev != length(beta.mu)) stop("The length of beta.mu isn't same as the number of explanatory variables: ", noev)
	if(noev != ncol(beta.Sigma)) stop("The nummber of columns of beta.Sigma isn't same as the number of explanatory variables:", noev)
	if(noev != nrow(beta.Sigma)) stop("The nummber of rows of beta.Sigma isn't same as the number of explanatory variables:", noev)

	np <- as.integer(noc + ncol(X)*(nok - 1))
	init.p <- seq(-2.0, 2.0, length.out = np)

	hp <- c(beta.mu, chol2inv(chol(beta.Sigma)))
	nhp <- length(hp)

	# サンプリングされた行列につける列名
	options <- levels(model.response(df02))
	colnames <- character((nok - 1)*ncol(X) + noc)
	k <- 0
	if(2<=nok && 1<=ncol(X)) for(i in 2:nok){
		j <- 1 +(i - 2)*ncol(X)
		k <- j + ncol(X) - 1
		colnames[j:k] <- sprintf("%s.%s", colnames(X), options[i])
	}
	if(1<=noc) colnames[(1 + k):(noc + k)] <- sprintf("(Condition).%d", 1:noc)

	r_optim <- .Fortran("optim_bayesian_mlogit",
		nrow(M), ncol(M), as.double(M), as.double(y),
		np, as.double(init.p),
		nhp, as.double(hp),
		noc, nok,
		as.double(0), as.double(matrix(0, np, np)), integer(1))
	r_optim <- list(p = r_optim[[6]], f = r_optim[[11]], h = matrix(r_optim[[12]], np), info = r_optim[[13]])

	init.p <- r_optim$p
	epsilons <- adjust_epsilons(epsilons, L, r_optim$h)

	st <- system.time({
		# 並行処理の準備
		library(parallel)
		library(doParallel)
		library(foreach)
		cl <- makeCluster(max(1, detectCores()))
		registerDoParallel(cl)

		r_sub <- foreach(i = 1:nchains) %dopar% {
			# dll <- paste("r_if_bayesian_logit", .Platform$dynlib.ext, sep="")
			# fs <- dyn.load(dll)

			r <- .Fortran("hmc_fit_bayesian_mlogit",
				nrow(M), ncol(M), M, as.double(y),
				as.integer(N), matrix(as.double(0), N, np),
				as.integer(np), as.double(init.p),
				as.integer(nhp), as.double(hp),
				as.integer(noc), as.integer(nok),
				epsilons, as.integer(adjustSpan),
				as.integer(L),
				FALSE,
				as.double(diag(np)),
				FALSE,
				as.integer(seed[i]),
				as.double(0))

			# dyn.unload(dll)

			list(s = r[[6]], ar = r[[20]])
		}

		# 並行処理の終了
		stopCluster(cl)
	})

	# 予約統計量などをつくっておく
	r <- makeResultObj(
		"mlogit", 
		r_sub, 
		BI, 
		colnames = colnames, 
		alpha = 0.05,
		r_optim,
		ar = mean(sapply(r_sub, "[[", 2)), 
		st = st)

	# 推定に使った変数
	r$input <- list(
		frml.mnl = frml.mnl,
		frml.cnd = frml.cnd,
		M = M,
		X = X,
		Z = Z,
		y = y,
		nok = nok,
		noc = noc,
		hp = hp,
		init.p = init.p
	)

	# MCMCで使ったパラメーター
	r$mcmc.param <- list(
		BI = BI,
		L = L,
		epsilons = epsilons,
		nchains = nchains,
		seed = seed)

	# 事後分布の要約統計量
	r$posterior <- list(beta = list(mu = mean(r), Sigma = vcov(r)))

	r
}

# 制約モデルのベイズファクター
log_BF.mlogit <- function(r, index, value){
	sample <- as.matrix(r$mcmc.list)
	with(r$input, .Fortran("bayes_factor_mlogit",
		nrow(M), ncol(M), as.double(M), as.double(y),
		nrow(sample), ncol(sample), sample,
		as.integer(length(hp)), as.double(hp),
        as.integer(noc), as.integer(nok),
		as.integer(length(index)), as.integer(index), as.double(value),
		as.double(0))[[15]])
}

predict.mlogit <- function(r, Q = NULL, P = NULL){
	if(is.null(Q)) Q <- r$input$M
	if(is.null(P)) P <- as.matrix(r$mcmc.list)
	if(!is.matrix(Q)) Q <- matrix(Q, 1, length(Q))
	y <- with(r$input, .Fortran("hmc_predict_mlogit",
		nrow(Q), ncol(Q), as.double(Q),
		nrow(P), ncol(P), as.double(P),
        as.integer(noc), as.integer(nok),
		double(nrow(Q) * nrow(P) * nok)))[[9]]
	matrix(y, nrow(Q) * nrow(P), r$input$nok)
}

### ポアソン回帰で使う関数群 ###

hmc.poisson_exp <- function(frml, beta.mu, beta.Sigma, data = NULL, 
	N = 3000, BI = as.integer(N*0.1), adjustSpan = 0, L = 20, epsilons = NULL, nchains = 2, seed = NULL){

	seed = set_seeds(seed, nchains)

	data <- model.frame(frml, data)
	X <- model.matrix(frml, data)
	y <- model.response(data)

	if(!is.integer(y)) stop("Dependent variables must be integer.")
	if(any(y < 0)) stop("Dependent variables must be equal to or greater than zero.")

	init.p <- runif(ncol(X))
	np <- length(init.p)

	hp <- c(beta.mu, chol2inv(chol(beta.Sigma)))
	nhp <- length(hp)

	r_optim <- .Fortran("optim_bayesian_poisson_exp",
		nrow(X), ncol(X), as.double(X), as.double(y),
		np, as.double(init.p),
		nhp, as.double(hp),
		as.double(0), as.double(matrix(0, np, np)), integer(1))
	r_optim <- list(p = r_optim[[6]], f = r_optim[[9]], h = matrix(r_optim[[10]], np), info = r_optim[[11]])

	init.p <- r_optim$p
	epsilons <- adjust_epsilons(epsilons, L, r_optim$h)

	st <- system.time({
		# 並行処理の準備
		library(parallel)
		library(doParallel)
		library(foreach)
		cl <- makeCluster(max(1, detectCores()))
		registerDoParallel(cl)

		r_sub <- foreach(i = 1:nchains) %dopar% {

			r <- .Fortran("hmc_fit_poisson_exp",
				nrow(X), ncol(X), X, as.double(y),
				as.integer(N + BI), matrix(as.double(0), N + BI, np),
				length(init.p), as.double(init.p),
				as.integer(nhp), as.double(hp),
				epsilons, as.integer(adjustSpan),
				as.integer(L),
				FALSE,
				as.double(diag(np)),
				FALSE,
				as.integer(seed[i]),
				as.double(0))

			list(s = r[[6]], ar = r[[18]])
		}

		# 並行処理の終了
		stopCluster(cl)
	})

	# 予約統計量などをつくっておく
	r <- makeResultObj(
		"poisson", 
		r_sub, 
		BI, 
		colnames = colnames(X), 
		alpha = 0.05,
		r_optim,
		ar = mean(sapply(r_sub, "[[", 2)), 
		st = st)

	# 推定に使った変数
	r$input <- list(
		frml = frml,
		X = X,
		y = y,
		hp = hp,
		init.p = init.p
	)

	# MCMCで使ったパラメーター
	r$mcmc.param <- list(
		BI = BI,
		L = L,
		epsilons = epsilons,
		nchains = nchains,
		seed = seed)

	# 事後分布の要約統計量
	r$posterior <- list(beta = list(mu = mean(r), Sigma = vcov(r)))

	r
}

log_BF.poisson <- function(r, index, value){
	sample <- as.matrix(r$mcmc.list)
	with(r$input, .Fortran("bayes_factor_poisson_exp",
		nrow(X), ncol(X), as.double(X), as.double(y),
		nrow(sample), ncol(sample), sample,
		as.integer(length(hp)), as.double(hp),
		as.integer(length(index)), as.integer(index), as.double(value),
		as.double(0))[[13]])
}

predict.poisson <- function(r, X = NULL, P = NULL){
	if(is.null(X)) X <- r$input$X
	if(is.null(P)) P <- as.matrix(r$mcmc.list)
	if(!is.matrix(X)) X <- matrix(X, 1, length(X))
	.Fortran("hmc_predict_poisson_exp",
		nrow(X), ncol(X), as.double(X),
		nrow(P), ncol(P), as.double(P),
		double(nrow(X) * nrow(P)))[[7]]
}

# ユーザー定義目的関数からのサンプリング
#
# func: ユーザー定義関数
# ...: ユーザー定義関数にリレーするパラメーター
# data: データフレーム
# mu: 線形予測部分のパラメーターの事前分布になる多変量正規分布の平均を表すハイパーパラメーター
# Sigma: 〃 の分散共分散行列を 〃
# N: サンプリングする数
# BI: Burin-in（捨てる初期部分）
# adjustSpan: epsilonsを調整する期間/-1を入れると初期値の調整のみ行う
# nchains: 生成する乱数列の本数
# L, epsilons, seed: Hamiltonian Monte Carlo法のパラメーター
#
hmc.usrfunc <- function(init.p, func0, grad0 = NULL, ..., h = 1e-6, N = 3000, BI = as.integer(N*0.1), adjustSpan = 0, L = 20, epsilons = NULL, nchains = 2, seed = NULL, optim.param = NULL){
	seed = set_seeds(seed, nchains)
	func1 <- function(p) as.double(func0(p, ...))
	grad1 <- NULL
	if(!is.null(grad0)) grad1 <- function(p) as.double(grad0(p, ...))
	np <- length(init.p)
	hp <- numeric(0)
	if(is.null(optim.param)) optim.param = list(nbd = rep(0, np), u = rep(0.0, np), l = rep(0.0, np))
	r_optim <- .Call("invoke_optim", as.double(init.p), 
		func1, environment(func1), grad1, environment(grad1), h, 
		as.integer(optim.param$nbd), as.double(optim.param$u), as.double(optim.param$l))
	names(r_optim) <- c("p", "f", "h", "info")
	init.p <- r_optim[[1]]
	epsilons <- adjust_epsilons(epsilons, L, r_optim[[3]])

	st <- system.time({
		# 並行処理の準備
		library(parallel)
		library(doParallel)
		library(foreach)
		cl <- makeCluster(max(1, detectCores()))
		registerDoParallel(cl)

		r_sub <- foreach(i = 1:nchains) %dopar% {
			# foreachループの外側で定義したfunc1を呼ぼうとすると、
			# <anonymous>: ... may be used in an incorrect context: ‘func0(p, ...)’
			# と警告が出るので、内側で同じ関数を再定義する
			func2 <- function(p) as.double(func0(p, ...))
			grad2 <- NULL
			if(!is.null(grad0)) grad2 <- function(p) as.double(grad0(p, ...))

			r <- .Call("invoker",
				as.double(init.p), as.integer(N), 
				func2, environment(func2), grad2, environment(grad2), as.double(h),
				epsilons, as.integer(adjustSpan), 
				as.integer(L),
				FALSE,
				as.double(diag(np)),
				FALSE,
				as.integer(seed[i]))

			list(s = r[[1]], ar = r[[2]])
		}

		# 並行処理の終了
		stopCluster(cl)
	})

	# 予約統計量などをつくっておく
	r <- makeResultObj(
		"usrfunc", 
		r_sub, 
		BI, 
		colnames = NULL, 
		alpha = 0.05,
		r_optim,
		ar = mean(sapply(r_sub, "[[", 2)), 
		st = st)

	# 推定に使った変数
	r$input <- list(
		func0 = func0,
		init.p = init.p
	)

	# MCMCで使ったパラメーター
	r$mcmc.param <- list(
		BI = BI,
		L = L,
		epsilons = epsilons,
		nchains = nchains,
		seed = seed)

	r
}
