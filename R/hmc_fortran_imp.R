### 補定処理あり ###
hmc.blm.imp <- function(frml, ig.alpha, ig.beta, beta.mu, beta.Sigma, X_mu, X_mu_Sigma, frml.aux = NULL, data = NULL, 
	N = 3000, BI = as.integer(N*0.1), adjustSpan =  0, L = 20, epsilons = NULL, nchains = 2, seed = NULL){

	seed <- set_seeds(seed, nchains)

	mf <- model.frame(frml, data, na.action = na.pass)
	M <- X <- model.matrix(frml, mf)
	y <- model.response(mf)
	nev <- ncol(X) # 回帰式の説明変数の数
	if(!is.null(frml.aux)){
		# 補定用の変数を足す
		frml.aux <- update(frml.aux, ~ . + 0)
		A <- model.matrix(frml.aux, data)
		M <- cbind(X, A)
	}
	n_X_mu <- ncol(M) - 1 # 切片項は補定に用いない
	if(!any(is.na(X))) stop("There is no missing value.")

	init.p <- c(0.5, runif(nev), X_mu)
	np <- length(init.p)

	if(length(beta.mu) != nev) stop("The length of beta.mu is wrong.")
	if(!is.matrix(beta.Sigma)) stop("beta.Sigma must be a matrix.")
	if(nev != ncol(beta.Sigma) || nev != nrow(beta.Sigma)) stop("The length of beta.Sigma is wrong.")
	if(length(X_mu) != n_X_mu) stop("The length of n_X_mu is wrong.")
	if(!is.matrix(X_mu_Sigma)) stop("X_mu_Sigma must be a matrix.")
	if(n_X_mu != ncol(X_mu_Sigma) || n_X_mu != nrow(X_mu_Sigma)) stop("The length of X_mu_Sigma is wrong.")

	hp <- c(ig.alpha, ig.alpha, beta.mu, chol2inv(chol(beta.Sigma)), X_mu, chol2inv(chol(X_mu_Sigma)))
	nhp <- length(hp)

	r_optim <- .Fortran("optim_bayesian_imp_lm",
		nrow(M), ncol(M), nev, as.double(M), as.double(y),
		np, as.double(init.p),
		nhp, as.double(hp),
		as.double(0), as.double(matrix(0, np, np)), integer(1), NAOK = TRUE)
	r_optim <- list(p = r_optim[[7]], f = r_optim[[10]], h = matrix(r_optim[[11]], np), info = r_optim[[12]])
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

			r <- .Fortran("hmc_fit_bayesian_imp_lm",
				nrow(X), ncol(X), nev, X, y,
				as.integer(N + BI), matrix(as.double(0), N + BI, np),
				np, as.double(init.p),
				nhp, as.double(hp),
				epsilons, as.integer(adjustSpan),
				as.integer(L),
				FALSE,
				as.double(diag(np)),
				FALSE,
				as.integer(seed[i]),
				as.double(0),
				NAOK = TRUE)

			list(s = r[[7]], ar = r[[19]])
		}

		# 並行処理の終了
		stopCluster(cl)
	})

	r <- makeResultObj(
		"blm.imp",
		r_sub, 
		BI, 
		colnames = c("(SD)", colnames(X), sprintf("imp.%s", colnames(M)[2:ncol(M)])), 
		alpha = 0.05,
		r_optim,
		ar = mean(sapply(r_sub, "[[", 2)), 
		st = st)

	# 推定に使った変数など
	r$input <- list(
		frml = frml,
		X = X,
		y = y,
		frml.aux = frml.aux,
		M = M,
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
		beta = list(mu = mean(r)[2:ncol(X)], Sigma = vcov(r)[2:ncol(X), 2:ncol(X)]),
		igamma = list(mean = mean(r)[[1]], var = max(1, vcov(r)[1, 1])),
		X = list(mu = mean(r)[(2 + ncol(X)):np], Sigma = vcov(r)[(2 + ncol(X)):np, (2 + ncol(X)):np])
	)
	r$posterior$igamma$alpha <- with(r$posterior$igamma, mean^2 / var + 2)
	r$posterior$igamma$beta <- with(r$posterior$igamma, mean*(alpha - 1))

	r
}

# 制約モデルのベイズファクター
log_BF.blm.imp <- function(r, index, value){
	sample <- as.matrix(r$mcmc.list)
	with(r$input, .Fortran("bayes_factor_imp_lm",
		nrow(M), ncol(M), ncol(X), as.double(M), as.double(y),
		nrow(sample), ncol(sample), sample,
		as.integer(length(hp)), as.double(hp),
		as.integer(length(index)), as.integer(index), as.double(value),
		as.double(0), NAOK = TRUE)[[14]])
}

# 予測値
predict.blm.imp <- function(r, X = NULL, P = NULL, A = NULL){
	if(is.null(X)) M <- r$input$M
	else{
		if(!is.matrix(X)) X <- matrix(X, 1, length(X))
		M <- cbind(X, A)
	}
	if(is.null(P)) P <- as.matrix(r$mcmc.list)
	.Fortran("hmc_predict_imp_lm",
		nrow(M), ncol(M), ncol(r$input$X), as.double(M),
		nrow(P), ncol(P), as.double(P),
		double(nrow(M) * nrow(P)),
		NAOK = TRUE)[[7]]
}

hmc.blogit.imp <- function(frml, beta.mu, beta.Sigma, X_mu, X_mu_Sigma, frml.aux = NULL, data = NULL, 
	N = 3000, BI = as.integer(N*0.1), adjustSpan =  0, L = 20, epsilons = NULL, nchains = 2, seed = NULL){

	seed = set_seeds(seed, nchains)

	mf <- model.frame(frml, data, na.action = na.pass)
	X <- model.matrix(frml, mf)
	M <- X
	y <- model.response(mf)
	nev <- ncol(X) # 回帰式の説明変数の数
	if(!is.null(frml.aux)){
		# 補定用の変数を足す
		frml.aux <- update(frml.aux, ~ . + 0)
		A <- model.matrix(frml.aux, data)
		M <- cbind(X, A)
	}
	n_X_mu <- ncol(M) - 1 # 切片項は補定に用いない
	if(!any(is.na(X))) stop("There is no missing value.")

	if(length(beta.mu) != nev) stop("The length of mu is wrong.")
	if(!is.matrix(beta.Sigma)) stop("Sigma must be a matrix.")
	if(nev != ncol(beta.Sigma) || nev != nrow(beta.Sigma)) stop("The length of Sigma is wrong.")
	if(length(X_mu) != n_X_mu) stop("The length of n_X_mu is wrong.")
	if(!is.matrix(X_mu_Sigma)) stop("X_mu_Sigma must be a matrix.")
	if(n_X_mu != ncol(X_mu_Sigma) || n_X_mu != nrow(X_mu_Sigma)) stop("The length of X_mu_Sigma is wrong.")

	init.p <- c(beta.mu, X_mu)
	np <- length(init.p)

	hp <- c(beta.mu, chol2inv(chol(beta.Sigma)), X_mu, chol2inv(chol(X_mu_Sigma)))
	nhp <- length(hp)
	r_optim <- .Fortran("optim_bayesian_imp_logit",
		nrow(M), ncol(M), nev, as.double(M), as.double(y),
		np, as.double(init.p),
		nhp, as.double(hp),
		as.double(0), as.double(matrix(0, np, np)), integer(1), NAOK = TRUE)
	r_optim <- list(p = r_optim[[7]], f = r_optim[[10]], h = matrix(r_optim[[11]], np), info = r_optim[[12]], X = X, X.imp = matrix(r_optim[[4]], nrow(X)))

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

			r <- .Fortran("hmc_fit_bayesian_imp_logit",
				nrow(M), ncol(M), nev, M, as.double(y),
				as.integer(N + BI), matrix(as.double(0), N + BI, np),
				length(init.p), as.double(init.p),
				as.integer(nhp), as.double(hp),
				epsilons, as.integer(adjustSpan),
				as.integer(L),
				FALSE,
				as.double(diag(np)),
				FALSE,
				as.integer(seed[i]),
				as.double(0),
				NAOK = TRUE)

			list(s = r[[7]], ar = r[[19]])
		}

		# 並行処理の終了
		stopCluster(cl)
	})

	# 予約統計量などをつくっておく
	r <- makeResultObj(
		"blogit.imp", 
		r_sub, 
		BI, 
		colnames = c(colnames(X), sprintf("imp.%s", colnames(M)[2:ncol(M)])), 
		alpha = 0.05,
		r_optim,
		ar = mean(sapply(r_sub, "[[", 2)), 
		st = st)

	# 推定に使った変数
	r$input <- list(
		frml = frml,
		X = X,
		y = y,
		frml.aux = frml.aux,
		M = M,
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
		X = list(mu = mean(r)[(1 + ncol(X)):np], Sigma = vcov(r)[(1 + ncol(X)):np, (1 + ncol(X)):np]) 
	)

	r
}

# 制約モデルのベイズファクター
log_BF.blogit.imp <- function(r, index, value){
	sample <- as.matrix(r$mcmc.list)
	with(r$input, .Fortran("bayes_factor_imp_logit",
		nrow(M), ncol(M), ncol(X), as.double(M), as.double(y),
		nrow(sample), ncol(sample), sample,
		as.integer(length(hp)), as.double(hp),
		as.integer(length(index)), as.integer(index), as.double(value),
		as.double(0), 
		NAOK = TRUE)[[14]])
}

predict.blogit.imp <- function(r, X = NULL, P = NULL, A = NULL){
	if(is.null(X)) M <- r$input$M
	else{
		if(!is.matrix(X)) X <- matrix(X, 1, length(X))
		M <- cbind(X, A)
	}
	if(is.null(P)) P <- as.matrix(r$mcmc.list)
	.Fortran("hmc_predict_imp_logit",
		nrow(M), ncol(M), ncol(r$input$X), as.double(M),
		nrow(P), ncol(P), as.double(P),
		double(nrow(M) * nrow(P)), 
		NAOK = TRUE)[[7]]
}

hmc.poisson_exp.imp <- function(frml, beta.mu, beta.Sigma, X_mu, X_mu_Sigma, frml.aux = NULL, data = NULL, 
	N = 3000, BI = as.integer(N*0.1), adjustSpan =  0, L = 20, epsilons = NULL, nchains = 2, seed = NULL){

	seed = set_seeds(seed, nchains)
	mf <- model.frame(frml, data, na.action = na.pass)
	X <- model.matrix(frml, mf)
	M <- X
	y <- model.response(mf)
	if(!is.integer(y)) stop("Dependent variables must be integer.")
	if(any(y < 0)) stop("Dependent variables must be equal to or greater than zero.")
	nev <- ncol(X) # 回帰式の説明変数の数
	if(!is.null(frml.aux)){
		# 補定用の変数を足す
		frml.aux <- update(frml.aux, ~ . + 0)
		A <- model.matrix(frml.aux, data)
		M <- cbind(X, A)
	}
	n_X_mu <- ncol(M) - 1 # 切片項は補定に用いない
	if(!any(is.na(X))) stop("There is no missing value.")
	if(length(beta.mu) != nev) stop("The length of mu is wrong.")
	if(!is.matrix(beta.Sigma)) stop("Sigma must be a matrix.")
	if(nev != ncol(beta.Sigma) || nev != nrow(beta.Sigma)) stop("The length of Sigma is wrong.")
	if(length(X_mu) != n_X_mu) stop("The length of n_X_mu is wrong.")
	if(!is.matrix(X_mu_Sigma)) stop("X_mu_Sigma must be a matrix.")
	if(n_X_mu != ncol(X_mu_Sigma) || n_X_mu != nrow(X_mu_Sigma)) stop("The length of X_mu_Sigma is wrong.")

	init.p <- c(runif(ncol(X)), X_mu)
	np <- length(init.p)

	hp <- c(beta.mu, chol2inv(chol(beta.Sigma)), X_mu, chol2inv(chol(X_mu_Sigma)))
	nhp <- length(hp)

	r_optim <- .Fortran("optim_bayesian_imp_poisson_exp",
		nrow(M), ncol(M), nev, as.double(M), as.double(y),
		np, as.double(init.p),
		nhp, as.double(hp),
		as.double(0), as.double(matrix(0, np, np)), integer(1), 
		NAOK = TRUE)
	r_optim <- list(p = r_optim[[7]], f = r_optim[[10]], h = matrix(r_optim[[11]], np), info = r_optim[[12]])

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

			r <- .Fortran("hmc_fit_imp_poisson_exp",
				nrow(M), ncol(M), nev, M, as.double(y),
				as.integer(N + BI), matrix(as.double(0), N + BI, np),
				length(init.p), as.double(init.p),
				as.integer(nhp), as.double(hp),
				epsilons, as.integer(adjustSpan),
				as.integer(L),
				FALSE,
				as.double(diag(np)),
				FALSE,
				as.integer(seed[i]),
				as.double(0), 
				NAOK = TRUE)

			list(s = r[[7]], ar = r[[19]])
		}

		# 並行処理の終了
		stopCluster(cl)
	})

	# 予約統計量などをつくっておく
	r <- makeResultObj(
		"poisson.imp", 
		r_sub, 
		BI, 
		colnames = c(colnames(X), sprintf("imp.%s", colnames(M)[2:ncol(M)])), 
		alpha = 0.05,
		r_optim,
		ar = mean(sapply(r_sub, "[[", 2)), 
		st = st)

	# 推定に使った変数
	r$input <- list(
		frml = frml,
		X = X,
		y = y,
		frml.aux = frml.aux,
		M = M,
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
		X = list(mu = mean(r)[(1 + ncol(X)):np], Sigma = vcov(r)[(1 + ncol(X)):np, (1 + ncol(X)):np])
	)

	r
}

log_BF.poisson.imp <- function(r, index, value){
	sample <- as.matrix(r$mcmc.list)
	with(r$input, .Fortran("bayes_factor_imp_poisson_exp",
		nrow(M), ncol(M), ncol(X), as.double(M), as.double(y),
		nrow(sample), ncol(sample), sample,
		as.integer(length(hp)), as.double(hp),
		as.integer(length(index)), as.integer(index), as.double(value),
		as.double(0), 
		NAOK = TRUE)[[14]])
}

predict.poisson.imp <- function(r, X = NULL, P = NULL){
	if(is.null(X)) M <- r$input$M
	else{
		if(!is.matrix(X)) X <- matrix(X, 1, length(X))
		M <- cbind(X, A)
	}
	if(is.null(P)) P <- as.matrix(r$mcmc.list)
	.Fortran("hmc_predict_imp_poisson_exp",
		nrow(M), ncol(M), ncol(r$input$X), as.double(M),
		nrow(P), ncol(P), as.double(P),
		double(nrow(M) * nrow(P)), 
		NAOK = TRUE)[[7]]
}

hmc.ologit.imp <- function(frml, beta.mu, beta.Sigma, gamma.mu, gamma.Sigma, X_mu, X_mu_Sigma, frml.aux = NULL, data = NULL, 
	N = 3000, BI = as.integer(N*0.1), adjustSpan =  0, L = 20, epsilons = NULL, nchains = 2, seed = NULL){

	seed = set_seeds(seed, nchains)

	frml <- update(frml, ~ . -1)
	mf <- model.frame(frml, data, na.action = na.pass)
	ans <- model.response(mf)
	nok_ans <- length(levels(ans))
	X <- model.matrix(frml, mf)
	M <- X
	nev <- ncol(X) # 回帰式の説明変数の数
	if(!is.null(frml.aux)){
		# 補定用の変数を足す
		frml.aux <- update(frml.aux, ~ . + 0)
		A <- model.matrix(frml.aux, data)
		M <- cbind(X, A)
	}
	n_X_mu <- ncol(M)

	if(!any(is.na(X))) stop("There is no missing value.")
	if(length(beta.mu) != nev) stop("The length of mu is wrong.")
	if(!is.matrix(beta.Sigma)) stop("Sigma must be a matrix.")
	if(nev != ncol(beta.Sigma) || nev != nrow(beta.Sigma)) stop("The length of Sigma is wrong.")
	if(length(gamma.mu) != nok_ans - 1)  stop("The length of gamma.mu is wrong.")
	if(length(gamma.mu) != ncol(gamma.Sigma) || length(gamma.mu) != nrow(gamma.Sigma)) stop("The length of X_mu_Sigma is wrong.")
	if(length(X_mu) != n_X_mu) stop("The length of n_X_mu is wrong.")
	if(!is.matrix(X_mu_Sigma)) stop("X_mu_Sigma must be a matrix.")
	if(n_X_mu != ncol(X_mu_Sigma) || n_X_mu != nrow(X_mu_Sigma)) stop("The length of X_mu_Sigma is wrong.")

	init.p <- as.double(c(runif(ncol(X), -1, 1), seq(1, 2, length.out = nok_ans - 1), X_mu))
	np <- length(init.p)

	hp <- c(beta.mu, chol2inv(chol(beta.Sigma)), gamma.mu, chol2inv(chol(gamma.Sigma)), X_mu, chol2inv(chol(X_mu_Sigma)))
	nhp <- length(hp)

	r_optim <- .Fortran("optim_bayesian_imp_ologit",
		nrow(M), ncol(M), nev, as.double(M), as.double(ans),
		np, as.double(init.p),
		nhp, as.double(hp),
		as.double(0), as.double(matrix(0, np, np)), integer(1), NAOK = TRUE)
	r_optim <- list(p = r_optim[[7]], f = r_optim[[10]], h = matrix(r_optim[[11]], np), info = r_optim[[12]])

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

			r <- .Fortran("hmc_fit_bayesian_imp_ologit",
				nrow(M), ncol(M), nev, as.double(M), as.double(ans),
				as.integer(N + BI), matrix(as.double(0), N + BI, np),
				length(init.p), as.double(init.p),
				as.integer(nhp), as.double(hp),
				epsilons, as.integer(adjustSpan),
				as.integer(L),
				FALSE,
				as.double(diag(np)),
				FALSE,
				as.integer(seed[i]),
				as.double(0),
				NAOK = TRUE)

			list(s = r[[7]], ar = r[[19]])
		}

		# 並行処理の終了
		stopCluster(cl)
	})

	colnames <- c(colnames(X), sprintf("gamma%d", 1:(nok_ans - 1)), sprintf("imp.%s", colnames(M)))
	r <- makeResultObj(
		"ologit.imp", 
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
		frml.aux = frml.aux,
		M = M,
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
		gamma = list(mu = mean(r)[(1 + ncol(X)):(ncol(X) + nok_ans - 1)], Sigma = vcov(r)[(1 + ncol(X)):(ncol(X) + nok_ans - 1), (1 + ncol(X)):(ncol(X) + nok_ans - 1)]),
		X_mu = list(mu = mean(r)[(ncol(X) + nok_ans):np], Sigma = vcov(r)[(ncol(X) + nok_ans):np, (ncol(X) + nok_ans):np])
	)

	r
}

log_BF.ologit.imp <- function(r, index, value){
	sample <- as.matrix(r$mcmc.list)
	with(r$input, .Fortran("bayes_factor_imp_ologit",
		nrow(M), ncol(M), ncol(X), as.double(M), as.double(ans),
		nrow(sample), ncol(sample), sample,
		as.integer(length(hp)), as.double(hp),
		as.integer(length(index)), as.integer(index), as.double(value),
		as.double(0), 
		NAOK = TRUE)[[14]])
}

predict.ologit.imp <- function(r, X = NULL, P = NULL){
	if(is.null(X)) M <- r$input$M
	else{
		if(!is.matrix(X)) X <- matrix(X, 1, length(X))
		M <- cbind(X, A)
	}
	if(is.null(P)) P <- as.matrix(r$mcmc.list)
	.Fortran("hmc_predict_imp_ologit",
		nrow(M), ncol(M), ncol(r$input$X), as.double(M),
		nrow(P), ncol(P), as.double(P),
		double(nrow(M) * nrow(P) * (length(levels(r$input$ans)) - 1)), 
		NAOK = TRUE)[[7]]
}

hmc.mlogit.imp <- function(frml.mnl = NULL, frml.cnd = NULL, beta.mu, beta.Sigma, 
	X_mu, X_mu_Sigma, frml.aux = NULL, data = NULL, 
	N = 3000, BI = as.integer(N*0.2), adjustSpan =  0, L = 20, epsilons = NULL, nchains = 2, seed = NULL){

	seed = set_seeds(seed, nchains)

	if(is.null(frml.mnl)) stop("Specify a multinominal logit model. If there is no variable, write as like 'dependent ~ 1'.")

	# 応答変数を抜き出す
	mf <- model.frame(frml.mnl, data, na.action = na.pass)
	y <- as.double(model.response(mf))
	nok <- as.integer(max(y)) # 選択肢の種類

	# 説明変数を整理
	X <- model.matrix(terms(mf), mf)
	if(!is.null(frml.cnd)){
		frml.cnd <- update(frml.cnd, ~ . + 0)
		mf_cnd <- model.frame(frml.cnd, data, na.action = na.pass)
		Z <- model.matrix(terms(mf_cnd), mf_cnd)
		if(0 != ncol(Z) %% nok) stop("The number of kinds of responses, ", nok,", doesn't match to the formula of conditions: ", ncol(Z))
		noc <- as.integer(ncol(Z)/nok) # conditionalな変数の数
	} else {
		noc <- 0
		Z <- NULL
	}
	M <- cbind(X, Z)
	nev <- ncol(M) # 回帰式の説明変数の数
	if(!is.null(frml.aux)){
		# 補定用の変数を足す
		frml.aux <- update(frml.aux, ~ . + 0)
		mf_aux <- model.frame(frml.aux, data, na.action = na.pass)
		A <- model.matrix(terms(mf_aux), mf_aux)
		M <- cbind(M, A)
	}
	n_X_mu <- ncol(M) - 1 # 切片項は使わない
	ncoef <- as.integer((nok - 1)*ncol(X) + noc) # 回帰式の係数の数 

	if(!any(is.na(M))) stop("There is no missing value.")
	if(ncoef != length(beta.mu)) stop("The length of beta.mu isn't same as the number of explanatory variables: ", ncoef)
	if(!is.matrix(beta.Sigma)) stop("beta.Sigma must be a matrix.")
	if(ncoef != ncol(beta.Sigma)) stop("The nummber of columns of beta.Sigma isn't same as the number of explanatory variables:", ncoef)
	if(ncoef != nrow(beta.Sigma)) stop("The nummber of rows of beta.Sigma isn't same as the number of explanatory variables:", ncoef)
	if(length(X_mu) != n_X_mu) stop("The length of n_X_mu is wrong.")
	if(!is.matrix(X_mu_Sigma)) stop("X_mu_Sigma must be a matrix.")
	if(n_X_mu != ncol(X_mu_Sigma) || n_X_mu != nrow(X_mu_Sigma)) stop("The length of X_mu_Sigma is wrong.")

	np <- as.integer(ncoef + n_X_mu)
	# init.p <- seq(-2.0, 2.0, length.out = np)
	# init.p <- runif(np, -1, 1)
	init.p <- rep(0, np)

	hp <- c(beta.mu, chol2inv(chol(beta.Sigma)), X_mu, chol2inv(chol(X_mu_Sigma)))
	nhp <- length(hp)

	# サンプリングされた行列につける列名
	options <- levels(model.response(mf))
	colnames <- character((nok - 1)*ncol(X) + noc)
	k <- 0
	if(2<=nok && 1<=ncol(X)) for(i in 2:nok){
		j <- 1 +(i - 2)*ncol(X)
		k <- j + ncol(X) - 1
		colnames[j:k] <- sprintf("%s.%s", colnames(X), options[i])
	}
	if(1<=noc) colnames[(1 + k):(noc + k)] <- sprintf("(Condition).%d", 1:noc)
	colnames <- c(colnames, sprintf("imp.%s", colnames(M)[-1]))

	r_optim <- .Fortran("optim_bayesian_imp_mlogit",
		nrow(M), ncol(M), nev, as.double(M), as.double(y),
		np, as.double(init.p),
		nhp, as.double(hp),
		noc, nok,
		as.double(0), as.double(matrix(0, np, np)), integer(1), NAOK = TRUE)
	impM <- matrix(r_optim[[4]], nrow(M))
	colnames(impM) <- colnames(M)
	r_optim <- list(p = r_optim[[7]], f = r_optim[[12]], h = matrix(r_optim[[13]], np), info = r_optim[[14]], M = impM)

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

			r <- .Fortran("hmc_fit_bayesian_imp_mlogit",
				nrow(M), ncol(M), nev, M, as.double(y),
				as.integer(N + BI), matrix(as.double(0), N + BI, np),
				as.integer(np), as.double(init.p),
				as.integer(nhp), as.double(hp),
				as.integer(noc), as.integer(nok),
				epsilons, as.integer(adjustSpan),
				as.integer(L),
				FALSE,
				as.double(diag(np)),
				FALSE,
				as.integer(seed[i]),
				as.double(0), NAOK = TRUE)

			list(s = r[[7]], ar = r[[21]])
		}

		# 並行処理の終了
		stopCluster(cl)
	})

	# 予約統計量などをつくっておく
	r <- makeResultObj(
		"mlogit.imp", 
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
		frml.aux = frml.aux,
		M = M,
		X = X,
		Z = Z,
		ans = model.response(mf),
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
	r$posterior <- list(
		beta = list(mu = mean(r)[1:nev], Sigma = vcov(r)[1:nev, 1:nev]),
		X_mu = list(mu = mean(r)[(1 + nev):np], Sigma = vcov(r)[(1 + nev):np, (1 + nev):np])
	)

	r
}

log_BF.mlogit.imp <- function(r, index, value){
	sample <- as.matrix(r$mcmc.list)
	with(r$input, .Fortran("bayes_factor_imp_mlogit",
		nrow(M), ncol(M), ncol(X) + ncol(Z), as.double(M), as.double(ans),
		nrow(sample), ncol(sample), sample,
		as.integer(length(hp)), as.double(hp),
		as.integer(noc), as.integer(nok),
		as.integer(length(index)), as.integer(index), as.double(value),
		as.double(0), 
		NAOK = TRUE)[[16]])
}

predict.mlogit.imp <- function(r, X = NULL, P = NULL){
	if(is.null(X)) M <- r$input$M
	else{
		if(!is.matrix(X)) X <- matrix(X, 1, length(X))
		M <- cbind(X, A)
	}
	if(is.null(P)) P <- as.matrix(r$mcmc.list)
	.Fortran("hmc_predict_imp_mlogit",
		nrow(M), ncol(M), ncol(r$input$X) + ncol(r$input$Z), as.double(M),
		nrow(P), ncol(P), as.double(P),
		as.integer(r$input$noc), as.integer(r$input$nok),
		double(nrow(M) * nrow(P) * length(levels(r$input$ans))), 
		NAOK = TRUE)[[10]]
}
