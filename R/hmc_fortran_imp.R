### 補定処理あり ###
hmc.blm.imp <- function(frml, ig.alpha, ig.beta, mu, Sigma, X_mu, X_mu_Sigma, frml.aux = NULL, data = NULL, 
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

	if(length(mu) != nev) stop("The length of mu is wrong.")
	if(!is.matrix(Sigma)) stop("Sigma must be a matrix.")
	if(nev != ncol(Sigma) || nev != nrow(Sigma)) stop("The length of Sigma is wrong.")
	if(length(X_mu) != n_X_mu) stop("The length of n_X_mu is wrong.")
	if(!is.matrix(X_mu_Sigma)) stop("X_mu_Sigma must be a matrix.")
	if(n_X_mu != ncol(X_mu_Sigma) || n_X_mu != nrow(X_mu_Sigma)) stop("The length of X_mu_Sigma is wrong.")

	hp <- c(ig.alpha, ig.alpha, mu, chol2inv(chol(Sigma)), X_mu, chol2inv(chol(X_mu_Sigma)))
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
