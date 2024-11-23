### すべての推定で共通に使う関数群 ###

set_seeds <- function(seed, nchains){
	if(is.null(seed)) return(round(runif(nchains, 9999.501, 99999.499)))
	if(length(seed) < nchains) stop(sprintf("specify %d values for seeds", nchains))
	seed
}

adjust_epsilons <- function(epsilons, L, H){
	if(!is.null(epsilons)) return(epsilons)
	w <- sqrt(diag(solve(-H)))
	as.double(w) / L
}

makeResultObj <- function(classname, r_sub, BI, 
	colnames = NULL, alpha = 0.05, 
	r_optim, 
	ar = NULL, st = NULL){
	# サブサンプルを一つに統合する
	mcmc.lst <- mcmc.list()
	nchains <- length(r_sub)
	for(i in 1:nchains){
		s <- r_sub[[i]]$s
		colnames(s) = colnames
		mcmc.lst[[i]] <- mcmc(s, start = 1 + BI, end = nrow(s))
	}
	s <- as.matrix(mcmc.lst)
	m <- apply(s, 2, mean)
	sd <- apply(s, 2, sd)
	q <- apply(s, 2, \(x) quantile(x, prob = c(0, alpha/2, 0.5, 1-alpha/2, 1), na.rm = TRUE)) 
	names(r_optim$p) <- colnames(q) <- names(m) <- colnames(s)
	calc_la <- function(f, h, nr, np){
		f + np/2*log(2*pi) - 1 / det(h) / 2 - np / 2 * log(nr)
	}
	ml <- calc_la(r_optim$f, r_optim$h, nrow(s), ncol(s))
	r <- list(
		mcmc.list = mcmc.lst,
		MAP = r_optim,
		mean = m,
		sd = sd,
		vcov = cov(s),
		quantile = t(q),
		accept.rate = ar, 
		system.time = st,
		marginal.likelihood = ml
	)
	class(r)  <- c(classname, "hmc_f", class(r))
	r
}

### S3メソッド ###
plot.hmc_f <- function(r, ...){
	plot(r$mcmc.list, ...)
}

mean.hmc_f <- function(r){
	r$mean
}

vcov.hmc_f <- function(r){
	r$vcov
}

median.hmc_f <- function(r){
	r$quantile[, 3]
}

quantile.hmc_f <- function(r, alpha = NULL){
	if(is.null(alpha)){
		return(r$quantile)
	} else {
		sample <- as.matrix(r$mcmc.list)
		q <- apply(sample, 2, \(x) quantile(x, prob = c(0, alpha/2, 0.5, 1-alpha/2, 1))) 
		colnames(q) <- names(r$mean)
		return(q)
	}
}

summary.hmc_f <- function(r, ...){
	summary(r$mcmc.list, ...)
}

HPDinterval.hmc_f <- function(r, ...){
	HPDinterval(r$mcmc.list)
}

log_ml <- function(obj) UseMethod("log_ml")
log_BF <- function(obj, ...) UseMethod("log_BF")
