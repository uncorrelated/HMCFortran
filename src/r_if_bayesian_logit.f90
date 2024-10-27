subroutine hmc_fit_bayesian_logit(nr, nc, X, y, N, theta_sample, np, theta_init, nhp, hp, &
	epsilons, adjustEpsilonsN, L, randlength, Sigma, constrain, seed, accept_r)
	use bayesian_logit
	implicit none
	integer, intent(in) :: nr, nc, np, nhp, N, adjustEpsilonsN, L, randlength, constrain
	double precision, dimension(np), intent(in) :: epsilons
	double precision, dimension(np, np), intent(in) :: sigma
	double precision, dimension(N, np), intent(inout) :: theta_sample
	double precision, dimension(np), intent(inout) :: theta_init
	double precision, dimension(nr), intent(in) :: y
	double precision, dimension(nr, nc),intent(in) :: X
	double precision, dimension(nhp) :: hp
	integer, intent(in) :: seed
	double precision, intent(out) :: accept_r
	type(logit) :: this

	call this%fit(nr, nc, X, y, N, theta_sample, np, theta_init, nhp, hp, &
		epsilons, adjustEpsilonsN, L, randlength, Sigma, constrain, seed, accept_r)
end subroutine

subroutine optim_bayesian_logit(nr, nc, X, y, np, p, nhp, hp, f, h, info)
	use bayesian_logit
	implicit none
	integer, intent(in) :: nr, nc, np, nhp
	double precision, dimension(nr), intent(in) :: y
	double precision, dimension(nr, nc), intent(in) :: X
	double precision, dimension(np), intent(inout) :: p
	double precision, dimension(nhp), intent(in) :: hp
	double precision, intent(out) :: f
	double precision, dimension(np, np), intent(out) :: h
	integer, intent(out) :: info
	type(logit) :: this

	call this%optim(nr, nc, X, y, np, p, nhp, hp, f, h, info)
end subroutine

subroutine log_ml_bayesian_logit(nr, nc, X, y, np, p, nhp, hp, r, info)
	use bayesian_logit
	implicit none
	integer, intent(in) :: nr, nc, np, nhp
	double precision, dimension(nr), intent(in) :: y
	double precision, dimension(nr, nc), intent(in) :: X
	double precision, dimension(np), intent(in) :: p
	double precision, dimension(nhp), intent(in) :: hp
	double precision, intent(out) :: r
	integer, intent(out) :: info
	type(logit) :: this

	this%h = 1d-5
	call this%la_log_ml(nr, nc, X, y, np, p, nhp, hp, r, info)

end subroutine

subroutine bayes_factor_logit(nr, nc, X, y, ns, np, sample, nhp, hp, ncnst, pcnst, cnst, r)
	use bayesian_logit
	implicit none
	! ns: サンプルサイズ, ncnst: 制約の数
	integer, intent(in) :: nr, nc, ns, np, nhp, ncnst
	double precision, dimension(nr), intent(in) :: y
	double precision, dimension(nr, nc), intent(in) :: X
	double precision, dimension(ns, np), intent(in) :: sample
	double precision, dimension(nhp), intent(in) :: hp
	integer, dimension(ncnst), intent(in) :: pcnst
	double precision, dimension(ncnst), intent(in) :: cnst ! 制約条件の位置, 成約
	double precision, intent(out) :: r
	type(logit) :: this

	r = this%log_sd_bf(nr, nc, X, y, ns, np, sample, nhp, hp, ncnst, pcnst, cnst)
end subroutine

subroutine hmc_predict_logit(nr, nc, X, ss, np, P, y)
	use bayesian_logit
	implicit none
	integer, intent(in) :: nr, nc, np, ss
	double precision, dimension(nr, nc), intent(in) :: X
	double precision, dimension(ss, np), intent(in) :: P
	double precision, dimension(nr * ss), intent(out) :: y
	type(logit) :: this

	call this%predict(nr, nc, X, ss, np, P, y)

end subroutine
