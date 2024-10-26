subroutine hmc_fit_poisson_exp(nr, nc, X, y, N, theta_sample, np, theta_init, nhp, hp, &
	epsilons, adjustEpsilonsN, L, randlength, Sigma, constrain, seed, accept_r)
	use poisson_exp
	implicit none
	integer, intent(in) :: nr, nc, np, nhp, N, adjustEpsilonsN, constrain, L, randlength
	double precision, dimension(np), intent(in) :: epsilons
	double precision, dimension(np, np), intent(in) :: sigma
	double precision, dimension(N, np), intent(out) :: theta_sample
	double precision, dimension(np), intent(inout) :: theta_init
	double precision, dimension(nr), intent(in) :: y
	double precision, dimension(nr, nc),intent(in) :: X
	double precision, dimension(nhp) :: hp
	integer, intent(in) :: seed
	double precision, intent(out) :: accept_r
	type(poisson) :: this

	call this%fit(nr, nc, X, y, N, theta_sample, np, theta_init, nhp, hp, &
		epsilons, adjustEpsilonsN, L, randlength, Sigma, constrain, seed, accept_r)
end subroutine

subroutine optim_bayesian_poisson_exp(nr, nc, X, y, np, p, nhp, hp, f, h)
	use poisson_exp
	implicit none
	integer, intent(in) :: nr, nc, np, nhp
	double precision, dimension(nr), intent(in) :: y
	double precision, dimension(nr, nc), intent(in) :: X
	double precision, dimension(np), intent(inout) :: p
	double precision, dimension(nhp), intent(in) :: hp
	double precision, intent(out) :: f
	double precision, dimension(np, np), intent(out) :: h
	type(poisson) :: this

	call this%optim(nr, nc, X, y, np, p, nhp, hp, f, h)
end subroutine

subroutine log_ml_poisson_exp(nr, nc, X, y, np, p, nhp, hp, r)
	use poisson_exp
	implicit none
	integer, intent(in) :: nr, nc, np, nhp
	double precision, dimension(nr), intent(in) :: y
	double precision, dimension(nr, nc), intent(in) :: X
	double precision, dimension(np), intent(in) :: p
	double precision, dimension(nhp), intent(in) :: hp
	double precision, intent(out) :: r
	type(poisson) :: this

	r = this%la_log_ml(nr, nc, X, y, np, p, nhp, hp)

end subroutine

subroutine bayes_factor_poisson_exp(nr, nc, X, y, ns, np, sample, nhp, hp, ncnst, pcnst, cnst, r)
	use poisson_exp
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
	type(poisson) :: this

	r = this%log_sd_bf(nr, nc, X, y, ns, np, sample, nhp, hp, ncnst, pcnst, cnst)
end subroutine

subroutine hmc_predict_poisson_exp(nr, nc, X, ss, np, P, y)
	use poisson_exp
	implicit none
	integer, intent(in) :: nr, nc, np, ss
	double precision, dimension(nr, nc), intent(in) :: X
	double precision, dimension(ss, np), intent(in) :: P
	double precision, dimension(nr * ss), intent(out) :: y
	type(poisson) :: this

	call this%predict(nr, nc, X, ss, np, P, y)

end subroutine
