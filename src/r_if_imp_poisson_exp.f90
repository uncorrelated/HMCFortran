subroutine hmc_fit_imp_poisson_exp(nr, nc, nev, X, y, N, theta_sample, np, theta_init, nhp, hp, &
	epsilons, adjustEpsilonsN, L, randlength, Sigma, constrain, seed, accept_r)
	use poisson_exp_imp
	implicit none
	integer, intent(in) :: nr, nc, nev, np, nhp, N, adjustEpsilonsN, constrain, L, randlength
	double precision, dimension(np), intent(in) :: epsilons
	double precision, dimension(np, np), intent(in) :: sigma
	double precision, dimension(N, np), intent(out) :: theta_sample
	double precision, dimension(np), intent(inout) :: theta_init
	double precision, dimension(nr), intent(in) :: y
	double precision, dimension(nr, nc),intent(inout) :: X
	double precision, dimension(nhp) :: hp
	integer, intent(in) :: seed
	double precision, intent(out) :: accept_r
	type(poisson_imp) :: this

	call this%initialize(nr, nc, nev, X)
	call this%fit(nr, nc, X, y, N, theta_sample, np, theta_init, nhp, hp, &
		epsilons, adjustEpsilonsN, L, randlength, Sigma, constrain, seed, accept_r)
	call this%imp%finalize

end subroutine

subroutine optim_bayesian_imp_poisson_exp(nr, nc, nev, X, y, np, p, nhp, hp, f, h, info)
	use poisson_exp_imp
	implicit none
	integer, intent(in) :: nr, nc, nev, np, nhp
	double precision, dimension(nr), intent(in) :: y
	double precision, dimension(nr, nc), intent(inout) :: X
	double precision, dimension(np), intent(inout) :: p
	double precision, dimension(nhp), intent(in) :: hp
	double precision, intent(out) :: f
	double precision, dimension(np, np), intent(out) :: h
	integer, intent(out) :: info
	type(poisson_imp) :: this

	call this%initialize(nr, nc, nev, X)
	call this%optim(nr, nc, X, y, np, p, nhp, hp, f, h, info)
	call this%imp%finalize

end subroutine

subroutine log_ml_imp_poisson_exp(nr, nc, nev, X, y, np, p, nhp, hp, r, info)
	use poisson_exp_imp
	implicit none
	integer, intent(in) :: nr, nc, nev, np, nhp
	double precision, dimension(nr), intent(in) :: y
	double precision, dimension(nr, nc), intent(inout) :: X
	double precision, dimension(np), intent(in) :: p
	double precision, dimension(nhp), intent(in) :: hp
	double precision, intent(out) :: r
	integer, intent(out) :: info
	type(poisson_imp) :: this

	call this%initialize(nr, nc, nev, X)
	call this%la_log_ml(nr, nc, X, y, np, p, nhp, hp, r, info)
	call this%imp%finalize

end subroutine

subroutine bayes_factor_imp_poisson_exp(nr, nc, nev, X, y, ns, np, sample, nhp, hp, ncnst, pcnst, cnst, r)
	use poisson_exp_imp
	implicit none
	! ns: サンプルサイズ, ncnst: 制約の数
	integer, intent(in) :: nr, nc, nev, ns, np, nhp, ncnst
	double precision, dimension(nr), intent(in) :: y
	double precision, dimension(nr, nc), intent(inout) :: X
	double precision, dimension(ns, np), intent(in) :: sample
	double precision, dimension(nhp), intent(in) :: hp
	integer, dimension(ncnst), intent(in) :: pcnst
	double precision, dimension(ncnst), intent(in) :: cnst ! 制約条件の位置, 成約
	double precision, intent(out) :: r
	type(poisson_imp) :: this

	call this%initialize(nr, nc, nev, X)
	r = this%log_sd_bf(nr, nc, X, y, ns, np, sample, nhp, hp, ncnst, pcnst, cnst)
	call this%imp%finalize

end subroutine

subroutine hmc_predict_imp_poisson_exp(nr, nc, nev, X, ss, np, P, y)
	use poisson_exp_imp
	implicit none
	integer, intent(in) :: nr, nc, nev, np, ss
	double precision, dimension(nr, nc), intent(inout) :: X
	double precision, dimension(ss, np), intent(in) :: P
	double precision, dimension(nr * ss), intent(out) :: y
	type(poisson_imp) :: this

	call this%initialize(nr, nc, nev, X)
	call this%predict(nr, nc, X, ss, np, P, y)
	call this%imp%finalize

end subroutine
