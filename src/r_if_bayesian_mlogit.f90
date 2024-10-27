subroutine hmc_fit_bayesian_mlogit(nr, nc, X, y, N, theta_sample, np, theta_init, nhp, hp, &
    noc, nok, &
    epsilons, adjustEpsilonsN, L, randlength, Sigma, constrain, seed, accept_r)
	use bayesian_mlogit
	implicit none
	integer, intent(in) :: nr, nc, np, nhp, N, adjustEpsilonsN, constrain, L, randlength, noc, nok
	double precision, dimension(np), intent(in) :: epsilons
	double precision, dimension(np, np), intent(in) :: sigma
	double precision, dimension(N, np), intent(inout) :: theta_sample
	double precision, dimension(np), intent(inout) :: theta_init
	double precision, dimension(nr), intent(in) :: y
	double precision, dimension(nr, nc),intent(in) :: X
	double precision, dimension(nhp) :: hp
	integer, intent(in) :: seed
	double precision, intent(out) :: accept_r
	type(mlogit) :: this

    this%noc = noc
    this%nok = nok

	call this%fit(nr, nc, X, y, N, theta_sample, np, theta_init, nhp, hp, &
		epsilons, adjustEpsilonsN, L, randlength, Sigma, constrain, seed, accept_r)
end subroutine

subroutine optim_bayesian_mlogit(nr, nc, X, y, np, p, nhp, hp, noc, nok, f, h, info)
	use bayesian_mlogit
	implicit none
	integer, intent(in) :: nr, nc, np, nhp, noc, nok
	double precision, dimension(nr), intent(in) :: y
	double precision, dimension(nr, nc), intent(in) :: X
	double precision, dimension(np), intent(inout) :: p
	double precision, dimension(nhp), intent(in) :: hp
	double precision, intent(out) :: f
	double precision, dimension(np, np), intent(out) :: h
	integer, intent(out) :: info
	type(mlogit) :: this

	this%noc = noc
    this%nok = nok
    this%h = 1d-5

	call this%optim(nr, nc, X, y, np, p, nhp, hp, f, h, info)
end subroutine

subroutine log_ml_bayesian_mlogit(nr, nc, X, y, np, p, nhp, hp, noc, nok, r, info)
	use bayesian_mlogit
	implicit none
	integer, intent(in) :: nr, nc, np, nhp, noc, nok
	double precision, dimension(nr), intent(in) :: y
	double precision, dimension(nr, nc), intent(in) :: X
	double precision, dimension(np), intent(in) :: p
	double precision, dimension(nhp), intent(in) :: hp
	double precision, intent(out) :: r
	integer, intent(out) :: info
	type(mlogit) :: this

    this%noc = noc
    this%nok = nok

	call this%la_log_ml(nr, nc, X, y, np, p, nhp, hp, r, info)

end subroutine

subroutine bayes_factor_mlogit(nr, nc, X, y, ns, np, sample, nhp, hp, noc, nok, ncnst, pcnst, cnst, r)
	use bayesian_mlogit
	implicit none
	! ns: サンプルサイズ, ncnst: 制約の数
	integer, intent(in) :: nr, nc, ns, np, nhp, ncnst, noc, nok
	double precision, dimension(nr), intent(in) :: y
	double precision, dimension(nr, nc), intent(in) :: X
	double precision, dimension(ns, np), intent(in) :: sample
	double precision, dimension(nhp), intent(in) :: hp
	integer, dimension(ncnst), intent(in) :: pcnst
	double precision, dimension(ncnst), intent(in) :: cnst ! 制約条件の位置, 成約
	double precision, intent(out) :: r
	type(mlogit) :: this

    this%noc = noc
    this%nok = nok

	r = this%log_sd_bf(nr, nc, X, y, ns, np, sample, nhp, hp, ncnst, pcnst, cnst)
end subroutine

subroutine hmc_predict_mlogit(nr, nc, X, ss, np, P, noc, nok, Y)
	use bayesian_mlogit
	implicit none
	integer, intent(in) :: nr, nc, np, ss, noc, nok
	double precision, dimension(nr, nc), intent(in) :: X
	double precision, dimension(ss, np), intent(in) :: P
	double precision, dimension(nr * ss, nok), intent(out) :: Y
	type(mlogit) :: this

	this%noc = noc
    this%nok = nok

	call this%predict(nr, nc, X, ss, np, P, Y)

end subroutine
