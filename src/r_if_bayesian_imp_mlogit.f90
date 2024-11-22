subroutine hmc_fit_bayesian_imp_mlogit(nr, nc, nev, X, y, N, theta_sample, np, theta_init, nhp, hp, &
    noc, nok, &
    epsilons, adjustEpsilonsN, L, randlength, Sigma, constrain, seed, accept_r)
	use bayesian_imp_mlogit
	implicit none
	integer, intent(in) :: nr, nc, nev, np, nhp, N, adjustEpsilonsN, constrain, L, randlength, noc, nok
	double precision, dimension(np), intent(in) :: epsilons
	double precision, dimension(np, np), intent(in) :: sigma
	double precision, dimension(N, np), intent(inout) :: theta_sample
	double precision, dimension(np), intent(inout) :: theta_init
	double precision, dimension(nr), intent(in) :: y
	double precision, dimension(nr, nc),intent(inout) :: X
	double precision, dimension(nhp) :: hp
	integer, intent(in) :: seed
	double precision, intent(out) :: accept_r
	type(mlogit_imp) :: this

	call this%initialize(nr, nc, nev, X, noc, nok)
	call this%fit(nr, nc, X, y, N, theta_sample, np, theta_init, nhp, hp, &
		epsilons, adjustEpsilonsN, L, randlength, Sigma, constrain, seed, accept_r)
	call this%imp%finalize

end subroutine

subroutine optim_bayesian_imp_mlogit(nr, nc, nev, X, y, np, p, nhp, hp, noc, nok, f, h, info)
	use bayesian_imp_mlogit
	implicit none
	integer, intent(in) :: nr, nc, nev, np, nhp, noc, nok
	double precision, dimension(nr), intent(in) :: y
	double precision, dimension(nr, nc), intent(inout) :: X
	double precision, dimension(np), intent(inout) :: p
	double precision, dimension(nhp), intent(in) :: hp
	double precision, intent(out) :: f
	double precision, dimension(np, np), intent(out) :: h
	integer, intent(out) :: info
	type(mlogit_imp) :: this

	call this%initialize(nr, nc, nev, X, noc, nok)
	call this%optim(nr, nc, X, y, np, p, nhp, hp, f, h, info)
	call this%imp%finalize

end subroutine

subroutine log_ml_bayesian_imp_mlogit(nr, nc, nev, X, y, np, p, nhp, hp, noc, nok, r, info)
	use bayesian_imp_mlogit
	implicit none
	integer, intent(in) :: nr, nc, nev, np, nhp, noc, nok
	double precision, dimension(nr), intent(in) :: y
	double precision, dimension(nr, nc), intent(inout) :: X
	double precision, dimension(np), intent(in) :: p
	double precision, dimension(nhp), intent(in) :: hp
	double precision, intent(out) :: r
	integer, intent(out) :: info
	type(mlogit_imp) :: this

	call this%initialize(nr, nc, nev, X, noc, nok)
	call this%la_log_ml(nr, nc, X, y, np, p, nhp, hp, r, info)
	call this%imp%finalize

end subroutine

subroutine bayes_factor_imp_mlogit(nr, nc, nev, X, y, ns, np, sample, nhp, hp, &
	noc, nok, ncnst, pcnst, cnst, r)
	use bayesian_imp_mlogit
	implicit none
	! ns: サンプルサイズ, ncnst: 制約の数
	integer, intent(in) :: nr, nc, nev, ns, np, nhp, ncnst, noc, nok
	double precision, dimension(nr), intent(in) :: y
	double precision, dimension(nr, nc), intent(inout) :: X
	double precision, dimension(ns, np), intent(in) :: sample
	double precision, dimension(nhp), intent(in) :: hp
	integer, dimension(ncnst), intent(in) :: pcnst
	double precision, dimension(ncnst), intent(in) :: cnst ! 制約条件の位置, 成約
	double precision, intent(out) :: r
	type(mlogit_imp) :: this

	call this%initialize(nr, nc, nev, X, noc, nok)
	r = this%log_sd_bf(nr, nc, X, y, ns, np, sample, nhp, hp, ncnst, pcnst, cnst)
	call this%imp%finalize

end subroutine

subroutine hmc_predict_imp_mlogit(nr, nc, nev, X, ss, np, P, noc, nok, Y)
	use bayesian_imp_mlogit
	implicit none
	integer, intent(in) :: nr, nc, nev, np, ss, noc, nok
	double precision, dimension(nr, nc), intent(inout) :: X
	double precision, dimension(ss, np), intent(in) :: P
	double precision, dimension(nr * ss, np - nc + 1), intent(out) :: Y
	type(mlogit_imp) :: this

	call this%initialize(nr, nc, nev, X, noc, nok)
	call this%predict(nr, nc, X, ss, np, P, Y)
	call this%imp%finalize

end subroutine
