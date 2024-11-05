subroutine hmc_fit_bayesian_ologit(nr, nc, X, y, N, theta_sample, np, theta_init, nhp, hp, &
	epsilons, adjustEpsilonsN, L, randlength, Sigma, constrain, seed, accept_r)
	use bayesian_ologit
	implicit none
	integer, intent(in) :: nr, nc, np, nhp, N, adjustEpsilonsN, constrain, L, randlength
	double precision, dimension(np), intent(in) :: epsilons
	double precision, dimension(np, np), intent(in) :: sigma
	double precision, dimension(N, np), intent(inout) :: theta_sample
	double precision, dimension(np), intent(inout) :: theta_init
	double precision, dimension(nr), intent(in) :: y
	double precision, dimension(nr, nc),intent(inout) :: X
	double precision, dimension(nhp) :: hp
	integer, intent(in) :: seed
	double precision, intent(out) :: accept_r
	type(ologit) :: this

!	write(*, *) "mu", hp(1:nc)
!	write(*, *) "sigma", reshape(hp(1 + nc:nc + nc**2), (/nc, nc/))

	call this%fit(nr, nc, X, y, N, theta_sample, np, theta_init, nhp, hp, &
		epsilons, adjustEpsilonsN, L, randlength, Sigma, constrain, seed, accept_r)

end subroutine

subroutine optim_bayesian_ologit(nr, nc, X, y, np, p, nhp, hp, f, h, info)
	use bayesian_ologit
	implicit none
	integer, intent(in) :: nr, nc, np, nhp
	double precision, dimension(nr), intent(in) :: y
	double precision, dimension(nr, nc), intent(inout) :: X
	double precision, dimension(np), intent(inout) :: p
	double precision, dimension(nhp), intent(in) :: hp
	double precision, intent(out) :: f
	double precision, dimension(np, np), intent(out) :: h
	integer, intent(out) :: info
	type(ologit) :: this

	call this%optim(nr, nc, X, y, np, p, nhp, hp, f, h, info)
end subroutine

subroutine log_ml_bayesian_ologit(nr, nc, X, y, np, p, nhp, hp, r, info)
	use bayesian_ologit
	implicit none
	integer, intent(in) :: nr, nc, np, nhp
	double precision, dimension(nr), intent(in) :: y
	double precision, dimension(nr, nc), intent(inout) :: X
	double precision, dimension(np), intent(in) :: p
	double precision, dimension(nhp), intent(in) :: hp
	double precision, intent(out) :: r
	integer, intent(out) :: info
	type(ologit) :: this

	call this%la_log_ml(nr, nc, X, y, np, p, nhp, hp, r, info)

end subroutine

subroutine bayes_factor_ologit(nr, nc, X, y, ns, np, sample, nhp, hp, ncnst, pcnst, cnst, r)
	use bayesian_ologit
	implicit none
	! ns: サンプルサイズ, ncnst: 制約の数
	integer, intent(in) :: nr, nc, ns, np, nhp, ncnst
	double precision, dimension(nr), intent(in) :: y
	double precision, dimension(nr, nc), intent(inout) :: X
	double precision, dimension(ns, np), intent(in) :: sample
	double precision, dimension(nhp), intent(in) :: hp
	integer, dimension(ncnst), intent(in) :: pcnst
	double precision, dimension(ncnst), intent(in) :: cnst ! 制約条件の位置, 成約
	double precision, intent(out) :: r
	type(ologit) :: this

	this%nok_ans = nhp - nc
	r = this%log_sd_bf(nr, nc, X, y, ns, np, sample, nhp, hp, ncnst, pcnst, cnst)
end subroutine

subroutine hmc_predict_ologit(nr, nc, X, ss, np, P, Y)
	use bayesian_ologit
	implicit none
	integer, intent(in) :: nr, nc, np, ss
	double precision, dimension(nr, nc), intent(inout) :: X
	double precision, dimension(ss, np), intent(in) :: P
	double precision, dimension(nr * ss, np - nc + 1), intent(out) :: Y
	type(ologit) :: this

	call this%predict(nr, nc, X, ss, np, P, Y)

end subroutine
