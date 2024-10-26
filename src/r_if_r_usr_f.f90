subroutine hmc_fit_rusrf(N, theta_sample, np, theta_init, &
	r_objf_ptr, r_objf_envir_ptr, r_objfg_ptr, r_objfg_envir_ptr, h, &
	epsilons, adjustEpsilonsN, L, randlength, Sigma, constrain, seed, accept_r) &
	bind(C, name = "hmc_fit_rusrf")

	use iso_c_binding
	use rusrf
	implicit none
	integer, intent(in) :: np, N, adjustEpsilonsN, constrain, L, randlength
	double precision, dimension(np), intent(in) :: epsilons
	double precision, dimension(np, np), intent(in) :: sigma
	double precision, dimension(N, np), intent(inout) :: theta_sample
	double precision, dimension(np), intent(inout) :: theta_init
    type(c_ptr), intent(in) :: r_objf_ptr, r_objf_envir_ptr, r_objfg_ptr, r_objfg_envir_ptr
	double precision, intent(in) :: h
	integer, intent(in) :: seed
	double precision, intent(out) :: accept_r
    integer, parameter :: nr = 0, nc = 0, nhp = 0
	double precision, dimension(nr) :: dummy_y
	double precision, dimension(nr, np) :: dummy_X
	double precision, dimension(nhp) :: dummy_hp
	type(rusrfnc) :: this

	this%r_objf_envir_ptr = r_objf_envir_ptr
	this%r_objf_ptr = r_objf_ptr
	this%r_objfg_envir_ptr = r_objfg_envir_ptr
	this%r_objfg_ptr = r_objfg_ptr
	this%h = h

	call this%fit(np, nc, dummy_X, dummy_y, N, theta_sample, np, theta_init, nhp, dummy_hp, &
		epsilons, adjustEpsilonsN, L, randlength, Sigma, constrain, seed, accept_r)
end subroutine

subroutine optim_rusrf(np, p, r_objf_ptr, r_objf_envir_ptr, r_objfg_ptr, r_objfg_envir_ptr, nbd, u, l, h, f, hessian) &
	bind(C, name = "optim_rusrf")
	use iso_c_binding
	use rusrf
	implicit none
	integer, intent(in) :: np
	double precision, dimension(np), intent(inout) :: p
    type(c_ptr), intent(in) :: r_objf_ptr, r_objf_envir_ptr, r_objfg_ptr, r_objfg_envir_ptr
	double precision, intent(in) :: h
	double precision, intent(out) :: f
	double precision, dimension(np, np), intent(out) :: hessian
	integer, dimension(np) :: nbd
	double precision, dimension(np) :: u, l
    integer, parameter :: nr = 1, nc = 1, nhp = 1
	double precision, dimension(nr) :: dummy_y
	double precision, dimension(nr, np) :: dummy_X
	double precision, dimension(nhp) :: dummy_hp
	type(rusrfnc) :: this

	this%r_objf_envir_ptr = r_objf_envir_ptr
	this%r_objf_ptr = r_objf_ptr
	this%r_objfg_envir_ptr = r_objfg_envir_ptr
	this%r_objfg_ptr = r_objfg_ptr
	this%h = h
	this%nbd(1:np) = nbd
	this%u(1:np) = u
	this%l(1:np) = l

	call this%optim(nr, nc, dummy_X, dummy_y, np, p, nhp, dummy_hp, f, hessian)
end subroutine

subroutine log_ml_rusrf(nr, nc, np, p, r_objf_ptr, r_objf_envir_ptr, r_objfg_ptr, r_objfg_envir_ptr, nbd, u, l, h, r) &
	bind(C, name = "log_ml_rusrf")
	use iso_c_binding
	use rusrf
	implicit none
	integer, intent(in) :: nr, nc, np
	double precision, dimension(np), intent(in) :: p
	double precision, intent(out) :: r
    type(c_ptr), intent(in) :: r_objf_ptr, r_objf_envir_ptr, r_objfg_ptr, r_objfg_envir_ptr
	double precision, intent(in) :: h
	integer, dimension(np) :: nbd
	double precision, dimension(np) :: u, l
    integer, parameter :: nhp = 0
	double precision, dimension(nr) :: dummy_y
	double precision, dimension(nr, nc) :: dummy_X
	double precision, dimension(nhp) :: dummy_hp
	type(rusrfnc) :: this

	this%r_objf_envir_ptr = r_objf_envir_ptr
	this%r_objf_ptr = r_objf_ptr
	this%r_objfg_envir_ptr = r_objfg_envir_ptr
	this%r_objfg_ptr = r_objfg_ptr
	this%h = h
	this%nbd(1:np) = nbd
	this%u(1:np) = u
	this%l(1:np) = l

	r = this%la_log_ml(nr, nc, dummy_X, dummy_y, np, p, nhp, dummy_hp)

end subroutine
