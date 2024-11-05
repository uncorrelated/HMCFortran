module rusrf
	use iso_c_binding
	use hmc
	integer, parameter :: max_np = 100
	type, extends(model) :: rusrfnc
		! ユーザー定義関数用
		! RからCに渡された配列を、そのままRに戻して高速化
		type(c_ptr) :: r_objf_ptr, r_objf_envir_ptr, r_objfg_ptr, r_objfg_envir_ptr
		integer, dimension(max_np) :: nbd
		double precision, dimension(max_np) :: u
		double precision, dimension(max_np) :: l
	contains
		procedure :: objf ! 目的関数
		procedure :: objfg ! objfのグラディエントを計算する関数
		procedure :: plimit ! パラメーターの正負の符号
		! procedure :: lg_marginal_prior ! 対数化事前確率
	end type
	contains
		double precision function objf(this, nr, nc, X, y, np, p, nhp, hp)
			implicit none
			class(rusrfnc), intent(inout) :: this
			interface
				REAL(c_double) function call_r_objf(np, p, r_fn_ptr, envir_ptr) bind(C, name = "call_r_objf")
					use iso_c_binding
					implicit none
					INTEGER(c_int), intent(in) :: np
					REAL(c_double), dimension(:), intent(in) :: p
					type(c_ptr), intent(in) :: r_fn_ptr, envir_ptr
				end function
			end interface
			integer, intent(in) :: nr, nc, np, nhp
			double precision, dimension(nr, nc), intent(inout) :: X
			double precision, dimension(nr), intent(in) :: y
			double precision, dimension(np), intent(in) :: p
			double precision, dimension(np), intent(in) :: hp

			objf = call_r_objf(np, p, this%r_objf_ptr, this%r_objf_envir_ptr)

       end function

		subroutine objfg(this, nr, nc, X, y, np, p, nhp, hp, g)
			implicit none
			class(rusrfnc), intent(inout) :: this
			interface
				subroutine call_r_objfg(np, p, r_fn_ptr, envir_ptr, g) bind(C, name = "call_r_objfg")
					use iso_c_binding
					implicit none
					INTEGER(c_int), intent(in) :: np
					REAL(c_double), dimension(:), intent(in) :: p
					type(c_ptr), intent(in) :: r_fn_ptr, envir_ptr
					REAL(c_double), dimension(:), intent(out) :: g
				end subroutine
			end interface
			integer, intent(in) :: nr, nc, np, nhp
			double precision, dimension(nr, nc), intent(inout) :: X
			double precision, dimension(nr), intent(in) :: y
			double precision, dimension(np), intent(in) :: p
			double precision, dimension(np), intent(in) :: hp
			double precision, dimension(np), intent(out) :: g
			double precision, dimension(np) :: q
			integer :: i
			double precision :: h, v_h, v_l

			if(C_ASSOCIATED(this%r_objfg_ptr)) then
				! 微分関数が設定されているときの処理
				call call_r_objfg(np, p, this%r_objfg_ptr, this%r_objfg_envir_ptr, g)
				return
			end if

			! 数値微分を行う
			h = this%h
			q(:) = p(:) ! 入力データはいじらない
			do i = 1, np
				q(i) = p(i) + h
				v_h = this%objf(nr, nc, X, y, np, q, nhp, hp)
				q(i) = p(i) - h
				v_l = this%objf(nr, nc, X, y, np, q, nhp, hp)
				g(i) = (v_h - v_l)/(2*h)
				q(i) = p(i) ! 戻す
			end do
	
        end subroutine

		! L-BFGS-Bの計算用
		subroutine plimit(this, np, nbd, l, u)
			implicit none
			class(rusrfnc), intent(inout) :: this
            integer, intent(in) :: np
            integer, dimension(np), intent(out) :: nbd
			double precision, dimension(np), intent(out) :: l, u

			nbd = this%nbd(1:np)
			u = this%u(1:np)
			l = this%l(1:np)
		end subroutine

		! double precision function lg_marginal_prior(this, np, nhp, hp, ncnst, pcnst, cnst)
		! 	implicit none
		! 	class(rusrfnc), intent(inout) :: this
		! 	integer, intent(in) :: np, nhp, ncnst
		! 	double precision, dimension(nhp), intent(in) :: hp
		! 	integer, dimension(ncnst), intent(in) :: pcnst ! 制約条件の位置
		! 	double precision, dimension(ncnst), intent(in) :: cnst ! 制約
		! 	! R側に事前分布を求める関数が要る。現状はスタブ。
		! 	lg_marginal_prior = 1d0
		! end function
end module
