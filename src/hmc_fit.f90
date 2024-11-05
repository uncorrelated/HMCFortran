module hmc
	use iso_c_binding
	implicit none
	! for sampling by MCMC
	type, abstract :: model
		double precision :: h = 1d-8
		procedure(hmc_fit), pointer, pass :: fit => hmc_fit
		procedure(maximum_likelihood_method), pointer, pass :: optim => maximum_likelihood_method
		procedure(laplace_log_marginal_likelihood), pointer, pass :: la_log_ml => laplace_log_marginal_likelihood
		! procedure(lg_savage_dickey_bayes_factors), pointer, pass :: log_sd_bf => lg_savage_dickey_bayes_factors
		contains
			procedure(objf_if), deferred, pass :: objf ! 目的関数
			procedure(objfg_if), deferred, pass :: objfg ! objfのグラディエントを計算する関数
			procedure(plimit_if), deferred, pass :: plimit ! パラメーターの正負の符号
			! procedure(lg_marginal_prior_if), deferred, pass :: lg_marginal_prior ! 対数化事前確率
	end type
	! for Bayes Factors/IWAMDE Savage-Dickey ratio
	type, abstract, extends(model) :: model_bf
		procedure(lg_savage_dickey_bayes_factors), pointer, pass :: log_sd_bf => lg_savage_dickey_bayes_factors
		contains
			procedure(lg_marginal_prior_if), deferred, pass :: lg_marginal_prior ! 対数化事前確率
	end type
    abstract interface
        double precision function objf_if(this, nr, nc, X, y, np, p, nhp, hp)
			import model
            implicit none
			class(model), intent(inout) :: this
            integer, intent(in) :: nr, nc, np, nhp
            double precision, dimension(nr, nc), intent(in) :: X
            double precision, dimension(nr), intent(in) :: y
            double precision, dimension(np), intent(in) :: p
            double precision, dimension(nhp), intent(in) :: hp
        end function
        subroutine objfg_if(this, nr, nc, X, y, np, p, nhp, hp, g)
			import model
			implicit none
			class(model), intent(inout) :: this
			integer, intent(in) :: nr, nc, np, nhp
			double precision, dimension(nr, nc), intent(in) :: X
			double precision, dimension(nr), intent(in) :: y
			double precision, dimension(np), intent(in) :: p
			double precision, dimension(nhp), intent(in) :: hp
			double precision, dimension(np), intent(out) :: g
        end subroutine 
        subroutine plimit_if(this, np, nbd, l, u)
			import model
            implicit none
			class(model), intent(inout) :: this
            integer, intent(in) :: np
            integer, dimension(np), intent(out) :: nbd
            double precision, dimension(np), intent(out) :: l, u
        end subroutine
		double precision function lg_marginal_prior_if(this, np, nhp, hp, ncnst, pcnst, cnst)
			import model_bf
			implicit none
			class(model_bf), intent(inout) :: this
			integer, intent(in) :: np, nhp, ncnst
			double precision, dimension(nhp), intent(in) :: hp
			integer, dimension(ncnst), intent(in) :: pcnst ! 制約条件の位置
			double precision, dimension(ncnst), intent(in) :: cnst ! 制約
		end function
	end interface
	contains
		! nr: Xの行数
		! nc: Xの列数
		! N: サンプリングする数
		! theta_sample: サンプル格納領域
		! np: パラメーターの数
		! theta_init: 初期パラメーター
		! nhp: ハイパーパラメーターの数
		! hp: ハイパーパラメーター
		! aux: 目的関数で使う補助情報
		! epsilons: ハミルトニアン・モンテ・カルロ法のステップ幅
		! adjustEpsilonsN: epsilonsを更新する回数 
		! L: 〃のステップ回数
		! randLength: Lを乱数で増減させる
		! Sigma: 〃の乱数生成用の分散共分散行列
		! constrain: 〃の非負制約
		! seed: 乱数のシード値
		! accept_r: 採択率
		subroutine hmc_fit(this, nr, nc, X, y, N, theta_sample, np, theta_init, nhp, hp, &
			epsilons, adjustEpsilonsN, L, randlength, Sigma, constrain, seed, accept_r)
		use IEEE_ARITHMETIC
		use dSFMT_interface ! 乱数生成モジュールの利用宣言
		implicit none
		class(model), intent(inout) :: this
		integer, intent(in) :: nr, nc, np, nhp, N, adjustEpsilonsN, L, randlength, constrain
		double precision, dimension(np), intent(in) :: epsilons
		double precision, dimension(np, np), intent(in) :: Sigma
		double precision, dimension(N, np), intent(inout) :: theta_sample
		double precision, dimension(np), intent(inout) :: theta_init
		double precision, dimension(nr), intent(in) :: y
		double precision, dimension(nr, nc),intent(in) :: X
		double precision, dimension(nhp), intent(in) :: hp
		double precision, dimension(np) :: g, mu_p
		integer, intent(in) :: seed
		double precision, intent(out) :: accept_r
		integer :: jj, r_L, noa = 0
		double precision, dimension(np, np) :: M_mx, Minv
		double precision, dimension(np) :: theta_new, r0, r_new
		type(dSFMT_t) :: rng ! 乱数生成器の状態をあらわす構造体
		integer :: min_size
		double precision, dimension(np) :: r_epsilons
		double precision :: u, log_numerator, log_denominator, log_accept_prob, accept_prob
		! for dual averagings
		double precision, dimension(np) :: da_eps_mean, da_eps, da_eps_init
		double precision :: H_mean
		integer :: da_t

		! 同時生成できる乱数の最小数
		min_size = dsfmt_get_min_array_size()

		! 最初の引数はシード値
		! 2番目の引数はバッファリングする乱数の数
		! 3番目の引数は乱数生成器の状態をあらわす構造体
		call dSFMT_init(seed, min_size * 16, rng)

		! # mass matrix
		mu_p = 0

		! # invert covariance Matrix for leapfrog
		! Σ = t(M_mx) %*% M_mx
		! Minv = M_mx⁻¹
		M_mx = Sigma
		call SQRT_S2(np, M_mx, Minv)

		jj = 1
		theta_sample(jj, :) = theta_init
		noa = 1

		! for dual averagings
		da_t  = 1
		H_mean = 0.0d0
		da_eps = epsilons

		! find reasonable epsilons
		! adjust the initial values of epsilons
		if(0 /= adjustEpsilonsN) then
			call find_reasonable_epsilons(da_eps)
		end if

		da_eps_init = da_eps
		da_eps_mean = da_eps

		do jj = 2, N

			call rchkusr()

			theta_new = theta_sample(jj - 1, :)
			theta_sample(jj, :) = theta_new
			call get_rand_arr_gaussian(rng, r0, np)
			r0 = MATMUL(M_mx, r0) + mu_p
			r_new = r0

			! # randomize epsilon and L
			if(0 /= randlength) then
				call get_rand_arr_open_open(rng, r_epsilons, np)
				r_epsilons = (2d-1 * r_epsilons + 9d-1) + epsilons
				r_L = nint((0.5d0 + 1.5d0 * get_rand_open_open(rng)) * L)
			else
				r_epsilons = da_eps
				r_L = L
			end if

			call leapfrog(this, nr, nc, X, y, np, theta_new, nhp, hp, r_new, r_epsilons, Minv, constrain, r_L)

			! if(4 < jj) then
			! 	if(sum(theta_sample(jj - 1, :)) == sum(theta_sample(jj - 2, :)) &
			! 		.and. sum(theta_sample(jj - 2, :)) == sum(theta_sample(jj - 3, :)) &
			! 		.and. sum(theta_sample(jj - 3, :)) == sum(theta_sample(jj - 4, :)) &
			! 		.and. sum(theta_new) /= sum(theta_new)) then
			! 		write(*, *) "jj", jj
			! 		write(*, *) "theta", theta_sample(jj - 1, :)
			! 		write(*, *) "theta_new", theta_new
			! 		write(*, *) "r_new", r_new
			! 		write(*, *) "r_epsilons", r_epsilons
			! 		write(*, *) "Minv", Minv
			! 		write(*, *) "constrain", constrain
			! 		write(*, *) "r_L", r_L
			! 		write(*, *) "this%objf(nr, nc, X, y, np, theta_sample(jj - 1, :), nhp, hp)", &
			! 			this%objf(nr, nc, X, y, np, theta_sample(jj - 1, :), nhp, hp)
			! 		write(*, *) "sum(MATMUL(r0, Minv) * r0)", &
			! 			sum(MATMUL(r0, Minv) * r0)
			! 		write(*, *) "da_eps", da_eps
			! 	end if 
			! end if

			! NaNやInfが戻った場合は、尤度関数の定義域外なので、採択確率ゼロとする
			accept_prob = 0d0

			log_numerator = this%objf(nr, nc, X, y, np, theta_new, nhp, hp) - sum(MATMUL(r_new, Minv) * r_new)/2
			if(IEEE_IS_FINITE(log_numerator)) then
				log_denominator = this%objf(nr, nc, X, y, np, theta_sample(jj - 1, :), nhp, hp) - sum(MATMUL(r0, Minv) * r0)/2
				log_accept_prob = min(0d0, log_numerator - log_denominator)
				accept_prob = dexp(log_accept_prob)

				u = get_rand_open_open(rng)

				! u < min(1, numerator/denominator)か？
				if(u < accept_prob) then
					theta_sample(jj, :) = theta_new(:)
					noa = noa + 1
				end if
			end if

			if(adjustEpsilonsN <= da_t) then
				da_eps = da_eps_mean
			else if(adjustEpsilonsN > da_t) then
				call dual_averaging(da_t, da_eps_init, da_eps, da_eps_mean, H_mean, accept_prob)
				da_t = da_t + 1
			end if
	end do

		accept_r = dble(noa) / dble(N)

		call dSFMT_end(rng)

		contains
			subroutine leapfrog(this, nr, nc, X, y, np, p, nhp, hp, r, epsilons, Minv, constrain, L)
				implicit none
				class(model), intent(inout) :: this
				integer, intent(in) :: nr, nc, np, nhp, constrain, L
				double precision, dimension(np), intent(in) :: epsilons
				double precision, dimension(np, np), intent(in) :: Minv
				double precision, dimension(np), intent(inout) :: p, r
				double precision, dimension(nr), intent(in) :: y
				double precision, dimension(nr, nc),intent(in) :: X
				double precision, dimension(nhp), intent(in) :: hp
				double precision, dimension(np) :: g
				double precision :: s
				integer :: i, j

				do j = 1, L
				!   # gradient of log posterior for old theta
					call this%objfg(nr, nc, X, y, np, p, nhp, hp, g)
					! # first momentum update
					r(:) = r(:) + epsilons(:) / 2 * g(:)
					! # theta update
					p(:) = p(:) + epsilons(:) * MATMUL(Minv, r)
					! # check positive
					if(0 /= constrain) then
						do i = 1, nc 
							if(0 > p(i)) then
								r = -r
								p = -p
							end if
						end do
					end if
					! 最後だけ2回目のrの更新をしない
					if(j < L) then
						! # gradient of log posterior for new theta
						call this%objfg(nr, nc, X, y, np, p, nhp, hp, g)
						r(:) = r(:) + epsilons(:) / 2 * g(:)
					end if
				end do
			end subroutine

			subroutine find_reasonable_epsilons(da_eps)
				use IEEE_ARITHMETIC
				implicit none
				double precision, dimension(np), intent(inout) :: da_eps
				integer, parameter :: maxit = 100
				integer :: it
				logical :: is_NaN
				double precision :: a
				double precision :: new_log_numerator
	
				is_NaN = 0 == 0
				it = 0
				da_eps = epsilons

				do while(is_NaN) 
					call get_rand_arr_gaussian(rng, r0, np)
					r0 = MATMUL(M_mx, r0) + mu_p
			
					theta_new = theta_init
					r_new =  r0
					
					call leapfrog(this, nr, nc, X, y, np, theta_new, nhp, hp, r_new, da_eps, Minv, constrain, L)
			
					log_numerator = this%objf(nr, nc, X, y, np, theta_new, nhp, hp) - sum(MATMUL(r_new, Minv) * r_new)/2
					log_denominator = this%objf(nr, nc, X, y, np, theta_init, nhp, hp) - sum(MATMUL(r0, Minv) * r0)/2

					is_NaN = (.not. IEEE_IS_FINITE(log_numerator)) .or. (.not. IEEE_IS_FINITE(log_denominator))
					it = it + 1
					if (maxit < it) then
						call rexit("Either log_numerator or log_denominator is NaN!")
					end if
					call rchkusr()
				end do

				a = -1d0
				if(log_numerator - log_denominator > dlog(0.5d0)) then
					a = a + 2
				end if


				it = 0
				do while(a*(log_numerator - log_denominator) > -a*dlog(2.0d0))

					! write(*, *) "a*(log_numerator - log_denominator)", a*(log_numerator - log_denominator), "-a*dlog(2.0d0)", -a*dlog(2.0d0)

					call rchkusr()

					r_new =  r0
					da_eps = da_eps*(2**(a))
					theta_new = theta_init

! write(*, *) "da_eps", da_eps
! write(*, *) "log_numerator", log_numerator, "log_denominator", log_denominator, &
! 	"log_numerator/log_denominator", log_numerator - log_denominator , "->", dexp(log_numerator - log_denominator)

					call leapfrog(this, nr, nc, X, y, np, theta_new, nhp, hp, r_new, da_eps, Minv, constrain, L)

					new_log_numerator = this%objf(nr, nc, X, y, np, theta_new, nhp, hp) - sum(MATMUL(r_new, Minv) * r_new)/2
					! log_denominator = this%objf(nr, nc, X, y, np, theta_init, nhp, hp) - sum(MATMUL(r0, Minv) * r0)/2

					if(IEEE_IS_FINITE(new_log_numerator)) then
						log_numerator = new_log_numerator
					else 
						da_eps = da_eps*(2**(-a))
					end if

					it = it + 1
					if(maxit <= it) then
						! write(*, *) "it", it, "da_eps", da_eps
						call rexit("it reached maxit!")
						! exit
					end if
				end do

				! write(*, *) "da_eps", da_eps
				! write(*, *) "log_numerator", log_numerator, "log_denominator", log_denominator, &
				! 	"log_numerator/log_denominator", log_numerator - log_denominator , "->", dexp(log_numerator - log_denominator)
				
			end subroutine

			subroutine dual_averaging(t, mu, eps, eps_mean, H_mean, accept_prob)
				implicit none
				double precision, parameter :: delta = 0.5d0, gamma = 0.05d0, t0 = 10.0d0, kappa = 0.75d0
				integer, intent(inout) :: t
				double precision, dimension(np), intent(in):: mu
				double precision, dimension(np), intent(inout):: eps, eps_mean ! \bar{x}に該当, xに該当
				double precision, dimension(np) :: log_eps
				double precision, intent(inout):: H_mean ! ΣH
				double precision, intent(in) :: accept_prob ! acceptance probability
				double precision :: eta, w

				w = 1d0/(dble(t) + t0)
				H_mean = (1 - w)*H_mean + w*(delta - accept_prob)

				log_eps = dlog(10d0*mu) - dsqrt(dble(t))/gamma * H_mean

				eta = t**(-kappa)
				eps_mean = exp(eta * log_eps + (1 - eta) * dlog(eps_mean))

				eps = exp(log_eps)

			end subroutine

			subroutine SQRT_S2(N, S, SI)
				implicit none
				integer, intent(in) :: N
				double precision, dimension(N, N), intent(inout) :: S
				double precision, dimension(N, N), intent(out) :: SI
				double precision, dimension(N) :: W
				integer :: LDA
				integer :: LWORK
				double precision, allocatable, dimension(:) :: WORK
				double precision, dimension(1) :: dsize
				integer :: INFO, i, j
				double precision, dimension(N, N) :: D, DI, B
				double precision, parameter :: minimum = 1d-10
				double precision :: sign
				LDA = N
				LWORK = MAX(1, 3*N)
				allocate(WORK(LWORK))
				call DSYEV('V', 'U', N, S, LDA, W, WORK, LWORK, INFO)
				sign = 1d0
				if(0 > W(1)) then
					sign = -1d0
					W = -W
				end if
				D(:, :) = 0
				do i = 1, N
					D(i, i) = DSQRT(W(i))
					DI(i, i) = 1/D(i, i)
				end do
				S = sign * S
				B = MATMUL(S, DI)
				S = MATMUL(D, TRANSPOSE(S))
				do j = 1, N
					do i = 1, N
						if(minimum > dabs(B(i, j))) then
							SI(i, j) = 0d0
						else
							SI(i, j) = B(i, j)
						end if
						if(minimum > dabs(S(i, j))) then
							S(i, j) = 0d0
						end if
					end do
				end do
			end subroutine
			
	end subroutine

	subroutine maximum_likelihood_method(this, nr, nc, X, y, np, p_lbfgsb, nhp, hp, f, h, info)
		use lbfgsb_module
		use iso_fortran_env, only: output_unit
		implicit none
		class(model), intent(inout) :: this
		integer, intent(in) :: nr, nc, np, nhp
		double precision, dimension(nr), intent(in) :: y
		double precision, dimension(nr, nc),intent(in) :: X
		! in: 初期値 out: 最適値
		double precision, dimension(np), intent(inout) :: p_lbfgsb
		double precision, dimension(nhp), intent(in) :: hp
		integer, intent(out) :: info
		integer :: i
		double precision, parameter :: pi = 3.141592653589793115998
	!!! L-BFGS-Bで利用する変数
		integer, parameter :: m = 5, iprint = -1
		double precision, parameter :: factr = 1.0d+7, pgtol = 1.0d-8
		character(len=60) :: task, csave
		logical :: lsave(4)
		integer :: isave(44)
		double precision :: f
		double precision :: dsave(29)
		integer, allocatable :: nbd(:), iwa(:)
		double precision, allocatable :: l(:), u(:), wa(:)
	!!! Hessianの計算用
		double precision, dimension(np) :: p_a, p_b
		double precision, dimension(np) :: g, g_a, g_b
		double precision, dimension(np, np) :: h

		! L-BFGS-Bの計算用の領域を確保
		allocate ( nbd(np), l(np), u(np) )
		allocate ( iwa(3*np) )
		allocate ( wa(2*m*np + 5*np + 11*m*m + 8*m) )

		! 上限と下限を設定
		call this%plimit(np, nbd, l, u)

		f = -this%objf(nr, nc, X, y, np, p_lbfgsb, nhp, hp)
		call this%objfg(nr, nc, X, y, np, p_lbfgsb, nhp, hp, g)
		g = -g

		task = 'START'
		do while(task(1:2)=='FG'.or.task=='NEW_X'.or.task=='START')
		call setulb ( np, m, p_lbfgsb, l, u, nbd, f, g, factr, pgtol, &
						wa, iwa, task, iprint,&
						csave, lsave, isave, dsave, &
						iteration_file = 'driver1_output.txt' )
		if (task(1:2) == 'FG') then
	! 符号を逆にしておく
			f = -this%objf(nr, nc, X, y, np, p_lbfgsb, nhp, hp)
			call this%objfg(nr, nc, X, y, np, p_lbfgsb, nhp, hp, g)
			g = -g
			end if
		end do

		if(task(1:4) == 'CONV') then
			info = 0
		else if(task(1:4) == 'ABNO') then
			info = 1
		else ! task(1:5) == 'ERROR'
			info = -1
		end if

		! 符号を反転
		f = -f

		! ヘッセ行列の計算
		p_a = p_lbfgsb
		p_b = p_lbfgsb
		h = 0
		do i = 1, np
			p_a(i) = p_lbfgsb(i) + this%h
			p_b(i) = p_lbfgsb(i) - this%h
			call this%objfg(nr, nc, X, y, np, p_a, nhp, hp, g_a)
			call this%objfg(nr, nc, X, y, np, p_b, nhp, hp, g_b)
			h(i, :) = (g_a - g_b) / (2*this%h)
			p_a(i) = p_lbfgsb(i)
			p_b(i) = p_lbfgsb(i)
		end do

	end subroutine

	! Laplace Approximated Log Marginal-Likelihood  
	subroutine laplace_log_marginal_likelihood(this, nr, nc, X, y, np, init_p, nhp, hp, r, error_code)
		implicit none
		class(model), intent(inout) :: this
		integer, intent(in) :: nr, nc, np, nhp
		double precision, dimension(nr), intent(in) :: y
		double precision, dimension(nr, nc),intent(in) :: X
		double precision, dimension(np), intent(in) :: init_p
		double precision, dimension(nhp), intent(in) :: hp
		double precision, parameter :: pi = 3.141592653589793115998
		double precision, intent(out) :: r
		integer, intent(out) :: error_code
	!!! f: 最適値での値
		double precision :: f
		double precision, dimension(np) :: p
	!!! Hessian
		double precision, dimension(np, np) :: h, tmp_h
	!!! 実対称行列の固有値計算用
		double precision, dimension(np) :: W
		double precision, allocatable, dimension(:) :: WORK
		integer :: LWORK, INFO

		p = init_p
		call maximum_likelihood_method(this, nr, nc, X, y, np, p, nhp, hp, f, h, error_code)

		! ヘッセ行列の固有値を求める
		LWORK = MAX(1, 3*np)
		allocate(WORK(LWORK))
		tmp_h = h 
		call DSYEV('N', 'U', np, tmp_h, np, W, WORK, LWORK, INFO)

		! 対数周辺尤度の計算
		! write(*, *) "p_lbfgsb", p_lbfgsb
		! write(*, *) "objf", objf(nr, nc, X, y, np, p_lbfgsb, nhp, hp, aux)
		! write(*, *) "np/2*dlog(2*pi)", np/2*dlog(2*pi)
		! write(*, *) "h", h
		! write(*, *) "W", W
		! write(*, *) "1/PRODUCT(W)/2", 1/PRODUCT(W)/2
		! write(*, *) "np/2*dlog(dble(nr)", np/2*dlog(dble(nr))
		! Σ⁻¹ = Hessian, then |Σ⁻¹| = 1/product(eigenvalues of Hessian)
		r = f + np/2*dlog(2*pi) - 1/PRODUCT(W)/2 - np/2*dlog(dble(nr))
	end subroutine

	! Bayes Factors/IWAMDE Savage-Dickey ratio
	double precision function lg_savage_dickey_bayes_factors(this, nr, nc, X, y, ns, np, sample, nhp, hp, ncnst, pcnst, cnst) result(r)
		use mvnorm
		implicit none
		class(model_bf), intent(inout) :: this
		! ns: サンプルサイズ, ncnst: 制約の数
		integer, intent(in) :: nr, nc, ns, np, nhp, ncnst
		double precision, dimension(nr), intent(in) :: y
		double precision, dimension(nr, nc), intent(in) :: X
		double precision, dimension(ns, np), intent(in) :: sample
		double precision, dimension(nhp), intent(in) :: hp
		integer, dimension(ncnst), intent(in) :: pcnst ! 制約条件の位置
		double precision, dimension(ncnst), intent(in) :: cnst ! 制約
		double precision, dimension(np) :: mu_sample, cnst_sample
		double precision, dimension(np, np) :: S_sample
		double precision :: pd_sum, log_pd, w_sum, log_w_n, log_w_d, pd
		double precision :: ig_a, ig_b
		double precision, allocatable, dimension(:) :: m_mu, w_mu
		double precision, allocatable, dimension(:, :) :: m_Sigma, w_Sigma
		integer :: i, j, INFO
		integer :: nncnst ! 制約なしの数
		integer, dimension(:), allocatable :: pncnst ! 制約なしパラメーター
		double precision, allocatable, dimension(:) :: vncnst

		! MCMCのサンプルから、多変量正規分布の平均と分散を求める
		call calc_mvnorm_prameters(ns, np, sample, mu_sample, S_sample)

		! weightsの計算用にパラメーターをまとめる
		nncnst = np - ncnst
		allocate(pncnst(nncnst))
		call setdiff(np, ncnst, pcnst, nncnst, pncnst)
		allocate(w_mu(nncnst), w_Sigma(nncnst, nncnst), vncnst(nncnst))
		call calc_marginal_normd_param(np, mu_sample, S_sample, nncnst, pncnst, w_mu, w_Sigma)

		! write(*, *) "ncnst", ncnst
		! write(*, *) "pcnst", pcnst
		! write(*, *) "nncnst", nncnst
		! write(*, *) "pncnst", pncnst
		! write(*, *) "mu_sample", mu_sample
		! write(*, *) "w_mu", w_mu
		! write(*, *) "S_sample", S_sample
		! write(*, *) "w_Sigma", w_Sigma

		! MCMCのサンプルから、事後周辺確率（条件付き確率）を求める
		w_sum = 0
		pd_sum = 0
		do i = 1, ns
			cnst_sample = sample(i, :)
			do j = 1, ncnst
				cnst_sample(pcnst(j)) = cnst(j)
			end do
			! w(b1|b2,D) = w(b1,b2|D)/w(b2|D)
			! precision matrixを使うようにしたい
			log_w_n = lg_dmvnorm(np, sample(i, :), mu_sample, S_sample)
			call subset(np, sample(i, :), nncnst, pncnst, vncnst)
			log_w_d = lg_dmvnorm(nncnst, vncnst, w_mu, w_Sigma)
			! なぜか条件付き確率が、条件無し確率より大きく計算されるときがある
			! if(log_w_n > log_w_d) then
			! 	write(*, *) "log_w_n", log_w_n, "log_w_d", log_w_d , "n - d", log_w_n - log_w_d, exp(log_w_n - log_w_d)
			! 	write(*, *) "sample", sample(i, :)
			! 	write(*, *) "mu_sample", mu_sample, "S_sample", S_sample
			! 	write(*, *) "vncnst", vncnst
			! 	write(*, *) "w_mu", w_mu, "w_Sigma", w_Sigma
			! 	log_w_d = log_w_n
			! end if
			w_sum = w_sum + exp(log_w_n - log_w_d)
			pd = exp( &
			log_w_n - log_w_d &
			+ this%objf(nr, nc, X, y, np, cnst_sample, nhp, hp) &
			- this%objf(nr, nc, X, y, np, sample(i, :), nhp, hp))
			if(pd /= pd) then
				pd = 0
			end if
			pd_sum = pd_sum + pd
		end do

		! write(*, *) "pd_sum", pd_sum
		! write(*, *) "dlog(pd_sum)", dlog(pd_sum)
		! write(*, *) "dlog(dble(ns))", dlog(dble(ns))
		! write(*, *) "np", np
		! write(*, *) "nhp", nhp
		! write(*, *) "hp", hp
		! write(*, *) "ncnst", ncnst
		! write(*, *) "pcnst", pcnst
		! write(*, *) "cnst", cnst
		! write(*, *) "lg_marginal_prior(np, nhp, hp, ncnst, pcnst, cnst)", this%lg_marginal_prior(np, nhp, hp, ncnst, pcnst, cnst)

		r = dlog(pd_sum) - dlog(dble(ns)) - this%lg_marginal_prior(np, nhp, hp, ncnst, pcnst, cnst)
	end function

	! lg_savage_dickey_bayes_factorsで使う制約の整理用関数
	! モジュールのメンバー関数ではあるが、クラスのメンバー関数ではない
	subroutine cnst_sort(n, unsort_index, unsort_value, sort_index, sort_value)
		implicit none
		integer, intent(in) :: n
		integer, dimension(n), intent(in) :: unsort_index
		integer, dimension(n), intent(out) :: sort_index
		double precision, dimension(n), intent(in) :: unsort_value
		double precision, dimension(n), intent(out) :: sort_value
		integer :: index_tmp
		double precision :: value_tmp
		integer i, j

		sort_index = unsort_index
		sort_value = unsort_value
		do j = 1, n - 1
			do i = 1, n - j
				if(sort_value(i) > sort_value(i + 1)) then
					index_tmp = sort_index(i) 
					value_tmp = sort_value(i)
					sort_index(i) = sort_index(i + 1) 
					sort_value(i) = sort_value(i + 1)
					sort_index(i + 1) = index_tmp 
					sort_value(i + 1) = value_tmp
				end if
			end do
		end do
	end subroutine

	! lg_savage_dickey_bayes_factorsで使う制約の整理用関数（補定用）
	! モジュールのメンバー関数ではあるが、クラスのメンバー関数ではない
	subroutine cnst_split(n, index_th, index_s, value_s, n_a, index_a, value_a, n_b, index_b, value_b)
		integer, intent(in) :: n, index_th
		integer, intent(out) :: n_a, n_b
		integer, dimension(n), intent(in) :: index_s
		integer, dimension(n), intent(out) :: index_a, index_b
		double precision, dimension(n), intent(in) :: value_s
		double precision, dimension(n), intent(out) :: value_a, value_b
		integer i, j

		n_a = 0
		n_b = 0
		do i = 1, n
			if(index_s(i) > index_th) then
				exit
			end if
			index_a(i) = index_s(i)
			value_a(i) = value_s(i)
			n_a = i
		end do
		i = n_a + 1
		do while(i <= n)
			index_b(i) = index_s(i) - index_th + 1
			value_b(i) = value_s(i)
			n_b = i
			i = i + 1
		end do
	end subroutine

! linear predictorの部分だけ事前分布を持つ場合に呼び出す
! e.g. logit/mlogit/poisson
! （モジュールのメンバーだが、抽象クラスのメンバー外）
	double precision function lg_marginal_prior_lp(np, nhp, hp, ncnst, pcnst, cnst) result(r)
		use mvnorm
		implicit none
		integer, intent(in) :: np, nhp, ncnst
		double precision, dimension(nhp), intent(in) :: hp
		integer, dimension(ncnst), intent(in) :: pcnst ! 制約条件の位置
		double precision, dimension(ncnst), intent(in) :: cnst ! 制約
		double precision :: log_cpd
		double precision, dimension(np) :: mu
		double precision, dimension(np, np) :: Sigma
		integer :: INFO, len_mu
		double precision, allocatable, dimension(:) :: m_mu
		double precision, allocatable, dimension(:, :) :: m_Sigma

		if(ncnst < 1) then
			r = 0
			return
		end if

		! ハイパーパラメーターを展開
		mu = hp(1:np)
		Sigma = reshape(hp(1 + np:np + np**2), (/np, np/))

		! precision matrixをvcovに戻す
		call DPOTRF('U', np, Sigma, np, INFO) ! コレスキー分解
		call DPOTRI('U', np, Sigma, np, INFO) ! 逆行列

		! 条件付き事前確率の計算
		allocate(m_mu(ncnst), m_Sigma(ncnst, ncnst))
		call calc_marginal_normd_param(np, mu, Sigma, ncnst, pcnst, m_mu, m_Sigma)

		r = lg_dmvnorm(ncnst, cnst, m_mu, m_Sigma)

	end function

end module
