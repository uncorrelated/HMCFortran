module bayesian_ologit
	use mvnorm
	use hmc
	implicit none
	type, extends(model_bf) :: ologit
	integer :: nok_ans ! Number Of Kinds of ANSwers, using only lg_marginal_prior
	contains
		procedure :: objf ! 目的関数
		procedure :: objfg ! objfのグラディエントを計算する関数
		procedure :: plimit ! パラメーターの正負の符号
		procedure :: lg_marginal_prior ! 対数化事前確率
		procedure :: predict ! 予測値
	end type

	contains
		double precision function objf(this, nr, nc, X, y, np, p, nhp, hp)
			implicit none
			class(ologit), intent(inout) :: this
			integer, intent(in) :: nr, nc, np, nhp
			double precision, dimension(nr, nc), intent(inout) :: X
			double precision, dimension(nr), intent(in) :: y
			double precision, dimension(np), intent(in) :: p
			double precision, dimension(nhp), intent(in) :: hp
			integer :: i, j, nk, r
			double precision :: Xb, s
			double precision, dimension(nc) :: beta_mu, beta
			double precision, dimension(np - nc) :: k, gamma_mu
			double precision, dimension(nc, nc) :: beta_sigma
			double precision, dimension(np - nc, np - nc) :: gamma_sigma

			nk = np - nc

			! ハイパーパラメーターを展開
			beta_mu = hp(1:nc)
			beta_sigma = reshape(hp(1 + nc:nc + nc**2), (/nc, nc/))
			gamma_mu = hp(1 + nc + nc**2:nc + nc**2 + nk)
			gamma_sigma = reshape(hp(nhp - nk**2 + 1:nhp), (/nk, nk/))

			! パラメーターを展開
			beta = p(1:nc)
			k(1:nk) = p(1 + nc:np)

			s = 0d0

			!$omp parallel
			!$omp do private(r, Xb) reduction(+: s)
			do i = 1, nr
				r = dint(y(i))
				Xb = sum(X(i, :) * beta)
				if(1 == r) then
					s = s - dlog(1 + dexp(Xb - k(r)))
				else if(nk + 1 == r) then
					s = s + dlog(1 - 1/(1 + dexp(Xb - k(r - 1))))
				else
					s = s + dlog(1/(1 + dexp(Xb - k(r))) - 1/(1 + dexp(Xb - k(r - 1))))
				end if
			end do
			!$omp end do
			!$omp end parallel

			objf = s &
				+ lg_dmvnorm_by_precision(nc, beta, beta_mu, beta_sigma) &
				+ lg_dmvnorm_by_precision(nk, k, gamma_mu, gamma_sigma)

		end function

		subroutine objfg(this, nr, nc, X, y, np, p, nhp, hp, g)
			implicit none
			class(ologit), intent(inout) :: this
			integer, intent(in) :: nr, nc, np, nhp
			double precision, dimension(nr, nc), intent(inout) :: X
			double precision, dimension(nr), intent(in) :: y
			double precision, dimension(np), intent(in) :: p
			double precision, dimension(nhp), intent(in) :: hp
			double precision, dimension(np), intent(out) :: g
			integer :: i, r, nk
			double precision :: Xb, A, B
			double precision, dimension(nc) :: beta_g_dmvnorm
			double precision, dimension(nc) :: beta_mu, beta
			double precision, dimension(nc, nc) :: beta_sigma
			double precision, dimension(np - nc) :: k, gamma_mu, gamma_g_dmvnorm
			double precision, dimension(np - nc, np - nc) :: gamma_sigma

			! write(*, *) "beta_mu", 1, nc
			! write(*, *) "beta_sigma", 1 + nc, nc + nc**2
			! write(*, *) "gammma_mu", 1 + nc + nc**2, nc + nc**2 + nk
			! write(*, *) "gammma_sigma", 1 + nc + nc**2 + nk, nc + nc**2 + nk + nk**2

			! ハイパーパラメーターを展開
			beta_mu = hp(1:nc)
			beta_sigma = reshape(hp(1 + nc:nc + nc**2), (/nc, nc/))

			nk = np - nc
			gamma_mu = hp(1 + nc + nc**2:nc + nc**2 + nk)
			gamma_sigma = reshape(hp(nhp - nk**2 + 1:nhp), (/nk, nk/))

			! write(*, *) "gamma_mu", gamma_mu
			! call rchkusr()

			! パラメーターを展開
			beta = p(1:nc)
			k(1:nk) = p(1 + nc:np)

			g = 0

			!$omp parallel
			!$omp do private(r, Xb, A, B) reduction(+: g )
			do i = 1, nr
				r = dint(y(i))
				Xb = sum(X(i, :) * beta)
				if(1 == r) then
					! (log(logit(Xb)))' = logit'(Xb)/logit(Xb)
					A = dexp(Xb - k(r))
					g(1:nc) = g(1:nc) - (1 + A)**(-2) * A * X(i, :) / (1 + A)**(-1)
					g(nc + r) = g(nc + r) + (1 + A)**(-2) * A / (1 + A)**(-1)
				else if(nk + 1 == r) then
					B = dexp(Xb - k(r - 1))
					g(1:nc) = g(1:nc) + (1 + B)**(-2) * B * X(i, :) / (1 - (1 + B)**(-1))
					g(nc + r - 1) = g(nc + r - 1) - (1 + B)**(-2) * B / (1 - (1 + B)**(-1))
				else
					A = dexp(Xb - k(r))
					B = dexp(Xb - k(r - 1))
					g(1:nc) = g(1:nc) & 
						- ((1 + A)**(-2) * A * X(i, :) - (1 + B)**(-2) * B * X(i, :)) & 
						/ ((1 + A)**(-1) - (1 + B)**(-1))
					g(nc + r) = g(nc + r) + (1 + A)**(-2) * A / ((1 + A)**(-1) - (1 + B)**(-1))
					g(nc + r - 1) = g(nc + r - 1) - (1 + B)**(-2) * B / ((1 + A)**(-1) - (1 + B)**(-1))
				end if
			end do
			!$omp end do
			!$omp end parallel

			call g_lg_dmvnorm_by_precision(nc, beta, beta_mu, beta_sigma, beta_g_dmvnorm)
			call g_lg_dmvnorm_by_precision(nk, k, gamma_mu, gamma_sigma, gamma_g_dmvnorm)
		
			g(1:nc) = g(1:nc)  + beta_g_dmvnorm
			g(1+nc:np) = g(1+nc:np)  + gamma_g_dmvnorm

		end subroutine

		! L-BFGS-Bの計算用
		subroutine plimit(this, np, nbd, l, u)
			implicit none
			class(ologit), intent(inout) :: this
			integer, intent(in) :: np
			integer, dimension(np), intent(out) :: nbd
			double precision, dimension(np), intent(out) :: l, u
			nbd = 0
		end subroutine

		! 事前確率の周辺分布を求める
		double precision function lg_marginal_prior(this, np, nhp, hp, ncnst, pcnst, cnst) result(r)
			implicit none
			class(ologit), intent(inout) :: this
			integer, intent(in) :: np, nhp, ncnst
			double precision, dimension(nhp), intent(in) :: hp
			integer, dimension(ncnst), intent(in) :: pcnst ! 制約条件の位置
			double precision, dimension(ncnst), intent(in) :: cnst ! 制約
			integer :: nk, nc, beta_ncnst, gamma_ncnst, i
			double precision, allocatable, dimension(:) :: beta_mu, gamma_mu, beta_cnst, gamma_cnst
			double precision, allocatable, dimension(:, :) :: beta_sigma, gamma_sigma
			integer, allocatable, dimension(:) :: beta_pcnst, gamma_pcnst

			nc = np - this%nok_ans + 1
			nk = np - nc
			allocate(beta_mu(nc), gamma_mu(nk), beta_sigma(nc, nc), gamma_sigma(nk, nk))

			! ハイパーパラメーターを展開
			beta_mu = hp(1:nc)
			beta_sigma = reshape(hp(1 + nc:nc + nc**2), (/nc, nc/))

			gamma_mu = hp(1 + nc + nc**2:nc + nc**2 + nk)
			gamma_sigma = reshape(hp(nhp - nk**2 + 1:nhp), (/nk, nk/))

			! 制約を分離
			allocate(beta_cnst(nc), beta_pcnst(nc), gamma_cnst(nk), gamma_pcnst(nk))
			beta_ncnst = 0
			gamma_ncnst = 0
			do i = 1, ncnst
				if(pcnst(i) <= nc) then
					beta_ncnst = beta_ncnst + 1
					beta_pcnst(beta_ncnst) = pcnst(i)
					beta_cnst(beta_ncnst) = cnst(i)
				else
					gamma_ncnst = gamma_ncnst + 1
					gamma_pcnst(gamma_ncnst) = pcnst(i) - nc
					gamma_cnst(gamma_ncnst) = cnst(i)
				end if
			end do

			r = 0d0
			if(0 < beta_ncnst) then
				r = r + lg_marginal_prior_lp(nc, nc + nc**2, hp(1:nc + nc**2), beta_ncnst, beta_pcnst, beta_cnst)
			end if 
			if(0 < gamma_ncnst) then
				! write(*, *) "nk", nk
				! write(*, *) "hp", hp(1+nc + nc**2:nhp)
				! write(*, *) "gamma_ncnst", gamma_ncnst
				! write(*, *) "gamma_pcnst", gamma_pcnst
				! write(*, *) "gamma_cnst", gamma_cnst
				r = r + lg_marginal_prior_lp(nk, nk + nk**2, hp(1+nc + nc**2:nhp), gamma_ncnst, gamma_pcnst, gamma_cnst)
			end if

		end function

		! 予測値を計算
		! nr: 行数
		! nc: 列数
		! X: 説明変数の行列
		! ss: パラメーターの集合の数（サンプルサイズ,行数）
		! np: パラメーターの数（列数）
		! p: 予測値を出すのに使う（サンプリングされた）パラメーター
		! Y: 予測値（出力）
		subroutine predict(this, nr, nc, X, ss, np, P, Y)
			implicit none
			class(ologit), intent(inout) :: this
			integer, intent(in) :: nr, nc, np, ss
			double precision, dimension(nr, nc), intent(inout) :: X
			double precision, dimension(ss, np), intent(in) :: P
			double precision, dimension(nr * ss, np - nc + 1), intent(out) :: Y
			integer :: i, j, r, nk
			double precision :: Xb
			double precision, dimension(nc) :: beta
			double precision, dimension(np - nc) :: kappa

			nk = np - nc

			!$omp parallel
			!$omp do private(i, beta, kappa, Xb)
			do j = 1, nr
				do i = 1, ss
					beta = P(i, 1:nc)
					kappa(1:nk) = P(i, 1 + nc:np)
					Xb = sum(X(j, :) * beta)
					do r = 1, nk + 1
						if(1 == r) then
							Y(i + nr*(j - 1), r) = 1/(1 + dexp(Xb - kappa(1)))
						else if(nk + 1 == r) then
							Y(i + nr*(j - 1), r) = 1 - 1/(1 + dexp(Xb - kappa(r - 1)))
						else 
							Y(i + nr*(j - 1), r) = 1/(1 + dexp(Xb - kappa(r))) - 1/(1 + dexp(Xb - kappa(r - 1)))
						end if
					end do
				end do
			end do
			!$omp end do
			!$omp end parallel
	
		end subroutine
	end module
