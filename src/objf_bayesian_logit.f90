module bayesian_logit
	use mvnorm
	use hmc
	implicit none
	type, extends(model_bf) :: logit
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
		class(logit), intent(inout) :: this
		integer, intent(in) :: nr, nc, np, nhp
		double precision, dimension(nr, nc), intent(inout) :: X
		double precision, dimension(nr), intent(in) :: y
		double precision, dimension(nc), intent(in) :: p
		double precision, dimension(nhp), intent(in) :: hp
		integer :: i, j
		double precision :: sxp, s
		double precision, dimension(nc) :: mu
		double precision, dimension(nc, nc) :: Sigma

		! ハイパーパラメーターを展開
		mu = hp(1:nc)
		Sigma = reshape(hp(1 + nc:nc + nc**2), (/nc, nc/))

		s = 0d0
		!$omp parallel
		!$omp do private(sxp) reduction(+: s)
		do i = 1, nr
			sxp = sum(X(i, :) * p)
			s = s + sxp * (y(i) - 1) - dlog(1.0d0 + dexp(-1 * sxp))
		end do
		!$omp end do
		!$omp end parallel
		objf = s + lg_dmvnorm_by_precision(nc, p, mu, Sigma)
	end function

	subroutine objfg(this, nr, nc, X, y, np, p, nhp, hp, g)
		implicit none
		class(logit), intent(inout) :: this
		integer, intent(in) :: nr, nc, np, nhp
		double precision, dimension(nr, nc), intent(inout) :: X
		double precision, dimension(nr), intent(in) :: y
		double precision, dimension(np), intent(in) :: p
		double precision, dimension(nhp), intent(in) :: hp
		double precision, dimension(np), intent(out) :: g
		integer :: i
		double precision :: sxp
		double precision, dimension(nc) :: g_dmvnorm
		double precision, dimension(nc) :: mu
		double precision, dimension(nc, nc) :: Sigma

		! ハイパーパラメーターを展開
		mu = hp(1:nc)
		Sigma = reshape(hp(1 + nc:nc + nc**2), (/nc, nc/))

		g = 0

		!$omp parallel
		!$omp do private(sxp) reduction(+: g )
		do i = 1, nr
			sxp = sum(X(i, :) * p)
			g = g + X(i, :) * (y(i) - 1 + dexp(-sxp)/(1 + dexp(-sxp)))
		end do
		!$omp end do
		!$omp end parallel

		call g_lg_dmvnorm_by_precision(nc, p, mu, Sigma, g_dmvnorm)
		g = g + g_dmvnorm

	end subroutine

	! L-BFGS-Bの計算用
	subroutine plimit(this, np, nbd, l, u)
		implicit none
		class(logit), intent(inout) :: this
		integer, intent(in) :: np
		integer, dimension(np), intent(out) :: nbd
		double precision, dimension(np), intent(out) :: l, u
		nbd = 2
		l = -100
		u = 100
	end subroutine

	! 事前確率の周辺分布を求める
	double precision function lg_marginal_prior(this, np, nhp, hp, ncnst, pcnst, cnst)
		implicit none
		class(logit), intent(inout) :: this
		integer, intent(in) :: np, nhp, ncnst
		double precision, dimension(nhp), intent(in) :: hp
		integer, dimension(ncnst), intent(in) :: pcnst ! 制約条件の位置
		double precision, dimension(ncnst), intent(in) :: cnst ! 制約

		lg_marginal_prior = lg_marginal_prior_lp(np, nhp, hp, ncnst, pcnst, cnst)

	end function

	! 予測値を計算
	! nr: 行数
	! nc: 列数
	! X: 説明変数の行列
	! ss: パラメーターの集合の数（サンプルサイズ,行数）
	! np: パラメーターの数（列数）
	! p: 予測値を出すのに使う（サンプリングされた）パラメーター
	! y: 予測値（出力）
	subroutine predict(this, nr, nc, X, ss, np, P, y)
		implicit none
		class(logit), intent(inout) :: this
		integer, intent(in) :: nr, nc, np, ss
		double precision, dimension(nr, nc), intent(inout) :: X
		double precision, dimension(ss, np), intent(in) :: P
		double precision, dimension(nr * ss), intent(out) :: y
		integer :: i, j
		double precision :: Xb

		!$omp parallel
		!$omp do private(i, Xb)
		do j = 1, nr
			do i = 1, ss
				Xb = sum(X(j, :) * P(i, :))
				y(i + nr*(j - 1)) = exp(Xb)/(1 + exp(Xb))
			end do
		end do
		!$omp end do
		!$omp end parallel

	end subroutine
end module
