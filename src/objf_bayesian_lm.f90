module bayesian_lm
	use mvnorm
	use hmc
	implicit none
	logical, parameter :: usePrecisionMatrix = .true.
	type, extends(model_bf) :: blm
	contains
		procedure :: objf ! 目的関数
		procedure :: objfg ! objfのグラディエントを計算する関数
		procedure :: plimit ! パラメーターの正負の符号
		procedure :: lg_marginal_prior ! 対数化事前確率
		procedure :: predict ! 予測値
	end type
	contains
		! nr: 行数
		! nc: 列数
		! X: 説明変数の行列
		! y: 被説明変数の列
		! p: nc+1個のパラメーター
		! hp: 2*nc + 2 + nc**2個のハイパー
		double precision function objf(this, nr, nc, X, y, np, p, nhp, hp)
			implicit none
			class(blm), intent(inout) :: this
			integer, intent(in) :: nr, nc, np, nhp
			double precision, dimension(nr, nc), intent(inout) :: X
			double precision, dimension(nr), intent(in) :: y
			double precision, dimension(nc + 1), intent(in) :: p
			double precision, dimension(2*nc + 2 + nc**2), intent(in) :: hp
			integer :: i, j
			double precision :: s
			double precision :: ig_a, ig_b
			double precision, dimension(nc) :: beta, mu
			double precision, dimension(nc, nc) :: Sigma

			! np: パラメーターの数
			! np = nc + 1
			! nhp: ハイパーパラメーターの数
			! nhp = 2*np + (np - 1)**2 ! = 2*nc + 2 + nc**2

			! ハイパーパラメーターを展開
			ig_a = hp(1)
			ig_b = hp(2)
			mu = hp(3:1 + np)
			Sigma = reshape(hp(2 + np:1 + np + nc**2), (/nc, nc/))

			beta = p(2: 1 + nc)
			s = 0

			!$omp parallel
			!$omp do reduction(+: s)
			do i = 1, nr
				s = s + ldnorm(y(i) - sum(X(i, :)*beta), p(1))
			end do
			!$omp end do
			!$omp end parallel

			if(.not. usePrecisionMatrix) then
				objf = s + lg_igamma(p(1), ig_a, ig_b) + lg_dmvnorm(nc, beta, mu, Sigma)
			else
				objf = s + lg_igamma(p(1), ig_a, ig_b) + lg_dmvnorm_by_precision(nc, beta, mu, Sigma)
			end if
		end function

		subroutine objfg(this, nr, nc, X, y, np, p, nhp, hp, g)
			implicit none
			class(blm), intent(inout) :: this
			integer, intent(in) :: nr, nc, np, nhp
			double precision, dimension(nr, nc), intent(inout) :: X
			double precision, dimension(nr), intent(in) :: y
			double precision, dimension(nc + 1), intent(in) :: p
			double precision, dimension(2*nc + 2 + nc**2), intent(in) :: hp
			double precision, dimension(nc + 1), intent(out) :: g
			double precision, dimension(nc + 1) :: q 
			double precision :: ssigma, dssigma, ig_a, ig_b, slev
			double precision, dimension(nc) :: beta, dbeta, mu, g_dmvnorm
			double precision, dimension(nr) :: res
			integer :: i
			double precision, dimension(nc, nc) :: Sigma

			! パラメーターを展開
			ssigma = p(1)
			beta(:) = p(2:nc + 1)

			! np: パラメーターの数
			! np = nc + 1
			! nhp: ハイパーパラメーターの数
			! nhp = 2*np + (np - 1)**2 ! = 2*nc + 2 + nc**2

			! ハイパーパラメーターを展開
			ig_a = hp(1)
			ig_b = hp(2)
			mu = hp(3:1 + np)
			Sigma = reshape(hp(2 + np:1 + np + nc**2), (/nc, nc/))

			! BLAS
			! res(:) = y(:)
			! call dgemm('n', 'n', nr, 1, nc, -1d0, X, nr, beta, nc, 1d0, res, nr)

			! No BLAS
			!$omp parallel
			!$omp do
			do i = 1, nr
				res(i) = y(i) - sum(X(i, :)*beta)
			end do
			!$omp end do
			!$omp end parallel

			dssigma = -1*dble(nr)/ssigma + sum(res(:)**2)/ssigma**3

			do i = 1, nc
				dbeta(i) = sum(X(:, i)*2*res)/ssigma**2/2
			end do

			if(.not. usePrecisionMatrix) then
				call g_lg_dmvnorm(nc, beta, mu, Sigma, g_dmvnorm)
			else
				call g_lg_dmvnorm_by_precision(nc, beta, mu, Sigma, g_dmvnorm)
			end if

			g(1) = dssigma + d_lg_igamma(ssigma, ig_a, ig_b)
			g(2:nc + 1) = dbeta + g_dmvnorm
		end subroutine

		double precision function ldnorm(x, sigma)
		implicit none
			double precision, intent(in) :: x, sigma
			ldnorm = -x**2/2/sigma**2 - dlog(dsqrt(2*pi)) - dlog(sigma)
		end function

		double precision function lg_igamma(x, a, b)
			double precision, intent(in) :: x, a, b
			lg_igamma = a*dlog(b) - DLGAMA(a) - b/x - (a + 1)*dlog(x)
		end function

		double precision function d_lg_igamma(x, a, b)
			double precision, intent(in) :: x, a, b
			d_lg_igamma = b / x**2 - (a + 1)/x
		end function

		! L-BFGS-Bの計算用
		subroutine plimit(this, np, nbd, l, u)
			implicit none
			class(blm), intent(inout) :: this
            integer, intent(in) :: np
            integer, dimension(np), intent(out) :: nbd
			double precision, dimension(np), intent(out) :: l, u
			nbd(1) = 1
			l(1) = 0d0
			nbd(2:np) = 0
		end subroutine

		! 事前確率の周辺分布を求める
		double precision function lg_marginal_prior(this, np, nhp, hp, ncnst, pcnst, cnst)
			implicit none
			class(blm), intent(inout) :: this
			integer, intent(in) :: np, nhp, ncnst
			double precision, dimension(nhp), intent(in) :: hp
			integer, dimension(ncnst), intent(in) :: pcnst ! 制約条件の位置
			double precision, dimension(ncnst), intent(in) :: cnst ! 制約
			double precision :: ig_a, ig_b, log_cpd
			double precision, dimension(np - 1) :: mu
			double precision, dimension(np - 1, np - 1) :: Sigma
			integer :: INFO, nc, len_mu
			double precision, allocatable, dimension(:) :: m_mu
			double precision, allocatable, dimension(:, :) :: m_Sigma
			integer, allocatable, dimension(:) :: mvnorm_pcnst
			nc = np - 1

			! ハイパーパラメーターを展開
			ig_a = hp(1)
			ig_b = hp(2)
			mu = hp(3:2 + nc)
			Sigma = reshape(hp(3 + nc:2 + nc + nc**2), (/nc, nc/))
			if(usePrecisionMatrix) then
				! precision matrixをvcovに戻す
				call DPOTRF('U', nc, Sigma, nc, INFO) ! コレスキー分解
				call DPOTRI('U', nc, Sigma, nc, INFO) ! 逆行列
			end if

			! 条件付き事前確率の計算
			! 分散に制約がかかっているかどうかで処理を変える
			if(1 == pcnst(1)) then
				log_cpd = lg_igamma(cnst(1), ig_a, ig_b)
				len_mu = ncnst - 1
				allocate(mvnorm_pcnst(len_mu))
				mvnorm_pcnst = pcnst(2:ncnst) - 1
			else 
				log_cpd = 0
				len_mu = ncnst
				allocate(mvnorm_pcnst(len_mu))
				mvnorm_pcnst = pcnst - 1
			end if
			if(0<len_mu) then
				allocate(m_mu(len_mu), m_Sigma(len_mu, len_mu))
				call calc_marginal_normd_param(nc, mu, Sigma, len_mu, mvnorm_pcnst, m_mu, m_Sigma)
				log_cpd = log_cpd + lg_dmvnorm(len_mu, cnst, m_mu, m_Sigma)
			end if

			! write(*, *) "log_cpd", log_cpd
			! write(*, *) "len_mu", len_mu
			! write(*, *) "m_mu", m_mu
			! write(*, *) "m_Sigma", m_Sigma
			lg_marginal_prior = log_cpd
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
			class(blm), intent(inout) :: this
			integer, intent(in) :: nr, nc, np, ss
			double precision, dimension(nr, nc), intent(inout) :: X
			double precision, dimension(ss, np), intent(in) :: P
			double precision, dimension(nr * ss), intent(out) :: y
			integer :: i, j

			!$omp parallel
			!$omp do private(i)
			do j = 1, ss
				do i = 1, nr
					y(i + nr*(j - 1)) = sum(X(i, :) * P(j, 2:(nc + 1)))
				end do
			end do 
			!$omp end do
			!$omp end parallel

		end subroutine

end module
