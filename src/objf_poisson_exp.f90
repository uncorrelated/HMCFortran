module poisson_exp
    use mvnorm
	use hmc
	implicit none
	type, extends(model_bf) :: poisson
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
            class(poisson), intent(inout) :: this
            integer, intent(in) :: nr, nc, np, nhp
            double precision, dimension(nr, nc), intent(in) :: X
            double precision, dimension(nr), intent(in) :: y
            double precision, dimension(np), intent(in) :: p
            double precision, dimension(nhp), intent(in) :: hp
            integer :: i, j
			double precision, dimension(np) :: p_mu
			double precision, dimension(np, np) :: p_sigma
            double precision, dimension(nr) :: lambda

			! ハイパーパラメーターを展開
			p_mu = hp(1:np)
			p_sigma = reshape(hp(1 + np:np + np**2), (/np, np/))

            lambda = dexp(MATMUL(X, p))

            objf = sum(y*dlog(lambda) - lambda - dlgama(y + 1)) &
                + lg_dmvnorm_by_precision(np, p, p_mu, p_sigma)
        end function

        subroutine objfg(this, nr, nc, X, y, np, p, nhp, hp, g)
			implicit none
            class(poisson), intent(inout) :: this
			integer, intent(in) :: nr, nc, np, nhp
			double precision, dimension(nr, nc), intent(in) :: X
			double precision, dimension(nr), intent(in) :: y
			double precision, dimension(np), intent(in) :: p
			double precision, dimension(nhp), intent(in) :: hp
			double precision, dimension(np), intent(out) :: g
			double precision, dimension(np) :: p_mu
			double precision, dimension(np, np) :: p_sigma
			double precision, dimension(np) :: p_g_dmvnorm
            integer :: i, j
            double precision, dimension(nr) :: lambda
            double precision, dimension(nr, nc) :: dlambda

			! ハイパーパラメーターを展開
			p_mu = hp(1:np)
			p_sigma = reshape(hp(1 + np:np + np**2), (/np, np/))

            ! lambda <- exp(X %*% p)
            lambda = dexp(MATMUL(X, p))
            ! dlambda <- apply(X, 2, \(x) x*lambda)
            ! apply(y*X - dlambda, 2, sum)
            do j = 1, nc
                g(j) = sum(y * X(:, j) - X(:, j)*lambda)
            end do

			call g_lg_dmvnorm_by_precision(np, p, p_mu, p_sigma, p_g_dmvnorm)
            g = g + p_g_dmvnorm
           
        end subroutine
	! L-BFGS-Bの計算用
        subroutine plimit(this, np, nbd, l, u)
            implicit none
            class(poisson), intent(inout) :: this
            integer, intent(in) :: np
            integer, dimension(np), intent(out) :: nbd
            double precision, dimension(np), intent(out) :: l, u
            nbd = 0
        end subroutine

        ! 事前確率の周辺分布を求める
        double precision function lg_marginal_prior(this, np, nhp, hp, ncnst, pcnst, cnst)
            implicit none
            class(poisson), intent(inout) :: this
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
            class(poisson), intent(inout) :: this
            integer, intent(in) :: nr, nc, np, ss
            double precision, dimension(nr, nc), intent(in) :: X
            double precision, dimension(ss, np), intent(in) :: P
            double precision, dimension(nr * ss), intent(out) :: y
            integer :: i, j
    
            !$omp parallel
            !$omp do private(i)
            do j = 1, nr
                do i = 1, ss
                    y(i + nr*(j - 1)) = dexp(sum(X(j, :) * P(i, :)))
                end do
            end do
            !$omp end do
            !$omp end parallel
    
        end subroutine
    end module
