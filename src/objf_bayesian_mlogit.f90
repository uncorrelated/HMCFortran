module bayesian_mlogit
	use mvnorm
    use hmc
	implicit none
	type, extends(model_bf) :: mlogit
    integer :: nok, noc
	contains
        procedure :: objf ! 目的関数
		procedure :: objfg ! objfのグラディエントを計算する関数
		procedure :: plimit ! パラメーターの正負の符号
		procedure :: lg_marginal_prior ! 対数化事前確率
		procedure :: predict ! 予測値
    end type
    contains
        ! 継承の都合でMの代わりにパラメーター名Xを使っている
        double precision function objf(this, nr, nc, X, y, np, p, nhp, hp)
            implicit none
            class(mlogit), intent(inout) :: this
            integer, intent(in) :: nr, nc, np, nhp
            double precision, dimension(nr, nc), intent(inout) :: X
            double precision, dimension(nr), intent(in) :: y
            double precision, dimension(np), intent(in) :: p
            double precision, dimension(nhp), intent(in) :: hp
            double precision, dimension(nr, this%nok) :: Za_Xb
            integer :: i, j
            double precision :: log_n, log_d
			double precision, dimension(np) :: p_mu
			double precision, dimension(np, np) :: p_sigma
            double precision :: a

			! ハイパーパラメーターを展開
			p_mu = hp(1:np)
			p_sigma = reshape(hp(1 + np:np + np**2), (/np, np/))

            ! # individualな部分とconditionな部分の線形予測を合計
            ! Za_Xb <- Xb + Za
            call calc_Za_Xb(this, nr, nc, X, np, p, Za_Xb)

            ! # 尤度関数の分子部分（各行で応答の値に応じた列の値を足している）
            ! log_n <- sum(Za_Xb[1:len + (y-1)*len])
            ! # 尤度関数の分母部分
            ! log_d <- sum(apply(Za_Xb, 1, \(x) log(sum(exp(x)))))
            log_n = 0
            log_d = 0
            do i = 1, nr
                j = int(y(i))
                log_n = log_n + Za_Xb(i, j)
                log_d = log_d + dlog(sum(exp(Za_Xb(i, :))))
            end do

            a = lg_dmvnorm_by_precision(np, p, p_mu, p_sigma)
            if(.not. a==a) then
                write(*, *) "lg_dmvnorm_by_precision(np, p, p_mu, p_sigma)", a
                write(*, *) "p", p
                call rexit("bayesian_mlogit")
            end if
            objf = log_n - log_d + lg_dmvnorm_by_precision(np, p, p_mu, p_sigma)

        end function

        subroutine objfg(this, nr, nc, X, y, np, p, nhp, hp, g)
			implicit none
            class(mlogit), intent(inout) :: this
			integer, intent(in) :: nr, nc, np, nhp
			double precision, dimension(nr, nc), intent(inout) :: X
			double precision, dimension(nr), intent(in) :: y
			double precision, dimension(np), intent(in) :: p
			double precision, dimension(nhp), intent(in) :: hp
			double precision, dimension(np), intent(out) :: g
            double precision, dimension(nr, this%nok) :: Za_Xb
            double precision, dimension(nr) :: denominator
            integer :: i, j, k, s, e, ps, pe, ncol_X
			double precision, dimension(np) :: p_mu
			double precision, dimension(np, np) :: p_sigma
			double precision, dimension(np) :: p_g_dmvnorm

			! ハイパーパラメーターを展開
			p_mu = hp(1:np)
			p_sigma = reshape(hp(1 + np:np + np**2), (/np, np/))

            call calc_Za_Xb(this, nr, nc, X, np, p, Za_Xb)

            ! # (log(A))' = 1/Aなので対数が消える
            ! denominator <- apply(Za_Xb, 1, \(x) sum(exp(x)))
            do i = 1, nr
                denominator(i) = sum(exp(Za_Xb(i, :)))
            end do

           ! conditonalな変数の係数の微分
            if(0 < this%noc) then
                do j = 1, this%noc
                    s = 1 + this%nok*(j - 1)
                    e = this%nok*j
                    g(j) = 0
                    do i = 1, nr
                        g(j) = g(j) + X(i, s - 1 + int(y(i))) - sum(X(i, s:e)*exp(Za_Xb(i, :))/denominator(i))
                    end do
                end do
            end if

            ! individualな変数の係数の微分
            if(this%noc < np) then 
                ! Mの中のXに該当する列の最初と最後
                s = 1 + this%nok * this%noc
                e = nc
                ncol_X = e - s + 1
                ! pの中のβに該当する部分の最初と最後
                ps = 1 + this%noc
                pe = np
                ! beta(i, j) = p(i + (j-2)*this%nok)
                ! j == 1のパラメーターは計算しないので省略
                do j = 2, this%nok
                    ! nrow(beta) == ncol_X
                    do i = 1, ncol_X
                        g(ps - 1 + i + (j - 2) * ncol_X) = &
                            + sumif(nr, X(:, s + i - 1), y(:), dble(j)) &
                            - sum(X(:, s + i - 1)*exp(Za_Xb(:, j))/denominator)
                    end do
                end do
            end if

			call g_lg_dmvnorm_by_precision(np, p, p_mu, p_sigma, p_g_dmvnorm)
            g = g + p_g_dmvnorm

        end subroutine

        ! x,a同じ長さのベクトル、bはスカラーで、a[i]==bのときだけ合計する
        double precision function sumif(n, x, a, b)
            integer, intent(in) :: n
            double precision, dimension(n), intent(in) :: x, a
            double precision, intent(in) :: b
            double precision :: s
            integer :: i
            ! first_term <- sum(X[, i]*(1 - pmin(abs(y - j), 1)))
            s = 0
            do i = 1, n
                s = s + x(i)*(1d0 - min(abs(a(i) - b), 1d0))
            end do
            sumif = s
        end function

        subroutine calc_Za_Xb(this, nr, nc, M, np, p, Za_Xb)
            implicit none
            type(mlogit) :: this
            integer, intent(in) :: nr, nc, np
            double precision, dimension(nr, nc), intent(in) :: M
            double precision, dimension(np), intent(in) :: p
            double precision, dimension(nr, this%nok), intent(out) :: Za_Xb
            double precision, dimension(nr, this%nok) :: Xb, ZaSum
            double precision, allocatable, dimension(:, :) :: beta ! , Za
            double precision, allocatable, dimension(:) :: alpha
            integer :: ncol_Z, ncol_X, cs_Z = 0, ce_Z = 0, cs_X = 0, ce_X = 0
            integer :: i, j, k

            ! noc: conditionalな説明変数の数 → 〃 な係数の数
            ! nok: 選択肢の数
            ! Zα1 + Zα2 + … + Xβを計算する
            ! β: nok列☓ncol(X)行 
            ! Xβ: nok列
            ! Zα1は内積ではなくて、列で対応した積で、nok列になる
            ! 引数pの前半はαのベクター後半はβのベクター
    
            ncol_Z = this%nok * this%noc
            ncol_X = nc - ncol_Z

            allocate(alpha(this%noc))
            allocate(beta(ncol_X, this%nok))

            ZaSum = 0

            if(0<ncol_Z) then
                ! allocate(Za(nr, this%nok * this%noc))
                do i = 1, this%noc
                    alpha(i) = p(i)
                end do
                cs_Z = 1
                ce_Z = this%nok * this%noc
                do j = 1, this%noc
                    do k = 1, this%nok
                        do i = 1, nr
                            ! Za(i, k + (1-j)*this%nok) = alpha(j) * M(i, k + (1-j)*this%nok)
                            ZaSum(i, k) = ZaSum(i, k) + alpha(j) * M(i, k + (1-j)*this%nok)
                        end do
                    end do
                end do 
            end if

            if(0<ncol_X) then
                do j = 2, this%nok
                    do i = 1, ncol_X
                        beta(i, j) = p(i + this%noc + (j - 2)*ncol_X)
                    end do
                end do
                beta(:, 1) = 0

                cs_X = ce_Z + 1
                ce_X = nc
                do j = 1, this%nok
                    do i = 1, nr
                        Xb(i, j) = sum(M(i, cs_X:ce_X) * beta(:, j))                   
                    end do
                end do 
            else
                Xb = 0
            end if

            Za_Xb = ZaSum + Xb

        end subroutine

        ! L-BFGS-Bの計算用
        subroutine plimit(this, np, nbd, l, u)
            implicit none
            class(mlogit), intent(inout) :: this
            integer, intent(in) :: np
            integer, dimension(np), intent(out) :: nbd
            double precision, dimension(np), intent(out) :: l, u
            nbd = 0
        end subroutine

        ! 事前確率の周辺分布を求める
        double precision function lg_marginal_prior(this, np, nhp, hp, ncnst, pcnst, cnst)
            implicit none
            class(mlogit), intent(inout) :: this
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
        ! Y: 予測値（出力）
        subroutine predict(this, nr, nc, X, ss, np, P, Y)
            implicit none
            class(mlogit), intent(inout) :: this
            integer, intent(in) :: nr, nc, np, ss
            double precision, dimension(nr, nc), intent(inout) :: X
            double precision, dimension(ss, np), intent(in) :: P
            double precision, dimension(nr*ss, this%nok), intent(out) :: Y
            integer :: i, j
            double precision, dimension(nr, this%nok) :: Za_Xb
            double precision, dimension(this%nok) :: exp_Xb
    
            !$omp parallel
            !$omp do private(i, Za_Xb, exp_Xb)
            do j = 1, ss
                call calc_Za_Xb(this, nr, nc, X, np, P(j, :), Za_Xb)
                do i = 1, nr
                    exp_Xb = exp(Za_Xb(i, :))
                    Y(i + nr*(j - 1), :) = exp_Xb/sum(exp_Xb)
                end do
                call rchkusr()
            end do
            !$omp end do
            !$omp end parallel
    
        end subroutine

end module
