module bayesian_imp_lm
    use mvnorm
    use bayesian_lm
    use imputation
    type, extends(blm) :: blm_imp
        type(impinfo) :: imp
        double precision, pointer, dimension(:, :) :: X ! 補定される用
        integer :: nev, nmu ! 回帰に用いる説明変数の数, 補定に用いる変数の数
        contains
            procedure :: objf => objf_imp
            procedure :: objfg => objfg_imp
            procedure :: plimit => plimit_imp
            procedure :: lg_marginal_prior => lg_marginal_prior_imp
            procedure :: predict => predict_imp
    end type
    contains
        double precision function objf_imp(this, nr, nc, X, y, np, p, nhp, hp) result(r)
            implicit none
            class(blm_imp), intent(inout) :: this
            integer, intent(in) :: nr, nc, np, nhp
            double precision, dimension(nr, nc), intent(in) :: X
            double precision, dimension(nr), intent(in) :: y
            double precision, dimension(np), intent(in) :: p
            double precision, dimension(2 + (nc + nc**2) + (nc - 1 + (nc - 1)**2)), intent(in) :: hp
            double precision, allocatable, dimension(:) :: mu, mu_X
            double precision, allocatable, dimension(:, :) :: iSigma_X
            double precision :: llf_x
            integer :: info, nmu

            ! 補定用のパラメーター
            nmu = (this%nev - 1) + (nc - this%nev)
            allocate(mu(nmu), mu_X(nmu), iSigma_X(nmu, nmu))
            mu = p(2 + this%nev:1 + this%nev + nmu)
            ! 補定用パラメーターのハイパーパラメーター
            mu_X = hp(3 + this%nev + this%nev**2:(2 + this%nev + this%nev**2) + nmu)
            iSigma_X = reshape(hp(1 + (2 + this%nev + this%nev**2) + nmu:(2 + this%nev + this%nev**2) + nmu +  nmu**2), &
                (/nmu, nmu/))
            ! 補定をする/補定に使ったパラメーターの対数尤度を得る
            call this%imp%impute(nr, nmu, this%X(:, 2:nmu), mu, llf_x, r, info)
            ! 線形回帰の目的関数を呼ぶ
            r = r + this%blm%objf(nr, this%nev, this%X, y, 1 + this%nev, p, 2 + (this%nev + this%nev**2), hp)
            ! 補定用パラメーターの対数化事前確率を足す
            r = r + lg_dmvnorm_by_precision(nmu, mu, mu_X, iSigma_X)
        end function

        subroutine objfg_imp(this, nr, nc, X, y, np, p, nhp, hp, g)
			implicit none
			class(blm_imp), intent(inout) :: this
			integer, intent(in) :: nr, nc, np, nhp
			double precision, dimension(nr, nc), intent(in) :: X
			double precision, dimension(nr), intent(in) :: y
            double precision, dimension(np), intent(in) :: p
            double precision, dimension(nhp), intent(in) :: hp
			double precision, dimension(nc + 1 + nc - 1), intent(out) :: g
			double precision, dimension(np) :: q
			integer :: i, j
			double precision :: h
			double precision, dimension(4) :: v

            ! 数値微分を行う/リチャードソンの外挿
			h = this%h
			q(:) = p(:) ! 入力データはいじらない
			do i = 1, np
                do j = 1, 4
                    q(i) = p(i) + 2*h - (j-1)*h
                    v(j) = this%objf(nr, nc, X, y, np, q, nhp, hp)
                end do
				g(i) = (- v(1) + 8*v(2) - 8*v(3) + v(4))/(12*h)
				q(i) = p(i) ! 戻す
			end do
        end subroutine

        ! L-BFGS-Bの計算用
		subroutine plimit_imp(this, np, nbd, l, u)
			implicit none
			class(blm_imp), intent(inout) :: this
            integer, intent(in) :: np
            integer, dimension(np), intent(out) :: nbd
			double precision, dimension(np), intent(out) :: l, u
			nbd(1) = 1
			l(1) = 0d0
			nbd(2:np) = 0
		end subroutine

		double precision function lg_marginal_prior_imp(this, np, nhp, hp, ncnst, pcnst, cnst) result(r)
            implicit none
			class(blm_imp), intent(inout) :: this
			integer, intent(in) :: np, nhp, ncnst
			double precision, dimension(nhp), intent(in) :: hp
			integer, dimension(ncnst), intent(in) :: pcnst ! 制約条件の位置
			double precision, dimension(ncnst), intent(in) :: cnst ! 制約
            integer np_anl
            integer :: ncnst_anl, ncnst_imp ! 線型回帰と補定部分の制約条件の数
            integer, dimension(ncnst):: pcnst_sorted, pcnst_anl, pcnst_imp ! 〃制約条件の位置
            double precision, dimension(ncnst) :: cnst_sorted, cnst_anl, cnst_imp ! 〃制約

            np_anl = this%nev + 1
            call cnst_sort(ncnst, pcnst, cnst, pcnst_sorted, cnst_sorted)
            call cnst_split(ncnst, np_anl, pcnst_sorted, cnst_sorted, &
                ncnst_anl, pcnst_anl, cnst_anl, &
                ncnst_imp, pcnst_imp, cnst_imp)
            r =  this%blm%lg_marginal_prior( & 
                    np_anl, 2 + this%nev + this%nev**2, hp(1:2 + this%nev + this%nev**2), ncnst_anl, pcnst_anl, cnst_anl & 
                ) + lg_marginal_prior_lp( &
                this%nmu, this%nmu + this%nmu**2, hp(3 + this%nev + this%nev**2:nhp), &
                ncnst_imp, pcnst_imp, cnst_imp)

        end function

		subroutine predict_imp(this, nr, nc, X, ss, np, P, y)
			implicit none
			class(blm_imp), intent(inout) :: this
			integer, intent(in) :: nr, nc, np, ss
			double precision, dimension(nr, nc), intent(in) :: X
			double precision, dimension(ss, np), intent(in) :: P
            double precision, dimension(nr * ss), intent(out) :: y
            double precision, allocatable, dimension(:) :: mu
            double precision :: llf_x, llf_p
            integer i, nmu, info

            ! パラメーターの行を一行ごとに処理していく
            ! 補定処理を行う
            nmu = (this%nev - 1) + (nc - this%nev)
            do i = 1, ss
                mu = P(i, 2 + this%nev:1 + this%nev + nmu)
                ! 補定をする/補定に使ったパラメーターの対数尤度を得る
                call this%imp%impute(nr, nmu, this%X(:, 2:nmu), mu, llf_x, llf_p, info)
                ! 線形回帰モデルの予測値をつくるコードを呼ぶ
                call this%blm%predict(nr, this%nev, X, 1, 1 + this%nev, P, y(1 + (i - 1)*nr : i*nr))
            end do

        end subroutine

end module
