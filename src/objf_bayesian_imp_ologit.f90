module bayesian_imp_ologit
	use mvnorm
	use bayesian_ologit
    use imputation
	implicit none
	type, extends(ologit) :: ologit_imp
	type(impinfo) :: imp
	integer :: nev, nmu ! 回帰に用いる説明変数の数, 補定に用いる変数の数
	contains
		procedure :: objf => objf_imp
		procedure :: objfg => objfg_imp
		procedure :: lg_marginal_prior => lg_marginal_prior_imp
		procedure :: predict => predict_imp ! 予測値
		procedure :: initialize => initialize_imp
	end type

	contains
	double precision function objf_imp(this, nr, nc, X, y, np, p, nhp, hp) result(r)
		implicit none
		class(ologit_imp), intent(inout) :: this
		integer, intent(in) :: nr, nc, np, nhp
		double precision, dimension(nr, nc), intent(inout) :: X
		double precision, dimension(nr), intent(in) :: y
		double precision, dimension(np), intent(in) :: p
		double precision, dimension(nhp), intent(in) :: hp
		double precision, allocatable, dimension(:) :: mu, mu_X
		double precision, allocatable, dimension(:, :) :: iSigma_X
		double precision :: llf_x
		integer :: info, not, ptr_mu, ptr_mu_X, ptr_iSigma_X, i

		! 補定用のパラメーター
		allocate(mu(this%nmu), mu_X(this%nmu), iSigma_X(this%nmu, this%nmu))
		not = this%nok_ans - 1 ! Number Of Threshold
		ptr_mu = this%nev + not + 1
		mu = p(ptr_mu:ptr_mu + this%nmu - 1)
		! 補定用パラメーターのハイパーパラメーター
		ptr_mu_X = this%nev + this%nev**2 + not + not**2 + 1
		mu_X = hp(ptr_mu_X:ptr_mu_X + this%nmu - 1)
		! write(0, *) "mu_X", mu_X
		ptr_iSigma_X = ptr_mu_X + this%nmu
		iSigma_X = reshape(hp(ptr_iSigma_X:ptr_iSigma_X +  this%nmu**2 - 1), &
			(/this%nmu, this%nmu/))
		! 補定をする/補定に使ったパラメーターの対数尤度を得る
		call this%imp%impute(nr, this%nmu, X(:, 1:this%nmu), mu, llf_x, r, info)
		! 順序ロジットの目的関数を呼ぶ
		r = r + this%ologit%objf(nr, this%nev, X, y, this%nev + not, p, this%nev + this%nev**2 + not + not**2, hp)
		! 補定用パラメーターの対数化事前確率を足す
		r = r + lg_dmvnorm_by_precision(this%nmu, mu, mu_X, iSigma_X)

	end function

	subroutine objfg_imp(this, nr, nc, X, y, np, p, nhp, hp, g)
		implicit none
		class(ologit_imp), intent(inout) :: this
		integer, intent(in) :: nr, nc, np, nhp
		double precision, dimension(nr, nc), intent(inout) :: X
		double precision, dimension(nr), intent(in) :: y
		double precision, dimension(np), intent(in) :: p
		double precision, dimension(nhp), intent(in) :: hp
		double precision, dimension(np), intent(out) :: g

		call this%gradient(nr, nc, X, y, np, p, nhp, hp, g)
	end subroutine

	! 事前確率の周辺分布を求める
	double precision function lg_marginal_prior_imp(this, np, nhp, hp, ncnst, pcnst, cnst) result(r)
		implicit none
		class(ologit_imp), intent(inout) :: this
		integer, intent(in) :: np, nhp, ncnst
		double precision, dimension(nhp), intent(in) :: hp
		integer, dimension(ncnst), intent(in) :: pcnst ! 制約条件の位置
		double precision, dimension(ncnst), intent(in) :: cnst ! 制約
		integer np_anl, nhp_anl
		integer :: ncnst_anl, ncnst_imp ! 線型回帰と補定部分の制約条件の数
		integer, dimension(ncnst):: pcnst_sorted, pcnst_anl, pcnst_imp ! 〃制約条件の位置
		double precision, dimension(ncnst) :: cnst_sorted, cnst_anl, cnst_imp ! 〃制約

		np_anl = this%nev + this%nok_ans - 1
		nhp_anl = np_anl + (this%nev)**2 + (this%nok_ans - 1)**2
		call cnst_sort(ncnst, pcnst, cnst, pcnst_sorted, cnst_sorted)
		call cnst_split(ncnst, np_anl, pcnst_sorted, cnst_sorted, &
			ncnst_anl, pcnst_anl, cnst_anl, &
			ncnst_imp, pcnst_imp, cnst_imp)
		r =  this%ologit%lg_marginal_prior( & 
				np_anl, nhp_anl, hp(1:nhp_anl), ncnst_anl, pcnst_anl, cnst_anl & 
			) + lg_marginal_prior_lp( &
			this%nmu, this%nmu + this%nmu**2, hp(nhp_anl + 1:nhp), &
			ncnst_imp, pcnst_imp, cnst_imp)

	end function

	subroutine predict_imp(this, nr, nc, X, ss, np, P, Y)
		implicit none
		class(ologit_imp), intent(inout) :: this
		integer, intent(in) :: nr, nc, np, ss
		double precision, dimension(nr, nc), intent(inout) :: X
		double precision, dimension(ss, np), intent(in) :: P
		double precision, dimension(nr * ss, this%nok_ans - 1), intent(out) :: Y
		double precision, allocatable, dimension(:) :: mu
		double precision :: llf_x, llf_p
		integer i, info

		! パラメーターの行を一行ごとに処理していく
		! 補定処理を行う
		do i = 1, ss
			mu = P(i, 1 + this%nev:this%nev + this%nmu)
			! 補定をする/補定に使ったパラメーターの対数尤度を得る
			call this%imp%impute(nr, this%nmu, X(:, 1:this%nmu), mu, llf_x, llf_p, info)
			! 線形回帰モデルの予測値をつくるコードを呼ぶ
			call this%ologit%predict(nr, this%nev, X, 1, this%nev, P(i, :), Y(1 + (i - 1)*nr : i*nr, :))
		end do

	end subroutine

	subroutine initialize_imp(this, nr, nc, nev, X, nok)
		implicit none
		class(ologit_imp), intent(inout) :: this
		integer, intent(in) :: nr, nc, nev, nok
		double precision, dimension(nr, nc), intent(inout) :: X
	
		this%h = 1d-5
		this%nev = nev ! 推定に用いる変数の数
		this%nok_ans = nok ! 披説明変数/応答の種類
		this%nmu = nc ! 補定に用いる変数の数
		call this%imp%initialize(nr, this%nmu, X)
	end subroutine

end module
