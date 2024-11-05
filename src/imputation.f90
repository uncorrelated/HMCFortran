module imputation
    use IEEE_ARITHMETIC
    use mvnorm
    use lbfgsb_module
    use iso_fortran_env, only: output_unit
    implicit none

    type impinfo
        ! n_na: NAの数 n_na_row: NAがある行の数 na_row: NAがある行
        integer :: n_na, n_na_row
        ! XとNAの位置
        integer, pointer, dimension(:, :) :: ptr_na, ptr_na_m
        logical, pointer, dimension(:) :: is_na_row
        integer, pointer, dimension(:) :: ptr_na_row
        contains
            procedure :: initialize
            procedure :: impute
            procedure :: finalize
    end type

    contains
        subroutine initialize(this, nr, nc, X)
            implicit none
            class(impinfo) :: this
            integer, intent(in) :: nr, nc
            double precision, dimension(nr, nc), intent(inout) :: X
            integer :: h, i, j, k

            ! NAの数と行を特定する
            allocate(this%is_na_row(nr))
            this%n_na = 0
            this%is_na_row = .false.
            this%n_na_row = 0
            do i = 1, nr
                do j = 1, nc
                    if(.not.IEEE_IS_FINITE(X(i, j))) then
                        this%is_na_row(i) = .true.
                        this%n_na = this%n_na + 1
                    end if
                end do
                if(this%is_na_row(i)) then
                    this%n_na_row = this%n_na_row + 1
                end if
            end do

            ! NAがある行だけのサブセットをつくる/同時に欠損値の位置を記録
            allocate(this%ptr_na(this%n_na, 2), this%ptr_na_row(this%n_na_row), this%ptr_na_m(this%n_na, 2))
            k = 1
            h = 1
            do i = 1, nr
                if(this%is_na_row(i)) then
                    do j = 1, nc
                        if(.not.IEEE_IS_FINITE(X(i, j))) then
                            this%ptr_na_m(h, 1) = k
                            this%ptr_na_m(h, 2) = j
                            this%ptr_na(h, 1) = i
                            this%ptr_na(h, 2) = j
                            h = h + 1
                        end if
                    end do
                    this%ptr_na_row(k) = i
                    k = k + 1
                end if
            end do
        end subroutine

        subroutine impute(this, nr, nc, X, mu, llf_x, llf_p, info)
            implicit none
            class(impinfo) :: this
            integer, intent(in) :: nr, nc
            double precision, dimension(nr, nc), intent(inout) :: X
            double precision, dimension(nc), intent(in) :: mu
            ! llf_x: 補定された値の対数尤度, llf_p: 補定に使ったパラメーターの対数尤度
            double precision, intent(out) :: llf_x, llf_p
            double precision, dimension(nc, nc) :: Sigma, iSigma
            double precision :: log_sqrt_det
            double precision, allocatable, dimension(:) :: p, g ! p: 補定する値, g: グラディエント 
            integer :: i
            !!! L-BFGS-Bで利用する変数
            integer, parameter :: m = 10, iprint = -1
            double precision, parameter :: factr = 1.0d+7, pgtol = 1.0d-8
            character(len=60) :: task, csave
            logical :: lsave(4)
            integer :: isave(44)
            double precision :: f
            double precision :: dsave(29)
            integer, allocatable :: nbd(:), iwa(:)
            double precision, allocatable :: l(:), u(:), wa(:)
            integer, intent(out) :: info

            ! 分散共分散行列を計算
            call colsCov(Sigma)
            ! 分散共分散行列の逆行列と固有値のルートを求める
            call invSymMatrix(nc, Sigma, iSigma, log_sqrt_det)
            ! 最尤法で補定する値を定める
            allocate(p(this%n_na), g(this%n_na))
            ! 初期値を代入
            p = 1
            ! L-BFGS-Bの計算用の領域を確保
            allocate ( nbd(this%n_na), l(this%n_na), u(this%n_na) )
            allocate ( iwa(3*this%n_na) )
            allocate ( wa(2*m*this%n_na + 5*this%n_na + 11*m*m + 8*m) )
            ! 上限と下限を設定
            nbd = 0
        
            f = -impf(p)
            call impfg(p, g)
            g = -g
        
            task = 'START'
            do while(task(1:2)=='FG'.or.task=='NEW_X'.or.task=='START')
                call setulb ( this%n_na, m, p, l, u, nbd, f, g, factr, pgtol, &
                                wa, iwa, task, iprint,&
                                csave, lsave, isave, dsave, &
                                iteration_file = 'driver1_output.txt' )
                if (task(1:2) == 'FG') then
                    f = -impf(p)
                    call impfg(p, g)
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

            ! 尤度をまとめる
            llf_x = impf(p)
            llf_p = 0
            do i = 1, nr 
                if(.not.this%is_na_row(i)) then
                    llf_p = llf_p + lg_dmvnorm_by_precision(nc, X(i, :), mu, iSigma, log_inv_sqrt_det = -log_sqrt_det) 
                end if
            end do

            ! 補定値を代入する
            call impX(p, X)
        
            contains
                subroutine impX(p, X)
                    implicit none
                    double precision, dimension(nr, nc), intent(inout) :: X
                    double precision, dimension(this%n_na), intent(in) :: p
                    integer i, j, k
        
                    do k = 1, this%n_na
                        i = this%ptr_na(k, 1)
                        j = this%ptr_na(k, 2)
                        X(i, j) = p(k)
                    end do
                end subroutine
        
                double precision function impf(p) result(r)
                    implicit none
                    double precision, dimension(this%n_na), intent(in) :: p
                    integer i
        
                    call impX(p, X)
        
                    r = 0
                    do i = 1, this%n_na_row
                        r = r + lg_dmvnorm_by_precision(nc, X(this%ptr_na_row(i), :), mu, iSigma, log_inv_sqrt_det = -log_sqrt_det)
                    end do
                end function
        
                subroutine impfg(p, g)
                    implicit none
                    double precision, dimension(this%n_na), intent(in) :: p
                    double precision, dimension(this%n_na), intent(out) :: g
                    double precision, dimension(this%n_na_row, nc) :: M
                    integer i, j, k
        
                    call impX(p, X)
        
                    do i = 1, this%n_na_row
                        M(i, :) = -1*MATMUL(iSigma, X(this%ptr_na_row(i), :) - mu)
                    end do
        
                    do k = 1, this%n_na
                        i = this%ptr_na_m(k, 1)
                        j = this%ptr_na_m(k, 2)
                        g(k) = M(i, j)
                    end do
                end subroutine
        
                subroutine colsCov(Sigma)
                    implicit none
                    double precision, dimension(nc, nc), intent(out) :: Sigma
                    double precision, dimension(nr, nc) :: X0
                    integer :: i, j, k
                    double precision :: s
        
                    do i = 1, nr
                        if(.not.this%is_na_row(i)) then
                            X0(i, :) = X(i, :) - mu
                        end if
                    end do
                    do j = 1, nc
                        do k = 1, j
                            s = 0
                            do i = 1, nr
                                if(.not.this%is_na_row(i)) then 
                                    s = s + X0(i, j) * X0(i, k)
                                end if
                            end do
                            Sigma(k, j) = s / (this%n_na_row - 1)
                            Sigma(j, k) = Sigma(k, j)
                        end do
                    end do
                end subroutine
        
                subroutine invSymMatrix(N, M, IM, log_sqrt_det)
                    implicit none
                    integer, intent(in) :: N
                    double precision, dimension(N, N), intent(in) :: M
                    double precision, dimension(N, N), intent(out) :: IM
                    double precision, dimension(N, N) :: C
                    integer :: INFO, i
                    double precision, intent(out) :: log_sqrt_det

                    ! 正定値対称行列のコレスキー分解
                    C = M
                    call DPOTRF('U', N, C, N, INFO)
        			! 三角行列Cの行列式の対数を計算
                    log_sqrt_det = 0
                    do i = 1, N
                        log_sqrt_det = log_sqrt_det + DLOG(C(i, i))
                    end do
                    ! DPOTRF で計算した結果を用いて、Sの逆行列を計算する。
                    IM = C
                    call DPOTRI('U', N, IM, N, INFO)
                    call copyUpperTriToLowerTri(N, IM)
                end subroutine
        
                subroutine copyUpperTriToLowerTri(N, M)
                    implicit none
                    integer, intent(in) :: N
                    double precision, dimension(N, N), intent(inout) :: M
                    integer :: i, j
                    do j = 1, N - 1
                        do i = j + 1, N
                            M(i, j) = M(j, i)
                        end do
                    end do
                end subroutine
       
        end subroutine

        subroutine finalize(this)
            implicit none
            class(impinfo) :: this
            deallocate(this%ptr_na, this%ptr_na_m, this%is_na_row, this%ptr_na_row)
        end subroutine

end module
