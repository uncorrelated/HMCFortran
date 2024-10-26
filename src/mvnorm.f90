module mvnorm
	implicit none
	double precision, parameter :: pi = 3.141592653589793115998

	contains
		double precision function lg_dmvnorm(N, y, mu, Sigma)
			integer, intent(in) :: N
			double precision, dimension(N, N), intent(in) :: Sigma
			double precision, dimension(N), intent(in) :: mu, y
			double precision, dimension(N, N) :: Chol
			double precision :: slev
			double precision, dimension(N) :: b
			integer :: INFO, i

			! コレスキー分解
			Chol = Sigma
			call DPOTRF('U', N, Chol, N, INFO)

			! 三角行列Cの行列式の対数を計算
			slev = 0
			do i = 1, N
				slev = slev + DLOG(Chol(i, i))
			end do

			b = y - mu
			call DTRTRS('U', 'T', 'N', N, 1, Chol, N, b, N, INFO)
			lg_dmvnorm = -slev - 0.5 * N * log(2 * pi) - 0.5 * sum(b**2)
		end function 

		subroutine g_lg_dmvnorm(N, y, mu, Sigma, b)
			integer, intent(in) :: N
			double precision, dimension(N, N), intent(in) :: Sigma
			double precision, dimension(N), intent(in) :: mu, y
			double precision, dimension(N), intent(inout) :: b
			double precision, dimension(N, N) :: Chol
			integer :: INFO

			! コレスキー分解
			Chol = Sigma
			call DPOTRF('U', N, Chol, N, INFO)

			b = y - mu
			call DPOTRS('U', N, 1, Chol, N, b, N, INFO)
			b = -1 * b
		end subroutine 

		double precision function lg_dmvnorm_by_precision(N, y, mu, PM)
			integer, intent(in) :: N
			double precision, dimension(N, N), intent(in) :: PM
			double precision, dimension(N), intent(in) :: mu, y
			double precision, dimension(N) :: b

			b = y - mu
			! PM⁻¹の行列式の対数など、サンプリング時に分子と分母でキャンセルアウトされる項は省略
			lg_dmvnorm_by_precision = - 0.5 * SUM(b * MATMUL(PM, b))
		end function 

		subroutine g_lg_dmvnorm_by_precision(N, y, mu, PM, b)
			integer, intent(in) :: N
			double precision, dimension(N, N), intent(in) :: PM
			double precision, dimension(N), intent(in) :: mu, y
			double precision, dimension(N), intent(out) :: b

			b = -1 * MATMUL(PM, y - mu)
		end subroutine 

		subroutine calc_mvnorm_prameters(ns, np, sample, mu, S)
			integer, intent(in) :: ns, np
			double precision, dimension(ns, np), intent(in) :: sample
			double precision, dimension(np), intent(out) :: mu
			double precision, dimension(np, np), intent(out) :: S
			integer :: i, j
	
			do j = 1, np
				mu(j) = sum(sample(:, j)) / ns
			end do
			do j = 1, np
				do i = j, np
					S(i, j) = sum((sample(:, j) - mu(j))*(sample(:, i) - mu(i)))
					S(j, i) = S(i, j)
				end do
			end do
			S = S / ns / (ns - 1)
		end subroutine

		! np: パラメーターの数
		! na: サブセットの要素数
		! a: サブセットのパラメーターの中の位置の配列
		! nb: 補集合の要素数
		! b: 補集合のパラメーターの中の位置の配列（戻り値）
		subroutine setdiff(np, na, a, nb, b)
			implicit none
			integer, intent(in) :: np, na, nb
			integer, dimension(na), intent(in) :: a
			integer, dimension(nb), intent(out) :: b
			integer, dimension(np) :: t
			integer i, j

			do i = 1, np
				t(i) = i
			end do

			do i = 1, na
				t(a(i)) = 0
			end do

			j = 1
			do i = 1, np
				if(.not. t(i) == 0) then
					b(j) = t(i)
					j = j + 1
				end if
			end do

		end subroutine

		! np: pの要素数
		! p: パラメーターの一次元配列
		! n_p_sub: パラメーターのサブセットの要素数
		! i_p_sub: pのどの位置のパラメーターでサブセットをつくるか示す配列
		! p_sub: パラメーターのサブセット（戻り値）
		subroutine subset(np, p, n_p_sub, i_p_sub, p_sub)
			implicit none
			integer, intent(in) :: np, n_p_sub
			double precision, dimension(np), intent(in) :: p
			integer, dimension(n_p_sub), intent(in) :: i_p_sub
			double precision, dimension(n_p_sub), intent(out) :: p_sub
			integer i, j

			do i = 1, n_p_sub
				p_sub(i) = p(i_p_sub(i))
			end do
		end subroutine

		! 正規分布の周辺分布のパラメーターを計算と言うか選ぶ
		! np: 元の分布のパラメーターの数
		! mu: 元の分布の期待値
		! S: 元の分布の分散共分散行列
		! nmp: 周辺分布のパラメーターの数
		! mp: 周辺分布で使うパラメーターの位置
		! m_mu: 周辺分布の期待値（戻り値）
		! m_S: 周辺分布の分散共分散行列（戻り値）
		subroutine calc_marginal_normd_param(np, mu, S, nmp, mp, m_mu, m_S)
			implicit none
			integer, intent(in) :: np, nmp
			double precision, dimension(np), intent(in) :: mu
			double precision, dimension(np, np), intent(in) :: S
			integer, dimension(nmp), intent(in) :: mp
			double precision, dimension(nmp), intent(out) :: m_mu
			double precision, dimension(nmp, nmp), intent(out) :: m_S
			integer i, j

			do j = 1, nmp
				m_mu(j) = mu(mp(j))
				do i = 1, nmp
						m_S(i, j) = S(mp(i), mp(j))
						! write(*, *) i, j, "<-", mp(i), mp(j), S(mp(i), mp(j))
				end do
			end do
		end subroutine

		subroutine calc_gammad_prameters(ns, sample, ig_a_sample, ig_b_sample)
			implicit none
			integer, intent(in) :: ns
			double precision, dimension(ns), intent(in) :: sample
			double precision, intent(out) :: ig_a_sample, ig_b_sample
			integer :: i
			double precision :: mean_sample, var_sample

			! MCMCのサンプルから、ガンマ分布のパラメーターを求める
			mean_sample = sum(sample) / ns
			var_sample = max(1.0d0, sum((sample - mean_sample)**2) / ns)
			ig_a_sample = mean_sample**2 / var_sample + 2
			ig_b_sample = mean_sample * (ig_a_sample - 1)
		end subroutine
end module
