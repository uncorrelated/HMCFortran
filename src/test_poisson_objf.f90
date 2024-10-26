subroutine test_poisson_objf(nr, nc, X, y, np, p, nhp, hp, r)
    use poisson_exp
    implicit none
    integer, intent(in) :: nr, nc, np, nhp
    double precision, dimension(nr, nc), intent(in) :: X
    double precision, dimension(nr), intent(in) :: y
    double precision, dimension(np), intent(in) :: p
    double precision, dimension(nhp), intent(in) :: hp
    double precision, intent(out) :: r
    type(poisson) :: this

    r = this%objf(nr, nc, X, y, np, p, nhp, hp)

end subroutine

subroutine test_poisson_objfg(nr, nc, X, y, np, p, nhp, hp, g)
    use poisson_exp
    implicit none
    integer, intent(in) :: nr, nc, np, nhp
    double precision, dimension(nr, nc), intent(in) :: X
    double precision, dimension(nr), intent(in) :: y
    double precision, dimension(np), intent(in) :: p
    double precision, dimension(nhp), intent(in) :: hp
    double precision, dimension(np), intent(out) :: g
    type(poisson) :: this

    call this%objfg(nr, nc, X, y, np, p, nhp, hp, g)

end subroutine
