subroutine test_ologit_objf(nr, nc, X, y, np, p, nhp, hp, r)
    use bayesian_ologit
    use hmc
    implicit none
    integer, intent(in) :: nr, nc, np, nhp
    double precision, dimension(nr, nc), intent(in) :: X
    double precision, dimension(nr), intent(in) :: y
    double precision, dimension(np), intent(in) :: p
    double precision, dimension(nhp), intent(in) :: hp
    double precision, intent(out) :: r
    type(ologit) :: this

    r = this%objf(nr, nc, X, y, np, p, nhp, hp)

end subroutine


subroutine test_ologit_objfg(nr, nc, X, y, np, p, nhp, hp, g)
    use bayesian_ologit
    use hmc
    implicit none
    integer, intent(in) :: nr, nc, np, nhp
    double precision, dimension(nr, nc), intent(in) :: X
    double precision, dimension(nr), intent(in) :: y
    double precision, dimension(np), intent(in) :: p
    double precision, dimension(nhp), intent(in) :: hp
    double precision, dimension(np), intent(out) :: g
    type(ologit) :: this

    call this%objfg(nr, nc, X, y, np, p, nhp, hp, g)

end subroutine
