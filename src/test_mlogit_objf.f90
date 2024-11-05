subroutine test_mlogit_objf(nr, nc, X, y, np, p, nhp, hp, noc, nok, r)
    use bayesian_mlogit
    implicit none
    integer, intent(in) :: nr, nc, np, nhp, noc, nok
    double precision, dimension(nr, nc), intent(inout) :: X
    double precision, dimension(nr), intent(in) :: y
    double precision, dimension(np), intent(in) :: p
    double precision, dimension(nhp), intent(in) :: hp
    double precision, intent(out) :: r
    type(mlogit) :: this

    this%noc = noc
    this%nok = nok

    r = this%objf(nr, nc, X, y, np, p, nhp, hp)

end subroutine


subroutine test_mlogit_objfg(nr, nc, X, y, np, p, nhp, hp, noc, nok, g)
    use bayesian_mlogit
    use hmc
    implicit none
    integer, intent(in) :: nr, nc, np, nhp, noc, nok
    double precision, dimension(nr, nc), intent(inout) :: X
    double precision, dimension(nr), intent(in) :: y
    double precision, dimension(np), intent(in) :: p
    double precision, dimension(nhp), intent(in) :: hp
    double precision, dimension(np), intent(out) :: g
    type(mlogit) :: this

    this%noc = noc
    this%nok = nok

    call this%objfg(nr, nc, X, y, np, p, nhp, hp, g)

end subroutine
