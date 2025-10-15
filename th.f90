program th

    use parameters
    use solver
    use io
    implicit none

    integer :: i

    do i = 1, ncells
        temp(i) = 300.0d0 + 10.0d0 * sin(2.0d0 * 3.14159d0 * i / ncells)
        pres(i) = 1.0e5 + 500.0d0 * cos(2.0d0 * 3.14159d0 * i / ncells)
    end do

    call comp_dens()

    print '(A4, 5A12)', "Cell", "Temp(K)", "Pres(Pa)", "Psat(Pa)", "Rho_liq", "Rho_vap"
    do i = 1, ncells
        print '(I4, 5F12.2)', i, temp(i), pres(i), psat(i), rhof(i), rhog(i)
    end do


end program th
