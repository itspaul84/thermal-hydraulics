module parameters

    implicit none

    integer, parameter :: ncells = 10

    real(8), dimension(ncells) :: temp
    real(8), dimension(ncells) :: pres
    real(8), dimension(ncells) :: rhof
    real(8), dimension(ncells) :: rhog
    real(8), dimension(ncells) :: psat


contains

    subroutine comp_dens()
        implicit none
        integer :: i
        real(8) :: T, Tc, Pc, rho_c, tau, R

        Tc = 647.096d0       ! Critical temperature (K)
        Pc = 22.064d6        ! Critical pressure (Pa)
        rho_c = 322.0d0      ! Critical density (kg/m^3)
        R = 461.526d0        ! Specific gas constant for water vapor (J/kg/K)

        do i = 1, ncells
            T = temp(i)
            if (T < 273.15d0 .or. T > Tc) then
                rhof(i) = -1.0d0
                rhog(i) = -1.0d0
            else
                tau = 1.0d0 - T / Tc

                ! --- Saturation pressure using Wagner equation ---
                psat(i) = Pc * exp( &
                    ( -7.85951783d0 * tau + &
                        1.84408259d0 * tau**1.5d0 - &
                    11.7866497d0 * tau**3.0d0 + &
                    22.6807411d0 * tau**3.5d0 - &
                    15.9618719d0 * tau**4.0d0 + &
                        1.80122502d0 * tau**7.5d0 ))

                ! --- Saturated liquid density (empirical fit) ---
                rhof(i) = rho_c * (1.992740d0 * tau**(1.0d0/3.0d0) + &
                                            1.099653d0 * tau**(2.0d0/3.0d0) - &
                                            0.510839d0 * tau + 1.0d0)

                ! --- Saturated vapor density using ideal gas law ---
                rhog(i) = psat(i) / (R * T)

            end if
        end do

    end subroutine comp_dens


end module parameters