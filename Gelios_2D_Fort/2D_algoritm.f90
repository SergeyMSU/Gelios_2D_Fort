module Algoritm
    USE STORAGE
    USE GEOMETRY
    USE Phys_parameter
    implicit none 

    contains

    subroutine Gas_dynamic_algoritm(SS, par)
        TYPE (Setka), intent(in out) :: SS
        TYPE (Phys_par), intent(in out) :: par

        if(par%init == .False.) then
            STOP "Gas_dynamic_algoritm  error init jegtetf94974"  
        end if

        if(SS%init_geo == .False.) then
            STOP "Gas_dynamic_algoritm  error init_geo erty67543"  
        end if

        call Bound_condition(SS, par)  ! Зададим граничные условия на внутренней сфере

    end subroutine Gas_dynamic_algoritm


    subroutine Bound_condition(SS, par)
        !! Задаём граничные условия в сетке (на внутренней сфере для газовой динамики)
        TYPE (Setka), intent(in out) :: SS
        TYPE (Phys_par), intent(in out) :: par
        integer(4) :: n1, i, cell
        real(8) :: x, y, rho, p, u, v

        n1 = size(SS%gl_Cell_A(1,:))
        do i = 1, n1
            cell = SS%gl_Cell_A(1,i)
            x = SS%gl_Cell_Centr(1, cell, 1)
            y = SS%gl_Cell_Centr(1, cell, 1)
            call Inner_Conditions(SS, x, y, rho, p, u, v)
            par%gd(1, cell, 1) = rho
            par%gd(2, cell, 1) = p
            par%gd(3, cell, 1) = u
            par%gd(4, cell, 1) = v
        end do

        n1 = size(SS%gl_Cell_B(1,:))
        do i = 1, n1
            cell = SS%gl_Cell_B(1,i)
            x = SS%gl_Cell_Centr(1, cell, 1)
            y = SS%gl_Cell_Centr(1, cell, 1)
            call Inner_Conditions(SS, x, y, rho, p, u, v)
            par%gd(1, cell, 1) = rho
            par%gd(2, cell, 1) = p
            par%gd(3, cell, 1) = u
            par%gd(4, cell, 1) = v
        end do


    end subroutine Bound_condition



end module Algoritm