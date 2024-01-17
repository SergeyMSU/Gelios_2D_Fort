module Phys_parameter
    use STORAGE 
    use GEOMETRY
    implicit none 

    contains

    subroutine Inner_Conditions(SS, x, y, par)
        !! Задаём внутренние граничные условия
        TYPE (Setka), intent(in) :: SS
        real(8), intent(in) :: x, y
        real(8), intent(out) :: par(:)
        real(8) :: r, p_0

        r = sqrt(x**2 + y**2)

        par(1) = 6.0 * (SS%par_R0/r)**2
        par(3) = x/r * 430.0
        par(4) = y/r * 430.0
        p_0 = 430.0**2 * 6.0**2/(SS%par_ggg * 10.0**2)
        par(2) = p_0 * (SS%par_R0/r)**(2.0 * SS%par_ggg)
        par(5) = 1.0

    end subroutine Inner_Conditions

    subroutine Phys_Innitial_Conditions(SS, x, y, par)
        !! Задаём начальные условия
        TYPE (Setka), intent(in) :: SS
        real(8), intent(in) :: x, y
        real(8), intent(out) :: par(:)
        real(8) :: r, p_0

        r = sqrt(x**2 + y**2)

        if(r < 50.0) then
            par(1) = 6.0 * (SS%par_R0/r)**2
            par(3) = x/r * 430.0
            par(4) = y/r * 430.0
            p_0 = 430.0**2 * 6.0**2/(SS%par_ggg * 10.0**2)
            par(2) = p_0 * (SS%par_R0/r)**(2.0 * SS%par_ggg)
            par(5) = 1.0
        else
            par(1) = 1.0
            par(3) = SS%par_Velosity_inf
            par(4) = 0.0
            par(2) = 1.0
            par(5) = 100.0
        end if
    end subroutine Phys_Innitial_Conditions

    subroutine Phys_input_flow(SS, par)
        !! Задаём набегающий поток
        TYPE (Setka), intent(in) :: SS
        real(8), intent(out) :: par(:)

        par(1) = 1.0
        par(3) = SS%par_Velosity_inf
        par(4) = 0.0
        par(2) = 1.0
        par(5) = 100.0

    end subroutine Phys_input_flow

end module Phys_parameter