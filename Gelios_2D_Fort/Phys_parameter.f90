module Phys_parameter
    use STORAGE 
    use GEOMETRY
    USE ieee_arithmetic
    implicit none 

    contains

    subroutine Inner_Conditions(SS, x, y, par)
        !! Задаём внутренние граничные условия
        TYPE (Setka), intent(in) :: SS
        real(8), intent(in) :: x, y
        real(8), intent(out) :: par(:)
        real(8) :: r, p_0

        r = sqrt(x**2 + y**2)

        if(size(par) < 5) then
            print*, "Error 18 Inner_Conditions size(par) nm5453fdbbdfcsd"
            print*, size(par)
            STOP
        end if

        ! par(1) = 1.0
        ! par(3) = SS%par_Velosity_inf
        ! par(4) = 0.0
        ! par(2) = 1.0
        ! par(5) = 100.0
        ! return

        par(1) = 150.0 * (SS%par_R0/r)**2
        par(3) = x/r * 41.0
        par(4) = y/r * 41.0
        p_0 = 41.0**2 * 150.0/(SS%par_ggg * 7.0**2)
        par(2) = p_0 * (SS%par_R0/r)**(2.0 * SS%par_ggg)
        par(5) = 150.0 * (SS%par_R0/r)**2

    end subroutine Inner_Conditions

    subroutine Phys_Innitial_Conditions(SS, x, y, par)
        !! Задаём начальные условия
        TYPE (Setka), intent(in) :: SS
        real(8), intent(in) :: x, y
        real(8), intent(out) :: par(:)
        real(8) :: r, p_0

        r = sqrt(x**2 + y**2)

        ! par(1) = 1.0
        ! par(3) = SS%par_Velosity_inf
        ! par(4) = 0.0
        ! par(2) = 1.0
        ! par(5) = 100.0
        ! return

        if(size(par) < 5) then
            print*, "Error 50 Phys_Innitial_Conditions size(par) 9876tfghjdkofo43w8udyh4f"
            print*, size(par)
            STOP
        end if

        if(r < 80.0) then
            par(1) = 150.0 * (SS%par_R0/r)**2
            par(3) = x/r * 41.0
            par(4) = y/r * 41.0
            p_0 = 41.0**2 * 150.0/(SS%par_ggg * 7.0**2)
            par(2) = p_0 * (SS%par_R0/r)**(2.0 * SS%par_ggg)
            par(5) = 150.0 * (SS%par_R0/r)**2
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

    subroutine Calc_sourse_MF(SS, cell, sourse, step)  ! Считаются мультифлюидные источники
        TYPE (Setka), intent(in) :: SS
        real(8), intent(out) :: sourse(:)  ! (масса, два импульса и энергия)
        integer(4), intent(in) :: cell, step
        
        integer(4) :: i
        real(8) :: U_M_H(SS%n_Hidrogen), UU_H(SS%n_Hidrogen), sigma(SS%n_Hidrogen), nu(SS%n_Hidrogen)
        real(8) ro, p, u, v, ro_H, p_H, u_H, v_H
        
        sourse = 0.0

        ro = SS%gd(1, cell, step)
        p = SS%gd(2, cell, step)
        u = SS%gd(3, cell, step)
        v = SS%gd(4, cell, step)
        
        ! Body of Calc_sourse_MF
        do i = 1, SS%n_Hidrogen
            ro_H = SS%hydrogen(1, i, cell, step)
            p_H = SS%hydrogen(2, i, cell, step)
            u_H = SS%hydrogen(3, i, cell, step)
            v_H = SS%hydrogen(4, i, cell, step)

            if(ro_H <= 0.0) then
                ro_H = 0.0000001
            end if

            if(p_H <= 0.0) then
                p_H = 0.0000001
            end if

            U_M_H(i) = sqrt( (u - u_H)**2 + (v - v_H)**2 + &
            (64.0 / (9.0 * par_pi)) * (0.5 * p / ro + p_H / ro_H) )
            UU_H(i) = sqrt( (u - u_H)**2 + (v - v_H)**2 + &
            (4.0 / par_pi) * (0.5 * p / ro + p_H / ro_H) )
            sigma(i) = (1.0 - SS%par_a_2 * log(U_M_H(i)))**2
            nu(i) = ro * ro_H * U_M_H(i) * sigma(i)
        end do
        
        do i = 1, 4
            ro_H = SS%hydrogen(1, i, cell, step)
            p_H = SS%hydrogen(2, i, cell, step)
            u_H = SS%hydrogen(3, i, cell, step)
            v_H = SS%hydrogen(4, i, cell, step)
            sourse(1) = 0.0
            sourse(2) =  sourse(2) + nu(i) * (u_H - u)
            sourse(3) =  sourse(3) + nu(i) * (v_H - v)
            sourse(4) = sourse(4) + nu(i) * ( (u_H**2 + v_H**2 - &
                u**2 - v**2)/2.0 + (UU_H(i)/U_M_H(i)) * (p_H/ro_H - 0.5 * p/ro ) )
        end do
        
        sourse =  sourse * (SS%par_n_H_LISM/SS%par_Kn)


        ! if(ieee_is_nan(sourse(2))) then
        !     print*, "error source nan re5brtnmujtymuntr"
		! 	ro_H = SS%hydrogen(1, 3, cell, step)
        !     p_H = SS%hydrogen(2, 3, cell, step)
        !     u_H = SS%hydrogen(3, 3, cell, step)
        !     v_H = SS%hydrogen(4, 3, cell, step)
        !     print*, ro, p, u, v
        !     print*, ro_H, p_H, u_H, v_H
        !     print*, "______________"
        !     print*, nu
		! 	print*, "______________"
        !     print*, sigma
		! 	print*, "______________"
        !     print*, U_M_H
        !     print*, "______________"
        !     print*, sourse
        !     STOP
        ! end if
        
	end subroutine Calc_sourse_MF

end module Phys_parameter