module Phys_parameter
    use STORAGE 
    use GEOMETRY
    USE ieee_arithmetic
    implicit none 

    contains

    subroutine Inner_Conditions(SS, x, y, par)
        !! ����� ���������� ��������� �������
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
        !! ����� ��������� �������
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
        !! ����� ���������� �����
        TYPE (Setka), intent(in) :: SS
        real(8), intent(out) :: par(:)

        par(1) = 1.0
        par(3) = SS%par_Velosity_inf
        par(4) = 0.0
        par(2) = 1.0
        par(5) = 100.0

    end subroutine Phys_input_flow

    subroutine Get_gran_parameter(SS, gran, cell, par1_, par2_, now)
        !! ��������� ���������� �� ����� � ���� ������
        !! ������� ��������� �������� � ��� �������������
        TYPE (Setka), intent(in) :: SS
        integer(4), intent(in) :: gran, cell, now
        real(8), intent(out) :: par1_(:), par2_(:)

        real(8) :: c1(2), c2(2), c3(2), c4(2), c5(2)
        real(8) :: r1, r2, r3, r4, r5
        real(8) :: phi1, phi2, phi3, phi4, phi5, Vr, Vphi
        real(8) :: par1(5), par2(5), par3(5), par4(5)
        integer(4) :: s1, s2, s3, s4, s5
        logical :: istoch1, istoch2

        s1 = SS%gl_Gran_neighbour(1, gran)
        s2 = SS%gl_Gran_neighbour(2, gran)
        s3 = SS%gl_Gran_neighbour_TVD(1, gran)
        s4 = SS%gl_Gran_neighbour_TVD(2, gran)
        if(cell /= s1) then
            s5 = s1
            s1 = s2
            s2 = s5
            s5 = s3
            s3 = s4
            s4 = s5
        end if

        !! ����� ����� ����� ������ ��������
        if(s1 /= cell) then
            print*, "ERROR 987yhfjik3f4ofw3vtg54wbbwqvctracqcfrarvw"
            STOP
        end if

        if(s2 == -1) then
            call Phys_input_flow(SS, par2)
            par1_ = par2
            par2_ = par2
            return
        else if(s2 == -2) then
            par1_ = SS%gd(1:5, cell, now)
            if(par1_(3) > SS%par_Velosity_inf/3.0) par1_(3) = SS%par_Velosity_inf
            par2_ = par1_
            return
        else if(s2 == -3) then
            par1_ = SS%gd(1:5, cell, now)
            par2_ = par1_
            return
        else if(s2 == -4) then  !!  ��� ���������
            par1_ = SS%gd(1:5, s1, now)
            if(norm2(par1(3:4))/sqrt(SS%par_ggg * par1(2)/par1(1)) > 2.5) then
                c1 = SS%gl_Cell_Centr(:, s1, now)
                c5 = SS%gl_Gran_Center(:, gran, now)
                r1 = norm2(c1)
                r5 = norm2(c5)
                phi1 = polar_angle(c1(1), c1(2))
                phi5 = polar_angle(c5(1), c5(2))
                call polyar_skorost(phi1, par1_(3), par1_(4), Vr, Vphi)
                call dekard_polyar_skorost(phi5, Vr, Vphi, par1_(3), par1_(4))
                par1_(1) = par1_(1) * r1**2 / r5**2
                par1_(5) = par1_(5) * r1**2 / r5**2
                par1_(2) = par1_(2) * r1**(2.0 * SS%par_ggg) / r5**(2.0 * SS%par_ggg)
            end if
            par2_ = par1_
            par2_(4) = -par1_(4)
        end if

        par1 = SS%gd(1:5, s1, now)
        par2 = SS%gd(1:5, s2, now)

        istoch1 = .False.
        istoch2 = .False.
        if(norm2(par1(3:4))/sqrt(SS%par_ggg * par1(2)/par1(1)) > 2.5) istoch1 = .True.
        if(norm2(par2(3:4))/sqrt(SS%par_ggg * par2(2)/par2(1)) > 2.5) istoch2 = .True.

        c1 = SS%gl_Cell_Centr(:, s1, now)
        c2 = SS%gl_Cell_Centr(:, s2, now)
        c5 = SS%gl_Gran_Center(:, gran, now)

        r1 = norm2(c1)
        r2 = norm2(c2)
        r5 = norm2(c5)

        phi1 = polar_angle(c1(1), c1(2))
        phi2 = polar_angle(c2(1), c2(2))
        phi5 = polar_angle(c5(1), c5(2))

        if(istoch1 == .True. .and. istoch2 == .False.) then
            call polyar_skorost(phi1, par1(3), par1(4), Vr, Vphi)
            call dekard_polyar_skorost(phi5, Vr, Vphi, par1(3), par1(4))
            par1(1) = par1(1) * r1**2 / r5**2
            par1(5) = par1(5) * r1**2 / r5**2
            par1(2) = par1(2) * r1**(2.0 * SS%par_ggg) / r5**(2.0 * SS%par_ggg)
        else if(istoch1 == .False. .and. istoch2 == .True.) then
            call polyar_skorost(phi2, par2(3), par2(4), Vr, Vphi)
            call dekard_polyar_skorost(phi5, Vr, Vphi, par2(3), par2(4))
            par2(1) = par2(1) * r2**2 / r5**2
            par2(5) = par2(5) * r2**2 / r5**2
            par2(2) = par2(2) * r2**(2.0 * SS%par_ggg) / r5**(2.0 * SS%par_ggg)

        else if(istoch1 == .True. .and. istoch2 == .True.) then
            call polyar_skorost(phi1, par1(3), par1(4), Vr, Vphi)
            call dekard_polyar_skorost(phi5, Vr, Vphi, par1(3), par1(4))
            par1(1) = par1(1) * r1**2 / r5**2
            par1(5) = par1(5) * r1**2 / r5**2
            par1(2) = par1(2) * r1**(2.0 * SS%par_ggg) / r5**(2.0 * SS%par_ggg)

            call polyar_skorost(phi2, par2(3), par2(4), Vr, Vphi)
            call dekard_polyar_skorost(phi5, Vr, Vphi, par2(3), par2(4))
            par2(1) = par2(1) * r2**2 / r5**2
            par2(5) = par2(5) * r2**2 / r5**2
            par2(2) = par2(2) * r2**(2.0 * SS%par_ggg) / r5**(2.0 * SS%par_ggg)
        else
            
        end if

        par1_ = par1
        par2_ = par2

    end subroutine Get_gran_parameter

    subroutine Calc_sourse_MF(SS, cell, sourse, step)  ! ��������� �������������� ���������
        TYPE (Setka), intent(in) :: SS
        real(8), intent(out) :: sourse(:)  ! (�����, ��� �������� � �������)
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