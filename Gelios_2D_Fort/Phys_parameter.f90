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

    subroutine Get_gran_parameter(SS, gran, cell, par1_, par2_, now)
        !! Получение параметров на грани с двух сторон
        !! Функция учитывает источник и ТВД распределения
        TYPE (Setka), intent(in) :: SS
        integer(4), intent(in) :: gran, cell, now
        real(8), intent(out) :: par1_(:), par2_(:)

        real(8) :: c1(2), c2(2), c3(2), c4(2), c5(2)
        real(8) :: r1, r2, r3, r4, r5
        real(8) :: phi1, phi2, phi3, phi4, phi5, Vr, Vphi, d1, d2, d3
        real(8) :: par1(5), par2(5), par3(5), par4(5)
        integer(4) :: s1, s2, s3, s4, s5, i
        logical :: istoch1, istoch2

        ! if(gran == 77) then
        !     continue
        ! end if

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

        !! Потом можно будет убрать проверку
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
        else if(s2 == -4) then  !!  Ось симметрии
            par1_ = SS%gd(1:5, s1, now)
            if(norm2(par1_(3:4))/sqrt(SS%par_ggg * par1_(2)/par1_(1)) > 2.5) then
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
                par2_ = par1_
                par2_(4) = -par1_(4)
                return
            end if
            par2_ = par1_

            ! c1 = SS%gl_Cell_Centr(:, s1, now)
            ! c5 = SS%gl_Gran_Center(:, gran, now)
            ! c3 = SS%gl_Cell_Centr(:, s3, now)
            ! d1 = -norm2(c5 - c1)
            ! d2 = -norm2(c5 - c3)
            ! par3 = SS%gd(1:5, s3, now)
            ! do i = 4, 4 
            !     par1_(i) = linear(d2, par3(i), d1, par1_(i), -d1, -par1_(i), 0.0_8)
            ! end do
            ! par2_(4) = -par1_(4)


            ! par1_(4) = 0.0
            ! par2_(4) = -par1_(4)


            !! Метод сноса в одну сторону

            c1 = SS%gl_Cell_Centr(:, s1, now)
            c5 = SS%gl_Gran_Center(:, gran, now)
            c3 = SS%gl_Cell_Centr(:, s3, now)
            d1 = -norm2(c5 - c1)
            d2 = -norm2(c5 - c3)
            par3 = SS%gd(1:5, s3, now)
            do i = 1, 5 
                par1_(i) = linear1(d2, par3(i), d1, par1_(i), 0.0_8)
            end do

            par2_ = par1_
            par2_(4) = -par1_(4)

            return
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
            if(SS%gl_all_Cell_zone(s1) == SS%gl_all_Cell_zone(s2)) then
                if(s3 > 0 .and. s4 > 0) then
                    c3 = SS%gl_Cell_Centr(:, s3, now)
                    c4 = SS%gl_Cell_Centr(:, s4, now)
                    d1 = -norm2(c5 - c1)
                    d2 = -norm2(c5 - c3)
                    d3 = norm2(c5 - c2)
                    par3 = SS%gd(1:5, s3, now)
                    par4 = SS%gd(1:5, s4, now)

                    do i = 1, 5 
                        par1_(i) = linear(d2, par3(i), d1, par1(i), d3, par2(i), 0.0_8)
                    end do

                    d2 = -d1
                    d1 = -d3
                    d3 = d2
                    d2 = -norm2(c4 - c5)

                    do i = 1, 5 
                        par2_(i) = linear(d2, par4(i), d1, par2(i), d3, par1(i), 0.0_8)
                    end do
                    
                    return
                else if(SS%gl_Cell_type(s1) == "A" .and. SS%gl_Cell_number(2, s1) == 1 .and. SS%gl_Cell_number(2, s2) == 2 .and. s4 > 0) then
                    !! Для ячеек рядом с осью симметрии (проекция вверх)
                    !print*, ")))))))))))))))))))"
                    !pause
                    par1_ = par1
                    !par1_(4) = par1_(4) * c5(2)/c1(2)
                    c4 = SS%gl_Cell_Centr(:, s4, now)
                    d1 = -norm2(c5 - c1)
                    d3 = norm2(c5 - c2)
                    d2 = 2.0 * d1! -norm2(c5 - (/c1(1), -c1(2)/) )
                    par4 = SS%gd(1:5, s4, now)

                    do i = 4, 4 
                        par1_(i) = linear(d2, 0.0_8, d1, par1(i), d3, par2(i), 0.0_8)
                        !par1_(i) = linear(d2, -par1(i), d1, par1(i), d3, par2(i), 0.0_8)
                    end do

                    c4 = SS%gl_Cell_Centr(:, s4, now)
                    d1 = -norm2(c5 - c1)
                    d3 = norm2(c5 - c2)
                    !par4 = SS%gd(1:5, s4, now)

                    d2 = -d1
                    d1 = -d3
                    d3 = d2
                    d2 = -norm2(c4 - c5)

                    do i = 1, 5 
                        par2_(i) = linear(d2, par4(i), d1, par2(i), d3, par1(i), 0.0_8)
                    end do
                    
                    return
                else if(SS%gl_Cell_type(s2) == "A" .and. SS%gl_Cell_number(2, s2) == 1 .and. SS%gl_Cell_number(2, s1) == 2 .and. s3 > 0) then
                    !! Для ячеек рядом с осью симметрии (проекция вверх)
                    !print*, ")))))))))))))))))))"
                    !pause
                    par2_ = par2
                    !par2_(4) = par2_(4) * c5(2)/c2(2)

                    c3 = SS%gl_Cell_Centr(:, s3, now)
                    d1 = -norm2(c5 - c2)
                    d2 = 2.0 * d1!-norm2(c5 - (/c2(1), -c2(2)/) )
                    d3 = norm2(c5 - c1)
                    !par4 = SS%gd(1:5, s4, now)

                    do i = 4, 4 
                        !par2_(i) = linear(d2, -par2(i), d1, par2(i), d3, par1(i), 0.0_8)
                        par2_(i) = linear(d2, 0.0_8, d1, par2(i), d3, par1(i), 0.0_8)
                    end do


                    !c3 = SS%gl_Cell_Centr(:, s3, now)
                    d1 = -norm2(c5 - c1)
                    d2 = -norm2(c5 - c3)
                    d3 = norm2(c5 - c2)
                    par3 = SS%gd(1:5, s3, now)

                    do i = 1, 5 
                        par1_(i) = linear(d2, par3(i), d1, par1(i), d3, par2(i), 0.0_8)
                    end do
                    
                    return

                end if
            end if
        end if

        ! if(c5(1) > 22 .and. c5(1) < 200 .and. c5(2) < 100 .and. SS%gl_all_Cell_zone(s1) == SS%gl_all_Cell_zone(s2)) then
        !     print*, "Center = "
        !     print*, c5
        !     print*, "_____________________________"
        !     print*, "Center = "
        !     print*, c1
        !     print*, "_____________________________"
        !     print*, "Center = "
        !     print*, c2
        !     print*, "_____________________________"
        !     pause
        ! end if

        par1_ = par1
        par2_ = par2

    end subroutine Get_gran_parameter

    subroutine Calc_sourse_MF(SS, cell, sourse, step, use_koeff_)  ! Считаются мультифлюидные источники
        TYPE (Setka), intent(in) :: SS
        real(8), intent(out) :: sourse(:)  ! (масса, два импульса и энергия)
        integer(4), intent(in) :: cell, step
        logical, intent(in), optional :: use_koeff_     ! Оспользовать ли умножение источников на коэффециенты, посчитанные Монте-Карло
        logical :: use_koeff
        integer(4) :: i
        real(8) :: U_M_H(SS%n_Hidrogen), UU_H(SS%n_Hidrogen), sigma(SS%n_Hidrogen), nu(SS%n_Hidrogen)
        real(8) ro, p, u, v, ro_H, p_H, u_H, v_H

        use_koeff = .False.
        if(PRESENT(use_koeff_))  use_koeff = use_koeff_
        
        sourse = 0.0

        ro = SS%gd(1, cell, step)
        p = SS%gd(2, cell, step)
        u = SS%gd(3, cell, step)
        v = SS%gd(4, cell, step)

        if(ro <= 0.0000001) then
            ro = 0.0000001
        end if

        if(p <= 0.0000001) then
            p = 0.0000001
        end if
        
        ! Body of Calc_sourse_MF
        do i = 1, SS%n_Hidrogen
            ro_H = SS%hydrogen(1, i, cell, step)
            p_H = SS%hydrogen(2, i, cell, step)
            u_H = SS%hydrogen(3, i, cell, step)
            v_H = SS%hydrogen(4, i, cell, step)

            if(ro_H <= 0.0000001) then
                ro_H = 0.0000001
            end if

            if(p_H <= 0.0000001) then
                p_H = 0.0000001
            end if

            U_M_H(i) = sqrt( (u - u_H)**2 + (v - v_H)**2 + &
            (64.0 / (9.0 * par_pi)) * (p / ro + 2.0 * p_H / ro_H) )
            UU_H(i) = sqrt( (u - u_H)**2 + (v - v_H)**2 + &
            (4.0 / par_pi) * (p / ro + 2.0 * p_H / ro_H) )
            sigma(i) = (1.0 - SS%par_a_2 * log(U_M_H(i)))**2
            nu(i) = ro * ro_H * U_M_H(i) * sigma(i)
        end do
        
        do i = 1, 4
            ro_H = SS%hydrogen(1, i, cell, step)
            p_H = SS%hydrogen(2, i, cell, step)
            u_H = SS%hydrogen(3, i, cell, step)
            v_H = SS%hydrogen(4, i, cell, step)

            if(ro_H <= 0.0000001) then
                ro_H = 0.0000001
            end if

            if(p_H <= 0.0000001) then
                p_H = 0.0000001
            end if

            sourse(1) = 0.0
            sourse(2) =  sourse(2) + nu(i) * (u_H - u)
            sourse(3) =  sourse(3) + nu(i) * (v_H - v)
            sourse(4) = sourse(4) + nu(i) * ( (u_H**2 + v_H**2 - &
                u**2 - v**2)/2.0 + (UU_H(i)/U_M_H(i)) * (2.0 * p_H/ro_H - p/ro ) )
        end do
        
        if (use_koeff == .False.) then
            sourse =  sourse * (SS%par_n_H_LISM/SS%par_Kn)
        else
            sourse(1) = SS%atom_source(4, cell)
            sourse(2) = sourse(2) * (SS%par_n_H_LISM/SS%par_Kn) * SS%atom_source(1, cell)
            sourse(4) = sourse(4) * (SS%par_n_H_LISM/SS%par_Kn) * SS%atom_source(3, cell)

            ! if(SS%gl_Cell_type(cell) == 'A' .and. SS%gl_Cell_number(2, cell) <= 2) then !! УБРАТЬ
            !     sourse(3) = sourse(3) * (SS%par_n_H_LISM/SS%par_Kn)
            ! else
            !     sourse(3) = sourse(3) * (SS%par_n_H_LISM/SS%par_Kn) * SS%atom_source(2, cell)
            ! end if
            sourse(3) = sourse(3) * (SS%par_n_H_LISM/SS%par_Kn) * SS%atom_source(2, cell)
        end if
        
	end subroutine Calc_sourse_MF

    subroutine Calc_Pogloshenie(SS)
        ! Печатаем поверхности, которые выделяем
        TYPE (Setka), intent(in) :: SS
        integer(4) :: i, j, i_luch, sort, cell
        real(8) :: alf, dl, r(2), u
        character(len=2) :: name
        real(8) :: pogl(SS%n_Hidrogen, SS%pogl_iter)
        LOGICAL :: inzone

        dl = 0.01
        
        do i_luch = 1, 1!19
            pogl = 0.0
            r = 0.0
            cell = 1
            write(unit=name,fmt='(i2.2)') i_luch
            open(1, file = "Pogloshenie_" // name // ".txt")

            alf = (i_luch - 1) * (par_pi/18)

            do while (.True.)
                r = r + dl * (/ cos(alf), sin(alf) /)
                call Geo_Find_Cell(SS, r(1), r(2), cell, inzone)
                if(inzone == .False.) EXIT

                do sort = 1, SS%n_Hidrogen
                    do i = 1, SS%pogl_iter
                        pogl(sort, i) = pogl(sort, i) + SS%par_n_H_LISM * dl * SS%pogloshenie(sort, i, cell)/SS%pogl_ddd
                    end do
                end do

            end do

            write(1,*) "TITLE = 'HP'  VARIABLES = 'u', 'f1', 'f2', 'f3', 'f4', 'ff'"
            pogl = pogl * par_poglosh

            do i = 1, SS%pogl_iter
                u = SS%pogl_v_min + (i + 0.5) * SS%pogl_ddd;
                write(1,*) u * 10.3804, exp(-pogl(1, i)), exp(-pogl(2, i)), exp(-pogl(3, i)), exp(-pogl(4, i)), exp(-sum(pogl(:, i)))
                !! Здесь скорость переведена в км/с
            end do

            close(1)
        end do
    end subroutine Calc_Pogloshenie

end module Phys_parameter