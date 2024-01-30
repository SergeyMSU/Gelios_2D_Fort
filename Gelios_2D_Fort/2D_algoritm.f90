module Algoritm
    USE STORAGE
    USE GEOMETRY
    USE Phys_parameter
	USE SURFACE
    USE Solvers
    USE ieee_arithmetic
    USE OMP_LIB
    USE Interpol
    USE cgod
    USE Monte_Karlo
    implicit none 

    contains

    subroutine Gas_dynamic_algoritm(SS)
        TYPE (Setka), intent(in out) :: SS
        integer(4) :: num, i, i_max
        real(8) :: par(5), parH(5, 4)

        call Read_setka_bin(gl_S3, "00002")   ! ДЛЯ ВОДОРОДА
        call Int_Init(gl_S2, gl_S3)

        call Dell_Setka(gl_S3)

        call Int_Print_Cell(gl_S2)

        call Read_setka_bin(SS, "00034")      ! ОСНОВНАЯ СЕТКА

        call Geo_Set_sxem(SS)
        



        call Algoritm_Reinterpol(SS, gl_S2)

        if(SS%init_geo == .False.) then
            STOP "Gas_dynamic_algoritm  error init_geo erty67543"  
        end if

        num = 1

        !call Algoritm_Initial_condition(SS)  ! Зададим начальные условия для всех ячеек сетки
        !call Geo_get_request(SS)
        call Algoritm_Bound_condition(SS)    ! Зададим граничные условия на внутренней сфере

        ! Проверка переменных
        print*, "____________________________________________"
        print*, "Proverka peremennix"
        call Geo_Find_Cell(SS, 250.0_8, 1.0_8, num)
        print*, "num = ", num
        print*, "PAR ", SS%gd(:, num, 1)
        print*, "____________________________________________"
        print*, "PAR ", SS%gd(:, num, 2)
        print*, "____________________________________________"
        print*, "PAR HYDROGEN 4", SS%hydrogen(:, 4, num, 1)
        print*, "____________________________________________"
        print*, "PAR HYDROGEN 4", SS%hydrogen(:, 4, num, 2)
        print*, "____________________________________________"

        ! call Int_Init(gl_S2, SS)
        ! call Int_Print_Cell(gl_S2)
        
        !call Int_Get_Parameter(gl_S2, 260.0_8, 19.5_8, num, PAR_gd = par, PAR_hydrogen = parH)

        ! call Calc_move_velosity(SS, 1)
        ! call Move_all(SS, 2, 1.0_8)



        i_max = 100!200!350
        do i = 1, i_max
            if (mod(i, 5) == 0) then
                print*, "Global step = ", i, "from ", i_max
            end if
            call Start_GD_algoritm(SS, 5000, 2) !5000

            call Culc_Cell_Centr(SS, 1)
            call Geo_Culc_normal(SS, 1) 
            call Geo_Culc_length_area(SS, 1)
            call Culc_Cell_Centr(SS, 2)
            call Geo_Culc_normal(SS, 2) 
            call Geo_Culc_length_area(SS, 2)

            call Start_GD_algoritm(SS, 500, 1) !500

            call Algoritm_Reinterpol(SS, gl_S2)
        end do

        call Print_GD(SS)
        call Geo_Print_Surface(SS)
        call Save_setka_bin(SS, "00035")
        call Print_Grans(SS)
        ! call Print_Cell_Centr(SS)
        call Print_GD_1D(SS)
        ! call Print_TVD_Sosed(SS)

        !pause
    end subroutine Gas_dynamic_algoritm

    subroutine MK_algoritm(SS)
        TYPE (Setka), intent(in out) :: SS

        call Read_setka_bin(SS, "00034")      ! ОСНОВНАЯ СЕТКА
		print*, "A1"
        call SUR_init(gl_surf1, SS)
		print*, "A2"
        call Int_Init(gl_S2, SS)
		print*, "A3"

        !! Задаём параметры мини-сетки
        gl_S3%par_m_A = 20! 30      ! Количество лучей A в плоскости
        gl_S3%par_m_BC = 10! 18      ! Количество лучей B/C в плоскости
        gl_S3%par_m_O = 10! 17      ! Количество лучей O в плоскости
        gl_S3%par_m_K = 8! 7      ! Количество лучей K в плоскости
        gl_S3%par_n_TS =  27! 26                    ! Количество точек до TS (TS включается)
        gl_S3%par_n_HP =  37! 40                 ! Количество точек HP (HP включается)  всё от 0 считается
        gl_S3%par_n_BS =  57! 60! 5                 ! Количество точек BS (BS включается)
        gl_S3%par_n_END = 65! 72! 6                ! Количество точек до конца сетки (конец включается)
        gl_S3%par_n_IA =  12! 12                   ! Количество точек, которые входят во внутреннюю область
        gl_S3%par_n_IB =  14! 14                   ! Количество точек, которые входят во внутреннюю область (с зазором)
        call Init_Setka(gl_S3)
		print*, "A4"
        call Build_Setka_start(gl_S3)
		print*, "A5"
	    call Algoritm_ReMove_Surface(gl_S3, gl_surf1)
		print*, "A6"

        !print*, gl_S2%gd(:, 1)
        !print*, gl_S2%gd(:, 2)
        !print*, gl_S2%gd(:, 3)

        print*, "A7"

        call Algoritm_Reinterpol(gl_S3, gl_S2)
		print*, gl_S3%gd(:, 1, 1)
        print*, gl_S3%gd(:, 2, 1)
        print*, gl_S3%gd(:, 3, 1)
		print*, "A8"

        call Print_Cell(gl_S3)

        call M_K_start(gl_S3, gl_S2)

        call Print_hydrogen(gl_S3)
        call Print_hydrogen_1D(gl_S3)

        call Save_setka_bin(gl_S3, "B0034")

        print*, "END"
        pause 

    end subroutine MK_algoritm

    subroutine Start_GD_algoritm(SS, all_step_, area)
        TYPE (Setka), intent(in out) :: SS
        integer(4), intent(in) :: all_step_
        integer(4), intent(in) :: area ! 1 или 2
        ! 1 - внутренняя сфера
        ! 2 - внешняя область

        integer(4) :: now, now2, step, cell, Ncell, gr, gran, sosed, all_step, TVD_sosed_1, TVD_sosed_2
        real(8) :: time, TT, center(2), par1(5), par2(5), ro, p, u, v, Q, normal(2), Sqv, Vol, Vol2
        real(8) :: ro2, p2, u2, v2, Q2, pp, par1_TVD(5), r, phi1, phi2, Vr, Vphi, phi3
        real(8) :: qqq1(9), qqq2(9), POTOK2(9), POTOK(5), source(4), r2, r3, par3(5), par4(5)
        real(8) :: dsl, dsc, dsp, loc_time, ALL_TIME, lenght, wc, gran_center(2), sosed_center(2)
        logical :: tvd1, tvd2
        integer(4) :: kdir, idgod, KOBL

		all_step = all_step_
        if(mod(all_step, 2) == 1) all_step = all_step + 1   ! Делаем общее число шагов чётным

        now = 2                        ! Какие параметры сейчас будут считаться (1 или 2). Они меняются по очереди
        time = 0.000002_8               ! Начальная инициализация шага по времени 
        Ncell = size(SS%gl_all_Cell(1, :))
        qqq1 = 0.0
        qqq2 = 0.0
        ALL_TIME = 0.0

        do step = 1, all_step

            if (mod(step, 12000) == 0) then
                print*, "Step = ", step, " dt = ", time, " All_time = ", ALL_TIME
				!pause
            end if

            now2 = now
            now = mod(now, 2) + 1

            TT = time
            time = 1000000.0
            ALL_TIME = ALL_TIME + TT

            if(area == 2) then
                call Calc_move_velosity(SS, now)
                call Move_all(SS, now, TT)
                call Culc_Cell_Centr(SS, now2)
                call Geo_Culc_normal(SS, now2) 
                call Geo_Culc_length_area(SS, now2)
            end if

            ! Предлагаю пробегаться сразу по ячейкам, а не по граням

            !$omp parallel


            !$omp do private(KOBL, kdir, idgod, sosed_center, phi3, gran_center, Vr, Vphi, phi1, phi2, r, par1_TVD, wc, &
            !$omp r3, r2, loc_time, gr, gran, sosed, center, par1, TVD_sosed_1, TVD_sosed_2, &
            !$omp par2, ro, p, u, v, Q, normal, Sqv, Vol, Vol2, ro2, p2, u2, v2, Q2, pp, qqq1, &
            !$omp qqq2, POTOK2, POTOK, source, dsl, dsc, dsp, lenght, par3, par4, tvd1, tvd2)
            do cell = 1, Ncell
                source = 0.0
                center = SS%gl_Cell_Centr(:, cell, now)
                r = norm2(center)
                phi1 = polar_angle(center(1), center(2))

                if(area == 2) then
                    if(r < 10.0) CYCLE
                else if(area == 1) then
                    if(r >= 10.0 .or. r <= 0.42) CYCLE
                end if

                par1 = SS%gd(1:5, cell, now)    ! Получили газодинамические параметры
                POTOK = 0.0
                
                if(par1(1) <= 0.000000001) then
                    print*, "Error rho 89 45465u76jhgrefrwcwdf4c, ", par1(1) 
                    pause
                end if


                ! Пробегаемся по граням
                do gr = 1, 4
                    ! tvd1 = .True.
                    ! tvd2 = .True.
					POTOK2 = 0.0
                    KOBL = 0
                    kdir = 0
                    idgod = 0
                    gran = SS%gl_Cell_gran(gr, cell)
                    sosed = SS%gl_Cell_neighbour(gr, cell)
                    if(sosed == 0) CYCLE
                    gran_center = SS%gl_Gran_Center(:, gran, now)


                    ! r2 = norm2(gran_center)
					! phi2 = polar_angle(gran_center(1), gran_center(2))

                    ! ! Определяем TVD-соседей (они не понадобятся в гиперзвуке)
                    ! TVD_sosed_1 = SS%gl_Gran_neighbour_TVD(1, gran)
                    ! if(SS%gl_Gran_neighbour(1, gran) /= cell) then
                    !     TVD_sosed_2 = TVD_sosed_1
                    !     TVD_sosed_1 = SS%gl_Gran_neighbour_TVD(2, gran)
                    ! else
                    !     TVD_sosed_2 = SS%gl_Gran_neighbour_TVD(2, gran)
                    ! end if

                    ! !! Условия в гиперзвуке
                    ! if(r < 60 .and. center(1) < 30.0 .and. norm2(par1(3:4))/sqrt(SS%par_ggg * par1(2)/par1(1)) > 2.5 ) then
                    !     par1_TVD = par1
                    !     call polyar_skorost(phi1, par1(3), par1(4), Vr, Vphi)
                    !     call dekard_polyar_skorost(phi2, Vr, Vphi, par1_TVD(3), par1_TVD(4))
                    !     par1_TVD(1) = par1(1) * r**2 / r2**2
                    !     par1_TVD(5) = par1(5) * r**2 / r2**2
                    !     par1_TVD(2) = par1(2) * r**(2.0 * SS%par_ggg) / r2**(2.0 * SS%par_ggg)
                    !     tvd1 = .False.
                    ! else
                    !     par1_TVD = par1
                    ! end if

                    ! ! Задаём параметры соседа для распада разрыва
                    ! if(sosed == 0) then
                    !     CYCLE
                    ! else if(sosed == -1) then
                    !     call Phys_input_flow(SS, par2)
                    ! else if(sosed == -2) then
                    !     par2 = par1
                    !     if(par2(3) > SS%par_Velosity_inf/2.0) par2(3) = SS%par_Velosity_inf
                    ! else if(sosed == -3) then
                    !     par2 = par1
                    ! else if(sosed == -4) then
                    !     par2 = par1_TVD
                    !     par2(4) = -par1_TVD(4)
                    ! else
                    !     par2 = SS%gd(1:5, sosed, now)    ! Получили газодинамические параметры
                    !     if(r < 60 .and. center(1) < 30.0 .and. norm2(par2(3:4))/sqrt(SS%par_ggg * par2(2)/par2(1)) > 2.5 ) then
                    !         tvd2 = .False.
                    !         sosed_center = SS%gl_Cell_Centr(:, sosed, now)
                    !         r3 = norm2(sosed_center)
                    !         phi3 = polar_angle(sosed_center(1), sosed_center(2))
                    !         call polyar_skorost(phi3, par2(3), par2(4), Vr, Vphi)
                    !         call dekard_polyar_skorost(phi2, Vr, Vphi, par2(3), par2(4))
                    !         par2(1) = par2(1) * r3**2 / r2**2
                    !         par2(5) = par2(5) * r3**2 / r2**2
                    !         par2(2) = par2(2) * r3**(2 * SS%par_ggg) / r2**(2 * SS%par_ggg)
                    !     else
                    !         !! Здесь надо делать ТВД

                    !     end if
                    ! end if

                    ! if(par2(1) <= 0.000000001) then
                    !     print*, "Error rho2 89 45465u76jhgrefrwcwdf4c, ", par2(1) 
                    !     print*, "num = ", cell, sosed
                    !     print*, "center = ", center
                    !     print*, "_______________"
                    !     print*, par2
                    !     print*, "_______________"
                    !     print*, par1
                    !     print*, "_______massiv________"
                    !     print*, SS%gd(1:5, sosed, 1)
                    !     print*, "_______________"
                    !     print*, SS%gd(1:5, sosed, 2)
                    !     pause
                    ! end if

                    normal = SS%gl_Gran_normal(:, gran, now)
                    Sqv = SS%gl_Gran_length(gran, now)
                    lenght = SS%gl_Cell_gran_dist(gr, cell, now)

                    if(SS%gl_Gran_neighbour(1, gran) /= cell) normal = -normal

                    ! Нужно вычислить скорость движения грани
                    wc = DOT_PRODUCT((SS%gl_Gran_Center(:, gran, now2) -  gran_center)/TT, normal)

                    call Get_gran_parameter(SS, gran, cell, par1_TVD, par2, now)

                    qqq1(1) = par1_TVD(1)
                    qqq1(5) = par1_TVD(2)
                    qqq1(2) = par1_TVD(3)
                    qqq1(3) = par1_TVD(4)
                    qqq1(4) = 0.0
                    qqq1(6:8) = 0.0
                    qqq1(9) = par1_TVD(5)

                    qqq2(1) = par2(1)
                    qqq2(5) = par2(2)
                    qqq2(2) = par2(3)
                    qqq2(3) = par2(4)
                    qqq2(4) = 0.0
                    qqq2(6:8) = 0.0
                    qqq2(9) = par2(5)

                    !if(.False.) then
                    !if(.True.) then
                    if(SS%gl_Gran_shem(gran) == 3) then
                        call cgod3d(KOBL, 0, 0, 0, kdir, idgod, &
                        normal(1), normal(2), 0.0_8, 1.0_8, &
                        wc, qqq1(1:8), qqq2(1:8), &
                        dsl, dsp, dsc, 1.0_8, 1.66666666666666_8, &
                        POTOK2)

                        if (idgod == 2) then
                            POTOK2 = 0.0
                            call chlld_Q(1, normal(1), normal(2), 0.0_8, &
                                wc, qqq1, qqq2, dsl, dsp, dsc, POTOK2, .False.)
                        else
                            if(dsc >= wc) then  !! Правильно считаем конвективный перенос
                                POTOK2(9) = POTOK2(1) / qqq1(1) * qqq1(9)
                            else
                                POTOK2(9) = POTOK2(1) / qqq2(1) * qqq2(9)
                            end if
                        end if
                    else
                        call chlld_Q(SS%gl_Gran_shem(gran), normal(1), normal(2), 0.0_8, &
                        wc, qqq1, qqq2, dsl, dsp, dsc, POTOK2, .False.)
                    end if

                    loc_time = 0.9 * lenght/( max(dabs(dsl), dabs(dsp)) + dabs(wc) )

                    POTOK(5) = POTOK(5) + POTOK2(9) * Sqv
                    POTOK(1) = POTOK(1) + POTOK2(1) * Sqv
                    POTOK(3) = POTOK(3) + POTOK2(2) * Sqv
                    POTOK(4) = POTOK(4) + POTOK2(3) * Sqv
                    POTOK(2) = POTOK(2) + POTOK2(5) * Sqv

                    if(loc_time < time) then
                        !$omp critical
                            time = loc_time
                        !$omp end critical
                    end if

                    ! if(cell == 2807) then
                    !    print*, "Gran = ", gr, sosed
                    !    print*, "sosed = ", SS%gl_Gran_neighbour(1, gran), SS%gl_Gran_neighbour(2, gran), cell
                    !    print*, "normal = ", normal
                    !    print*, "Sqv = ", Sqv
                    !    print*, "POTOK2 = ", POTOK2
                    ! end if

                end do

                Vol = SS%gl_Cell_square(cell, now)
                Vol2 = SS%gl_Cell_square(cell, now2)

                call Calc_sourse_MF(SS, cell, source, now)

                if(ieee_is_nan(source(2))) then
                    print*, "error source nan 186 tyujhwgeftywfwf"
                    print*, "centr = ", center
                    print*, "______________ 0"
                    print*, par1
                    print*, "______________ 1"
                    print*, SS%hydrogen(:, 1, cell, now)
                    print*, "______________ 2"
                    print*, SS%hydrogen(:, 2, cell, now)
                    print*, "______________ 3"
                    print*, SS%hydrogen(:, 3, cell, now)
                    print*, "______________ 4"
                    print*, SS%hydrogen(:, 4, cell, now)
                    print*, "______________"
                    print*, source
                    STOP
                end if

                ro = par1(1)
                p = par1(2)
                u = par1(3)
                v = par1(4)
                Q = par1(5)

                ! Законы сохранения в ячейке
                ro2 = ro * Vol/Vol2 - TT * (POTOK(1) / Vol2 + ro * v/center(2))
                if(ro2 <= 0.0) then
                    print*, "Ro < 0", ro2, ro, TT, Vol, cell
                    print*, "centr = ", center
                    print*, "_____________"
                    print*, "POTOK = ", POTOK
                    print*, "_____________"
                    print*, par1, "||||||||| ", par2
                    print*, "_____________"
                    stop
                end if
                Q2 = Q * Vol/Vol2 - TT * (POTOK(5) / Vol2 + Q * v/center(2))
                u2 = (ro * u * Vol/Vol2 - TT * ( POTOK(3) / Vol2 + ro * v * u/center(2) - source(2) )) / ro2
                v2 = (ro * v * Vol/Vol2 - TT * ( POTOK(4) / Vol2 + ro * v * v/center(2) - source(3) )) / ro2

                pp = v * (SS%par_ggg * p / (SS%par_ggg - 1.0) + ro * (u * u + v * v) * 0.5) / center(2)
                p2 = ((  ( p / (SS%par_ggg - 1.0) + 0.5 * ro * (u**2 + v**2))  * Vol/Vol2   &
                        - TT * (POTOK(2)/ Vol2 + pp - source(4)) ) - 0.5 * ro2 * (u2**2 + v2**2) ) * (SS%par_ggg - 1.0)

                if(p2 <= 0.0) then
                    p2 = 0.000001
                end if

                SS%gd(1, cell, now2) = ro2
                SS%gd(2, cell, now2) = p2
                SS%gd(3, cell, now2) = u2
                SS%gd(4, cell, now2) = v2
                SS%gd(5, cell, now2) = Q2

			end do
			!$omp end do
            !$omp end parallel

        end do

        ! print*, "ALL_TIME = ", ALL_TIME

    end subroutine Start_GD_algoritm

    subroutine Calc_move_velosity(SS, step)
        !! Вычисляем скорости движения граней
        TYPE (Setka), intent(in out) :: SS
        integer(4), intent(in) :: step
        integer(4) :: Num, i, s1, s2, gran, Num2, Num3, s3
        real(8) :: normal(2), qqq1(8), qqq2(8), POTOK(8)
        real(8) :: dsl, dsc, dsp
        real(8) :: koeff_TS, koeff_HP, koeff_BS, c1(2), c2(2), c3(2), wc
        integer(4) :: kdir, idgod, KOBL

        koeff_TS = 0.002
        koeff_HP = 0.1
        koeff_BS = 0.1

        Num = size(SS%gl_TS)

        SS%gl_Point_num = 0
        SS%gl_yzel_Vel = 0.0

        Num = size(SS%gl_TS)
        Num2 = size(SS%gl_HP)
        Num3 = size(SS%gl_BS)

        !$omp parallel
        !$omp do private(KOBL, wc, kdir, idgod, gran, normal, s1, s2, s3, qqq1, qqq2, POTOK, dsl, dsc, dsp, c1, c2, c3)
        do i = 1, Num
            qqq1 = 0.0
            qqq2 = 0.0
            wc = 0.0
            KOBL = 0
            kdir = 0
            idgod = 0
            gran = SS%gl_TS(i)
            normal = SS%gl_Gran_normal(:, gran, step)
            s1 = SS%gl_Gran_neighbour(1, gran)
            s2 = SS%gl_Gran_neighbour(2, gran)

            qqq1(1) = SS%gd(1, s1, step)
            qqq1(2) = SS%gd(3, s1, step)
            qqq1(3) = SS%gd(4, s1, step)
            qqq1(5) = SS%gd(2, s1, step)

            qqq2(1) = SS%gd(1, s2, step)
            qqq2(2) = SS%gd(3, s2, step)
            qqq2(3) = SS%gd(4, s2, step)
            qqq2(5) = SS%gd(2, s2, step)


            call cgod3d(KOBL, 0, 0, 0, kdir, idgod, &
                normal(1), normal(2), 0.0_8, 1.0_8, &
                wc, qqq1(1:8), qqq2(1:8), &
                dsl, dsp, dsc, 1.0_8, 1.66666666666666_8, &
                POTOK)

            if (idgod == 2) then
                call chlld(2, normal(1), normal(2), 0.0_8, &
				0.0_8, qqq1, qqq2, dsl, dsp, dsc, POTOK)
            end if

            normal = normal * dsl * koeff_TS
            s1 = SS%gl_all_Gran(1, gran)
            s2 = SS%gl_all_Gran(2, gran)

            !SS%gl_yzel_Vel(:, s1) = SS%gl_yzel_Vel(:, s1) + normal
            SS%gl_yzel_Vel(:, s2) = SS%gl_yzel_Vel(:, s2) + normal

            !SS%gl_Point_num(s1) = SS%gl_Point_num(s1) + 1
            SS%gl_Point_num(s2) = SS%gl_Point_num(s2) + 1

            !! Сюда же напишем поверхностное натяжение
            !! Это поверхностое натяжение работает только в случае, если скорость движения узла везда как среднее движения его граней
            s1 = SS%gl_all_Gran(1, gran)
            s2 = SS%gl_all_Gran(2, gran)

            c1 = SS%gl_yzel(:, s1, step)
            c2 = SS%gl_yzel(:, s2, step)

            if(i > 1) then
                gran = SS%gl_TS(i - 1)
                s3 = SS%gl_all_Gran(2, gran)
                c3 = SS%gl_yzel(:, s3, step)
                c3 = (c2 + c3)/2.0
                SS%gl_yzel_Vel(:, s1) = SS%gl_yzel_Vel(:, s1) + 1.0 * (c3 - c1) * SS%par_nat_TS
            end if

            if(i < Num) then
                gran = SS%gl_TS(i + 1)
                s3 = SS%gl_all_Gran(1, gran)
                c3 = SS%gl_yzel(:, s3, step)
                c3 = (c1 + c3)/2.0
                SS%gl_yzel_Vel(:, s2) = SS%gl_yzel_Vel(:, s2) + 1.0 * (c3 - c2) * SS%par_nat_TS
            end if


        end do
        !$omp end do

        !$omp do private(KOBL, wc, kdir, idgod, gran, normal, s1, s2, s3, qqq1, qqq2, POTOK, dsl, dsc, dsp, c1, c2, c3)
        do i = 1, Num2
            qqq1 = 0.0
            qqq2 = 0.0
            wc = 0.0
            KOBL = 0
            kdir = 0
            idgod = 0
            gran = SS%gl_HP(i)
            normal = SS%gl_Gran_normal(:, gran, step)
            s1 = SS%gl_Gran_neighbour(1, gran)
            s2 = SS%gl_Gran_neighbour(2, gran)

            qqq1(1) = SS%gd(1, s1, step)
            qqq1(2) = SS%gd(3, s1, step)
            qqq1(3) = SS%gd(4, s1, step)
            qqq1(5) = SS%gd(2, s1, step)

            qqq2(1) = SS%gd(1, s2, step)
            qqq2(2) = SS%gd(3, s2, step)
            qqq2(3) = SS%gd(4, s2, step)
            qqq2(5) = SS%gd(2, s2, step)

            ! call chlld(2, normal(1), normal(2), 0.0_8, &
			! 	0.0_8, qqq1, qqq2, dsl, dsp, dsc, POTOK)

            call cgod3d(KOBL, 0, 0, 0, kdir, idgod, &
                normal(1), normal(2), 0.0_8, 1.0_8, &
                wc, qqq1, qqq2, &
                dsl, dsp, dsc, 1.0_8, 1.66666666666666_8, &
                POTOK)

            if (idgod == 2) then
                call chlld(2, normal(1), normal(2), 0.0_8, &
				0.0_8, qqq1, qqq2, dsl, dsp, dsc, POTOK)
            end if

            normal = normal * dsc * koeff_HP
            s1 = SS%gl_all_Gran(1, gran)
            s2 = SS%gl_all_Gran(2, gran)

            c1 = SS%gl_yzel(:, s1, step)
            c2 = SS%gl_yzel(:, s2, step)


            SS%gl_yzel_Vel(:, s1) = SS%gl_yzel_Vel(:, s1) + normal
            SS%gl_yzel_Vel(:, s2) = SS%gl_yzel_Vel(:, s2) + normal

            SS%gl_Point_num(s1) = SS%gl_Point_num(s1) + 1
            SS%gl_Point_num(s2) = SS%gl_Point_num(s2) + 1

            !! Сюда же напишем поверхностное натяжение
            !! Это поверхностое натяжение работает только в случае, если скорость движения узла везда как среднее движения его граней
            ! s1 = SS%gl_all_Gran(1, gran)
            ! s2 = SS%gl_all_Gran(2, gran)

            ! c1 = SS%gl_yzel(:, s1, step)
            ! c2 = SS%gl_yzel(:, s2, step)

            ! if(c1(1) < -80.0) then
			! 	continue
            !     CYCLE
            ! end if

            if(i > 1) then
                gran = SS%gl_HP(i - 1)
                s3 = SS%gl_all_Gran(2, gran)
                c3 = SS%gl_yzel(:, s3, step)
                c3 = (c2 + c3)/2.0
                SS%gl_yzel_Vel(:, s1) = SS%gl_yzel_Vel(:, s1) + 2.0 * (c3 - c1) * SS%par_nat_HP
            end if

            if(i < Num2) then
                gran = SS%gl_HP(i + 1)
                s3 = SS%gl_all_Gran(1, gran)
                c3 = SS%gl_yzel(:, s3, step)
                c3 = (c1 + c3)/2.0
                SS%gl_yzel_Vel(:, s2) = SS%gl_yzel_Vel(:, s2) + 2.0 * (c3 - c2) * SS%par_nat_HP
            end if

        end do
        !$omp end do

        !$omp do private(KOBL, wc, kdir, idgod, gran, normal, s1, s2, s3, qqq1, qqq2, POTOK, dsl, dsc, dsp, c1, c2, c3)
        do i = 1, Num3
            qqq1 = 0.0
            qqq2 = 0.0
            wc = 0.0
            KOBL = 0
            kdir = 0
            idgod = 0
            gran = SS%gl_BS(i)
            normal = SS%gl_Gran_normal(:, gran, step)
            s1 = SS%gl_Gran_neighbour(1, gran)
            s2 = SS%gl_Gran_neighbour(2, gran)

            qqq1(1) = SS%gd(1, s1, step)
            qqq1(2) = SS%gd(3, s1, step)
            qqq1(3) = SS%gd(4, s1, step)
            qqq1(5) = SS%gd(2, s1, step)

            qqq2(1) = SS%gd(1, s2, step)
            qqq2(2) = SS%gd(3, s2, step)
            qqq2(3) = SS%gd(4, s2, step)
            qqq2(5) = SS%gd(2, s2, step)

            ! call chlld(2, normal(1), normal(2), 0.0_8, &
			! 	0.0_8, qqq1, qqq2, dsl, dsp, dsc, POTOK)

            call cgod3d(KOBL, 0, 0, 0, kdir, idgod, &
                normal(1), normal(2), 0.0_8, 1.0_8, &
                wc, qqq1, qqq2, &
                dsl, dsp, dsc, 1.0_8, 1.66666666666666_8, &
                POTOK)

            if (idgod == 2) then
                call chlld(2, normal(1), normal(2), 0.0_8, &
				0.0_8, qqq1, qqq2, dsl, dsp, dsc, POTOK)
            end if

            normal = normal * dsp * koeff_BS
            s1 = SS%gl_all_Gran(1, gran)
            s2 = SS%gl_all_Gran(2, gran)

            ! SS%gl_yzel_Vel(:, s1) = SS%gl_yzel_Vel(:, s1) + normal
            SS%gl_yzel_Vel(:, s2) = SS%gl_yzel_Vel(:, s2) + normal

            ! SS%gl_Point_num(s1) = SS%gl_Point_num(s1) + 1
            SS%gl_Point_num(s2) = SS%gl_Point_num(s2) + 1

            !! Сюда же напишем поверхностное натяжение
            !! Это поверхностое натяжение работает только в случае, если скорость движения узла везда как среднее движения его граней
            s1 = SS%gl_all_Gran(1, gran)
            s2 = SS%gl_all_Gran(2, gran)

            c1 = SS%gl_yzel(:, s1, step)
            c2 = SS%gl_yzel(:, s2, step)

            if(i > 1) then
                gran = SS%gl_BS(i - 1)
                s3 = SS%gl_all_Gran(2, gran)
                c3 = SS%gl_yzel(:, s3, step)
                c3 = (c2 + c3)/2.0
                SS%gl_yzel_Vel(:, s1) = SS%gl_yzel_Vel(:, s1) + 1.0 * (c3 - c1) * SS%par_nat_BS
            end if

            if(i < Num3) then
                gran = SS%gl_BS(i + 1)
                s3 = SS%gl_all_Gran(1, gran)
                c3 = SS%gl_yzel(:, s3, step)
                c3 = (c1 + c3)/2.0
                SS%gl_yzel_Vel(:, s2) = SS%gl_yzel_Vel(:, s2) + 1.0 * (c3 - c2) * SS%par_nat_BS
            end if
        end do
        !$omp end do

        !$omp end parallel

    end subroutine Calc_move_velosity

    subroutine Move_all(SS, step, TT)
        !! Передвигаем узлы сетки
        TYPE (Setka), intent(in out) :: SS
        integer(4), intent(in) :: step
        real(8), intent(in) :: TT

        integer(4) :: step2, N1, N2, j, i, node, del, ii, yz
        real(8) :: the, R_TS, R_HP, R_BS, coord(2), norma, vel(2), the2, coord2(2)

        step2 = mod(step, 2) + 1 

        !! Движение сетки
        ! A - лучи ************************************************************
        !TODO Нужно отдельно сделать движение точки на оси симметрии для A и K лучей
        N2 = size(SS%gl_RAY_A(1, :))
        N1 = size(SS%gl_RAY_A(:, 1))
        do j = 1, N2

            the = (j - 1) * par_pi/2.0/(N2 - 1)

            node = SS%gl_RAY_A(SS%par_n_TS, j)
            coord = SS%gl_yzel(:, node, step)
            norma = norm2(coord)
            vel = SS%gl_yzel_Vel(:, node)
            del = SS%gl_Point_num(node)
            if(del > 1) vel = vel/del
            R_TS = norm2(coord * (1.0 + DOT_PRODUCT(vel * TT, coord/norma)/norma))

            node = SS%gl_RAY_A(SS%par_n_HP, j)
            coord = SS%gl_yzel(:, node, step)
            norma = norm2(coord)
            vel = SS%gl_yzel_Vel(:, node)
            del = SS%gl_Point_num(node)
            if(del > 1) vel = vel/del
            R_HP = norm2(coord * (1.0 + DOT_PRODUCT(vel * TT, coord/norma)/norma))
            !R_HP = norm2(SS%gl_yzel(:, node, step))

            node = SS%gl_RAY_A(SS%par_n_BS, j)
            coord = SS%gl_yzel(:, node, step)
            norma = norm2(coord)
            vel = SS%gl_yzel_Vel(:, node)
            del = SS%gl_Point_num(node)
            if(del > 1) vel = vel/del

            if(j == N2) vel = vel * 5

            R_BS = norm2(coord * (1.0 + DOT_PRODUCT(vel * TT, coord/norma)/norma))
            !R_BS = norm2(SS%gl_yzel(:, node, step))

            do i = 1, N1
                if (i == 1) then
                    CYCLE
                end if

                call Set_Ray_A(SS, i, j, R_TS, R_HP, R_BS, step2)
            end do
        end do

        do j = 1, 1
            the = (j - 1) * par_pi/2.0/(N2 - 1)

            node = SS%gl_RAY_A(SS%par_n_TS, 2)
            coord = SS%gl_yzel(:, node, step)
            norma = norm2(coord)
            R_TS = norma

            node = SS%gl_RAY_A(SS%par_n_HP, 2)
            coord = SS%gl_yzel(:, node, step)
            norma = norm2(coord)
            R_HP = norma

            node = SS%gl_RAY_A(SS%par_n_BS, 2)
            coord = SS%gl_yzel(:, node, step)
            norma = norm2(coord)
            R_BS = norma

            do i = 1, N1
                call Set_Ray_A(SS, i, j, R_TS, R_HP, R_BS, step2)
            end do
        end do

        ! B - лучи ************************************************************
        N2 = size(SS%gl_RAY_B(1, :))
        N1 = size(SS%gl_RAY_B(:, 1))
        do j = 1, N2

            the = par_pi/2.0 + (j) * SS%par_triple_point/(N2)
            the2 = par_pi/2.0 + (j) * SS%par_triple_point_2/(N2)

            node = SS%gl_RAY_B(SS%par_n_TS, j)
            coord = SS%gl_yzel(:, node, step)
            norma = norm2(coord)
            vel = SS%gl_yzel_Vel(:, node)
            del = SS%gl_Point_num(node)
            if(del > 1) vel = vel/del
            R_TS = norm2(coord * (1.0 + DOT_PRODUCT(vel * TT, coord/norma)/norma))
            
            node = SS%gl_RAY_B(SS%par_n_HP, j)
            coord2 = SS%gl_yzel(:, node, step)
            vel = SS%gl_yzel_Vel(:, node)
            del = SS%gl_Point_num(node)
            if(del > 1) vel = vel/del
            coord2 = coord2 + vel * TT
            R_HP = (coord2(2) - R_TS * sin(the))/sin(the2)

            do i = 1, N1
                if (i == 1) then
                    CYCLE
                end if

                call Set_Ray_B(SS, i, j, R_TS, R_HP, step2)
            end do
        end do

        ! C - лучи ************************************************************
        N2 = size(SS%gl_RAY_C(1, :))
        N1 = size(SS%gl_RAY_C(:, 1))
        do j = 1, N2
            do i = 1, N1
                if (i == 1) then
                    CYCLE
                end if

                call Set_Ray_C(SS, i, j, step2)

            end do
        end do

        ! O - лучи ************************************************************
        N2 = size(SS%gl_RAY_O(1, :))
        N1 = size(SS%gl_RAY_O(:, 1))

        do j = 1, N2

            node = SS%gl_RAY_O(1, j)
            coord = SS%gl_yzel(:, node, step)
            vel = SS%gl_yzel_Vel(:, node)
            del = SS%gl_Point_num(node)
            if(del > 1) vel = vel/del
            R_HP = coord(2) + vel(2) * TT

            do i = 1, N1
                call Set_Ray_O(SS, i, j, R_HP, step2)
            end do
        end do

        ! K - лучи ************************************************************
        N2 = size(SS%gl_RAY_K(1, :))
        N1 = size(SS%gl_RAY_K(:, 1))
        do j = 1, N2
            the = par_pi/2.0 + SS%par_triple_point + (N2 - j + 1) * (par_pi/2.0 - SS%par_triple_point)/(N2)  !TODO небезопасно, лучше организовать функцию, которая по номеру луча будет выдавать его угол
            node = SS%gl_RAY_K(SS%par_n_TS, j)
            coord = SS%gl_yzel(:, node, step)
            norma = norm2(coord)
            vel = SS%gl_yzel_Vel(:, node)
            del = SS%gl_Point_num(node)
            if(del > 1) vel = vel/del
            R_TS = norm2(coord * (1.0 + DOT_PRODUCT(vel * TT, coord/norma)/norma))

            do i = 1, N1
                if (i == 1) then
                    CYCLE
                end if

                call Set_Ray_K(SS, i, j, R_TS, step2)
            end do
        end do

        do j = 1, 2
            the = par_pi/2.0 + SS%par_triple_point + (N2 - j + 1) * (par_pi/2.0 - SS%par_triple_point)/(N2)  !TODO небезопасно, лучше организовать функцию, которая по номеру луча будет выдавать его угол
            node = SS%gl_RAY_K(SS%par_n_TS, 3)
            coord = SS%gl_yzel(:, node, step)
            norma = norm2(coord)
            R_TS = norma

            do i = 1, N1
                if (i == 1) then
                    CYCLE
                end if

                call Set_Ray_K(SS, i, j, R_TS, step2)
            end do
        end do

        ! D - лучи ************************************************************
        N2 = size(SS%gl_RAY_D(1, :))
        N1 = size(SS%gl_RAY_D(:, 1))
        do j = 1, N2
            do i = 1, N1
                if(i == 1) then
					CYCLE
                end if

                call Set_Ray_D(SS, i, j, step2)
            end do
        end do

        ! E - лучи ************************************************************
        N2 = size(SS%gl_RAY_E(1, :))
        N1 = size(SS%gl_RAY_E(:, 1))
        do j = 1, N2
            do i = 1, N1
                if(i == 1) then
                    CYCLE
                end if

                if(i == N1) then
                    CYCLE
                end if

                call Set_Ray_E(SS, i, j, step2)
            end do
        end do
    end subroutine Move_all

    subroutine Move_all_play(SS)
        !! Функция для демонстрации движения сетки (произволньное движение поверхностей)
        !! Не-физичная функция, просто для презентации
        !! Изменяет первые координаты на основе вторых
        TYPE (Setka), intent(in out) :: SS

        integer(4) :: N1, N2, j, i, node, del, ii, yz, st
        real(8) :: the, R_TS, R_HP, R_BS, coord(2), norma, vel(2), the2, coord2(2)
        real(8) :: TT


        do st = 1, 200
            TT = (4.0 * 2.0 * par_pi / 200) * (st - 1)
            !! Движение сетки
            ! A - лучи ************************************************************
            !TODO Нужно отдельно сделать движение точки на оси симметрии для A и K лучей
            N2 = size(SS%gl_RAY_A(1, :))
            N1 = size(SS%gl_RAY_A(:, 1))
            do j = 1, N2

                the = (j - 1) * par_pi/2.0/(N2 - 1)

                node = SS%gl_RAY_A(SS%par_n_TS, j)
                coord = SS%gl_yzel(:, node, 2)
                R_TS = norm2(coord * (1.0 + 0.05 * sin(TT) + 0.05 * sin(the * 4.0 + TT) * sin(TT)))

                node = SS%gl_RAY_A(SS%par_n_HP, j)
                coord = SS%gl_yzel(:, node, 2)
                R_HP = norm2(coord * (1.0 + 0.1 * sin(TT/2.0)))

                node = SS%gl_RAY_A(SS%par_n_BS, j)
                coord = SS%gl_yzel(:, node, 2)
                R_BS = norm2(coord * (1.0 + 0.2 * sin(TT/4.0)))
                !R_BS = norm2(SS%gl_yzel(:, node, step))

                do i = 1, N1
                    if (i == 1) then
                        CYCLE
                    end if

                    call Set_Ray_A(SS, i, j, R_TS, R_HP, R_BS, 1)
                end do
            end do

            ! B - лучи ************************************************************
            N2 = size(SS%gl_RAY_B(1, :))
            N1 = size(SS%gl_RAY_B(:, 1))
            do j = 1, N2

                the = par_pi/2.0 + (j) * SS%par_triple_point/(N2)
                the2 = par_pi/2.0 + (j) * SS%par_triple_point_2/(N2)

                node = SS%gl_RAY_B(SS%par_n_TS, j)
                coord = SS%gl_yzel(:, node, 2)
                ! R_TS = norm2(coord * (1.0 + 0.1 * sin(TT)))
                R_TS = norm2(coord * (1.0 + 0.05 * sin(TT) + 0.05 * sin(the * 4.0 + TT) * sin(TT)))
                
                node = SS%gl_RAY_B(SS%par_n_HP, j)
                coord2 = SS%gl_yzel(:, node, 2)
                coord2(2) = coord2(2) * (1.0 + 0.1 * sin(TT/2.0))
                R_HP = (coord2(2) - R_TS * sin(the))/sin(the2)

                do i = 1, N1
                    if (i == 1) then
                        CYCLE
                    end if

                    call Set_Ray_B(SS, i, j, R_TS, R_HP, 1)
                end do
            end do

            ! C - лучи ************************************************************
            N2 = size(SS%gl_RAY_C(1, :))
            N1 = size(SS%gl_RAY_C(:, 1))
            do j = 1, N2
                do i = 1, N1
                    if (i == 1) then
                        CYCLE
                    end if

                    call Set_Ray_C(SS, i, j, 1)

                end do
            end do

            ! O - лучи ************************************************************
            N2 = size(SS%gl_RAY_O(1, :))
            N1 = size(SS%gl_RAY_O(:, 1))

            do j = 1, N2

                yz = SS%gl_RAY_O(1, j)
                coord = SS%gl_yzel(:, yz, 2)
                R_HP = coord(2) * (1.0 + 0.1 * sin(TT/2.0))

                do i = 1, N1
                    call Set_Ray_O(SS, i, j, R_HP, 1)
                end do
            end do

            ! K - лучи ************************************************************
            N2 = size(SS%gl_RAY_K(1, :))
            N1 = size(SS%gl_RAY_K(:, 1))
            do j = 1, N2
                the = par_pi/2.0 + SS%par_triple_point + (N2 - j + 1) * (par_pi/2.0 - SS%par_triple_point)/(N2)  !TODO небезопасно, лучше организовать функцию, которая по номеру луча будет выдавать его угол
                node = SS%gl_RAY_K(SS%par_n_TS, j)
                coord = SS%gl_yzel(:, node, 2)
                ! R_TS = norm2(coord * (1.0 + 0.1 * sin(TT)))
                R_TS = norm2(coord * (1.0 + 0.05 * sin(TT) + 0.05 * sin(the * 4.0 + TT) * sin(TT)))

                do i = 1, N1
                    if (i == 1) then
                        CYCLE
                    end if

                    call Set_Ray_K(SS, i, j, R_TS, 1)
                end do
            end do

            ! D - лучи ************************************************************
            N2 = size(SS%gl_RAY_D(1, :))
            N1 = size(SS%gl_RAY_D(:, 1))
            do j = 1, N2
                do i = 1, N1
                    if(i == 1) then
                        CYCLE
                    end if

                    call Set_Ray_D(SS, i, j, 1)
                end do
            end do

            ! E - лучи ************************************************************
            N2 = size(SS%gl_RAY_E(1, :))
            N1 = size(SS%gl_RAY_E(:, 1))
            do j = 1, N2
                do i = 1, N1
                    if(i == 1) then
                        CYCLE
                    end if

                    if(i == N1) then
                        CYCLE
                    end if

                    call Set_Ray_E(SS, i, j, 1)
                end do
            end do

            ! call Geo_Print_Surface(SS, st)
            call Print_Grans(SS, st)
        end do

    end subroutine Move_all_play

    subroutine Algoritm_Bound_condition(SS)
        !! Задаём граничные условия в сетке (на внутренней сфере для газовой динамики)
        TYPE (Setka), intent(in out) :: SS
        integer(4) :: n1, i, cell, j
        real(8) :: x, y, par(SS%n_par)

        n1 = size(SS%gl_Cell_A(1,:))
        do i = 1, n1
            do j = 1, 2
                cell = SS%gl_Cell_A(j, i)
                x = SS%gl_Cell_Centr(1, cell, 1)
                y = SS%gl_Cell_Centr(2, cell, 1)
                call Inner_Conditions(SS, x, y, par)
                SS%gd(1:5, cell, 1) = par(1:5)
                SS%gd(1:5, cell, 2) = SS%gd(1:5, cell, 1)
            end do
        end do

        n1 = size(SS%gl_Cell_B(1,:))
        do i = 1, n1
            do j = 1, 2
                cell = SS%gl_Cell_B(j, i)
                x = SS%gl_Cell_Centr(1, cell, 1)
                y = SS%gl_Cell_Centr(2, cell, 1)
                call Inner_Conditions(SS, x, y, par)
                SS%gd(1:5, cell, 1) = par(1:5)
                SS%gd(1:5, cell, 2) = SS%gd(1:5, cell, 1)
            end do
        end do

	end subroutine Algoritm_Bound_condition

    subroutine Algoritm_Initial_condition(SS)
        !! Задаём начальные условия в сетке
        TYPE (Setka), intent(in out) :: SS
        integer(4) :: N, i
        real(8) :: x, y, par(SS%n_par)

        N = size(SS%gl_Cell_Centr(1, :, 1))

        do i = 1, N
            x = SS%gl_Cell_Centr(1, i, 1)
            y = SS%gl_Cell_Centr(2, i, 1)

            call Phys_Innitial_Conditions(SS, x, y, par)

            SS%gd(1:5, i, 1) = par(1:5)
            SS%gd(1:5, i, 2) = SS%gd(1:5, i, 1)
        end do

    end subroutine Algoritm_Initial_condition
	
    subroutine Algoritm_ReMove_Surface(SS, SURF)
	    ! Передвигает поверхности сетки согласно поверхностям в SURF
	    TYPE (Setka), intent(in out) :: SS
		TYPE (Surfaces), intent(in out) :: SURF

        real(8) R_TS, the, R_HP, R_BS, yy, xx, the2
		integer(4) :: N1, N2, i, j, yz, ii
		
        ! A - лучи ************************************************************
        N2 = size(SS%gl_RAY_A(1, :))
        N1 = size(SS%gl_RAY_A(:, 1))
        do j = 1, N2

            the = (j - 1) * par_pi/2.0/(N2 - 1)

            R_TS = SUR_GET_TS(SURF, the)
            R_HP = SUR_GET_HP_alp(SURF, the)
            R_BS = SUR_GET_BS_alp(SURF, the)

            do i = 1, N1
                if (i == 1) then
                    CYCLE
                end if

                call Set_Ray_A(SS, i, j, R_TS, R_HP, R_BS, 1)
            end do
        end do

        ! B - лучи ************************************************************
        N2 = size(SS%gl_RAY_B(1, :))
        N1 = size(SS%gl_RAY_B(:, 1))
		
		do ii = 1, 5  !! Этот цикл нужен, т.к. HP ищется по старой x - координате, нужно делать итерации
            do j = 1, N2

                the = par_pi/2.0 + (j) * SS%par_triple_point/(N2)
		        the2 = par_pi/2.0 + (j) * SS%par_triple_point_2/(N2)

                R_TS = SUR_GET_TS(SURF, the)
                yz = SS%gl_RAY_B(SS%par_n_HP, j)
                xx = SS%gl_yzel(1, yz, 1)
                yy = SUR_GET_HP_x(SURF, xx)
                ! Здесь для HP нужна её y - координата (т.е. yy - это просто высота)
                R_HP = (yy - R_TS * sin(the))/sin(the2)

                do i = 1, N1
                    if (i == 1) then
                        CYCLE
                    end if

                    call Set_Ray_B(SS, i, j, R_TS, R_HP, 1)
                end do
		    end do
		end do
		
        ! C - лучи ************************************************************
        N2 = size(SS%gl_RAY_C(1, :))
        N1 = size(SS%gl_RAY_C(:, 1))
        do j = 1, N2
            do i = 1, N1
                if (i == 1) then
                    CYCLE
                end if

                call Set_Ray_C(SS, i, j, 1)

            end do
        end do

        ! O - лучи ************************************************************
        N2 = size(SS%gl_RAY_O(1, :))
        N1 = size(SS%gl_RAY_O(:, 1))

        do ii = 1, 5  !! Этот цикл нужен, т.к. HP ищется по старой x - координате, нужно делать итерации
            do j = 1, N2

                yz = SS%gl_RAY_O(1, j)
                xx = SS%gl_yzel(1, yz, 1)
                R_HP = SUR_GET_HP_x(SURF, xx)

                do i = 1, N1

                    call Set_Ray_O(SS, i, j, R_HP, 1)

                end do
            end do
        end do

        ! K - лучи ************************************************************
        N2 = size(SS%gl_RAY_K(1, :))
        N1 = size(SS%gl_RAY_K(:, 1))
        do j = 1, N2
            the = par_pi/2.0 + SS%par_triple_point + (N2 - j + 1) * (par_pi/2.0 - SS%par_triple_point)/(N2)  !TODO небезопасно, лучше организовать функцию, которая по номеру луча будет выдавать его угол
            R_TS = SUR_GET_TS(SURF, the)
            do i = 1, N1
                if (i == 1) then
                    CYCLE
                end if

                call Set_Ray_K(SS, i, j, R_TS, 1)

            end do
        end do

        ! D - лучи ************************************************************
        N2 = size(SS%gl_RAY_D(1, :))
        N1 = size(SS%gl_RAY_D(:, 1))
        do j = 1, N2
            do i = 1, N1
                if(i == 1) then
					CYCLE
                end if

                call Set_Ray_D(SS, i, j, 1)
            end do
        end do

        ! E - лучи ************************************************************
        N2 = size(SS%gl_RAY_E(1, :))
        N1 = size(SS%gl_RAY_E(:, 1))
        do j = 1, N2
            do i = 1, N1
                if(i == 1) then
                    CYCLE
                end if

                if(i == N1) then
                    CYCLE
                end if

                call Set_Ray_E(SS, i, j, 1)
            end do
        end do

        SS%gl_yzel(:, :, 2) = SS%gl_yzel(:, :, 1)

        call Culc_Cell_Centr(SS, 1)
        call Culc_Cell_Centr(SS, 2)  
        call Geo_Culc_normal(SS, 1) 
        call Geo_Culc_normal(SS, 2)
        

        call Belong_Init(SS)

        call Geo_Culc_length_area(SS, 1)
        call Geo_Culc_length_area(SS, 2)

        call Geo_Find_Surface(SS)  ! Находим поверхности, которые выделяем

        call Geo_culc_TVD_sosed(SS)
        call Geo_Culc_zone(SS)
		
	end subroutine Algoritm_ReMove_Surface

    subroutine Algoritm_Reinterpol(SS, XX)
        !! Переинтерполирует значения атомарного водорода из сетки интерполяции в центры ячеек новой сетки
	    ! Передвигает поверхности сетки согласно поверхностям в SURF
	    TYPE (Setka), intent(in out) :: SS
		TYPE (Inter_Setka), intent(in out) :: XX
        integer(4) :: N, i, num
        real(8) :: center(2)
        real(8) :: parH(5, 4)
        real(8) :: par(5)

        if(size(SS%hydrogen(:, 1, 1, 2)) /= size(parH(:, 1))) STOP "ERROR 1 size Algoritm_Reinterpol 5y65gw4fervsgf "
        if(size(SS%hydrogen(1, :, 1, 2)) /= size(parH(1, :))) STOP "ERROR 2 size Algoritm_Reinterpol 98y8t4uhvtiewvgtssfvgs "
    

        N = size(SS%gl_Cell_Centr(1, :, 1))
        num = 1

        do i = 1, N
            center = SS%gl_Cell_Centr(:, i, 1)
            call Int_Get_Parameter(XX, center(1), center(2), num, PAR_hydrogen = parH, PAR_gd = par)
            SS%hydrogen(:, :, i, 1) = parH
            SS%hydrogen(:, :, i, 2) = parH
            if(par(1) <= 0.0) par(1) = 0.0000001
            if(par(2) <= 0.0) par(2) = 0.000001
            SS%gd(:, i, 1) = par
            SS%gd(:, i, 2) = par
        end do
    end subroutine Algoritm_Reinterpol

end module Algoritm