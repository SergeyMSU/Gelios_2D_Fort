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
    USE PUI
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

    subroutine Gas_dynamic_algoritm2(SS)
        ! Это временно, чтобы не удалять предыдущую функцию
        TYPE (Setka), intent(in out) :: SS
        integer(4) :: num, i, i_max
        real(8) :: par(5), parH(5, 4), r(2)

        print*, "A"
        call Read_setka_bin(gl_S3, "DD018")   ! ДЛЯ ВОДОРОДА
        print*, "B"
        call Int_Init(gl_I1, gl_S3)
        print*, "C"

        
        print*, "D"
        call Int_Print_Cell(gl_I1)
        print*, "E"

        call Read_setka_bin(SS, "CC019")      ! ОСНОВНАЯ СЕТКА
        print*, "F"
        call Geo_Set_sxem(SS)
        print*, "G"

        ! call Print_hydrogen(gl_S3)
        call Dell_Setka(gl_S3)
        

        !call Algoritm_Reinterpol(SS, gl_S2)

        if(SS%init_geo == .False.) STOP "Gas_dynamic_algoritm  error init_geo erty67543"  


        ! Параметры МАЛАМЫ
        ! SS%par_n_H_LISM = 3.0
        ! SS%par_Velosity_inf = -2.54385
        ! SS%par_Kn = 34.516
        ! SS%par_nu_ph = 11.7745 
        ! SS%par_E_ph = 0.0725533
        ! SS%par_R0 = 0.20445
        ! SS%par_chi = 43.1202
        ! SS%par_rho_e = 128.667
        ! SS%par_Max_e = 6.21523

        ! Параметры Модели, которые сейчас строим
        SS%par_n_H_LISM = 3.0
        SS%par_Velosity_inf = -2.54278_8
        SS%par_Kn = 50.3858
        SS%par_nu_ph = 12.0969 
        SS%par_E_ph = 0.10878
        SS%par_R0 = 0.198956
        SS%par_chi = 41.0391
        SS%par_rho_e = 150.0
        SS%par_Max_e = 5.91662
        SS%par_a_2 = 0.130735_8
        SS%culc_pui = .True.


        !call Algoritm_Initial_condition(SS)  ! Зададим начальные условия для всех ячеек сетки
        !call Geo_get_request(SS)
        call Algoritm_Bound_condition(SS)    ! Зададим граничные условия на внутренней сфере

        ! call Int_Init(gl_S2, SS)
        ! call Int_Print_Cell(gl_S2)
        
        !call Int_Get_Parameter(gl_S2, 260.0_8, 19.5_8, num, PAR_gd = par, PAR_hydrogen = parH)

        ! call Calc_move_velosity(SS, 1)
        ! call Move_all(SS, 2, 1.0_8)

        ! do i = 1, size(SS%gl_Cell_Centr(1, :, 1))  !! УДАЛИТЬ
        !     if(SS%gl_Cell_Centr(1, i, 1) > 120.0) then
        !         SS%gd(1, i, 1) = 1.0
        !         SS%gd(2, i, 1) = 1.0
        !         SS%gd(3, i, 1) = SS%par_Velosity_inf
        !         SS%gd(4, i, 1) = 0.0
        !         SS%gd(:, i, 2) = SS%gd(:, i, 1)
        !     end if
        ! end do

        call Print_GD(SS)
        call Geo_Print_Surface(SS)
        call Print_Grans(SS)
        call Print_GD_1D(SS)
        call Algoritm_Reinterpol(SS, gl_I1, gd_ = .False.)

        !call Print_hydrogen(SS)

        ALLOCATE(SS%gl_Gran_POTOK(size(SS%gl_all_Gran(1, :))))

        SS%gl_Gran_POTOK = -100000.0

        SS%par_nat_TS = 0.04_8 ! 0.03_8
        SS%par_nat_HP = 0.003_8  ! 0.04  0.06
        SS%par_nat_BS = 0.006_8 !0.004_8

        SS%par_koeff_HP = 0.03_8  

    

        print*, "Proverim parametry"
        print*, SS%par_n_H_LISM
        print*, SS%par_Kn
        print*, SS%par_a_2
        print*, SS%par_nu_ph
        print*, "________________________"


        print*, "H"
        i_max = 700!200!350   100 - 7 минут
        do i = 1, i_max
            !SS%par_kk2 = SS%par_kk2 + 0.2/300
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

            call Algoritm_Reinterpol(SS, gl_I1, gd_ = .False.)
        end do



        call Print_GD(SS)
        call Geo_Print_Surface(SS, 20)
        call Save_setka_bin(SS, "CC020")
        call Print_Grans(SS)
        ! call Print_Cell_Centr(SS)
        call Print_GD_1D(SS)
        ! call Print_TVD_Sosed(SS)

        ! call Print_hydrogen(SS)
        ! call Print_hydrogen_1D(SS)

        !pause
    end subroutine Gas_dynamic_algoritm2

    subroutine Perestroika_algoritm(SS)
        ! Пересетройка основной сетки
        TYPE (Setka), intent(in out) :: SS

        call Read_setka_bin(SS, "A0047")      ! ОСНОВНАЯ СЕТКА
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
        gl_S3%par_n_TS =  33! 26                    ! Количество точек до TS (TS включается)
        gl_S3%par_n_HP =  63! 40                 ! Количество точек HP (HP включается)  всё от 0 считается
        gl_S3%par_n_BS =  89! 60! 5                 ! Количество точек BS (BS включается)
        gl_S3%par_n_END = 101! 72! 6                ! Количество точек до конца сетки (конец включается)
        gl_S3%par_n_IA =  20! 12                   ! Количество точек, которые входят во внутреннюю область
        gl_S3%par_n_IB =  22! 14                   ! Количество точек, которые входят во внутреннюю область (с зазором)
        gl_S3%par_kk2 = 2.08
        call Init_Setka(gl_S3)
		print*, "A4"
        call Build_Setka_start(gl_S3)
		print*, "A5"
	    call Algoritm_ReMove_Surface(gl_S3, gl_surf1)
		print*, "A6"
        call Algoritm_Reinterpol(gl_S3, gl_S2)

        call Print_Cell(gl_S3)
        call Print_GD(gl_S3)
        call Geo_Print_Surface(gl_S3)
        call Print_Grans(gl_S3)
        call Print_GD_1D(gl_S3)

        call Save_setka_bin(gl_S3, "A0048")

    end subroutine Perestroika_algoritm

    subroutine MK_algoritm(SS)
        TYPE (Setka), intent(in out) :: SS
        integer(4) :: i, j, cell

        ! Сетка водорода нужно только в случае использования ПИКАПОВ   culc_pui == True
        call Read_setka_bin(gl_S4, "DD017")   ! ДЛЯ ВОДОРОДА (Предыдущий расчёт)

        call Read_setka_bin(SS, "CC018")      ! ОСНОВНАЯ СЕТКА

        ! call Print_GD(SS)
        !call Geo_Print_Surface(SS)
        !call Print_Grans(SS)
        ! call Print_GD_1D(SS)
        !pause

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
        gl_S3%par_n_END = 68! 72! 6                ! Количество точек до конца сетки (конец включается)
        gl_S3%par_n_IA =  12! 12                   ! Количество точек, которые входят во внутреннюю область
        gl_S3%par_n_IB =  14! 14                   ! Количество точек, которые входят во внутреннюю область (с зазором)
        gl_S3%par_kk13 =  1.6!
        gl_S3%par_kk14 =  0.8!
        gl_S3%par_kk113 = 1.2

        !! Физические параметры модели
        gl_S3%par_n_H_LISM = 3.0
        gl_S3%par_Velosity_inf = -2.54278_8
        gl_S3%par_Kn = 50.3858
        gl_S3%par_nu_ph = 12.0969 
        gl_S3%par_E_ph = 0.10878
        gl_S3%par_R0 = 0.198956
        gl_S3%par_chi = 41.0391
        gl_S3%par_rho_e = 150.0
        gl_S3%par_Max_e = 5.91662
        gl_S3%par_a_2 = 0.130735_8
        gl_S3%par_poglosh = 0.389274
        gl_S3%culc_pui = .True.

        call Init_Setka(gl_S3)
		print*, "A4"
        call Build_Setka_start(gl_S3)

        if(gl_S3%culc_pui == .True.) then
            call PUI_SET(gl_S3)
        end if

		print*, "A5"
	    call Algoritm_ReMove_Surface(gl_S3, gl_surf1)
		print*, "A6"

        ! call Print_Cell(gl_S3)
        ! call Geo_Print_Surface(gl_S3)
        ! call Print_Grans(gl_S3)
        ! return

        !print*, gl_S2%gd(:, 1)
        !print*, gl_S2%gd(:, 2)
        !print*, gl_S2%gd(:, 3)

        print*, "A7"

        call Algoritm_Reinterpol(gl_S3, gl_S2)
        
		print*, "A8"

        if(gl_S3%culc_pui == .True.) then
            call Algoritm_Reinterpol_S_pui(gl_S3, gl_S4)  ! Берём S+ S- с предыдущего расчёта
            call Culc_f_pui(gl_S3, gl_S2)
            call PUI_F_integr_Set(gl_S3)
            call PUI_Culc_h0(gl_S3)
            call PUI_F_integr_Culc(gl_S3)
            call PUI_n_T_culc(gl_S3)
        end if
        call Dell_Setka(gl_S4)


        print*, "Proverim parametry"
        print*, gl_S3%par_n_H_LISM
        print*, gl_S3%par_Kn
        print*, gl_S3%par_a_2
        print*, gl_S3%par_nu_ph
        print*, "________________________"

        call Print_GD_PUI(gl_S3)
        !call PUI_print_pui(gl_S3, -56.2_8, 12.96_8)
        !call PUI_print_pui(gl_S3, -74.5_8, 21.35_8)
        !return 

        call M_K_start(gl_S3, gl_S2)

        call Print_hydrogen(gl_S3)
        call Print_hydrogen_1D(gl_S3)
        call Calc_Pogloshenie(gl_S3)

        call Print_GD_PUI(gl_S3)

        if(gl_S3%culc_pui == .True.) then
            call Culc_f_pui(gl_S3, gl_S2)
            call PUI_print_pui(gl_S3, 15.0_8, 0.0001_8)
            call PUI_print_pui(gl_S3, -15.0_8, 0.0001_8)
            call PUI_print_pui(gl_S3, 0.0_8, 15.0001_8)
            call PUI_print_pui(gl_S3, 25.0_8, 0.0001_8)
            call PUI_print_pui(gl_S3, -50.0_8, 0.0001_8)
        end if


        call Save_setka_bin(gl_S3, "DD018")
        call Print_PUI_1D(gl_S3)

        print*, "END"

    end subroutine MK_algoritm

    subroutine Print_PUI_algoritm(SS)
        TYPE (Setka), intent(in out) :: SS

         ! Сетка водорода нужно только в случае использования ПИКАПОВ   culc_pui == True
        call Read_setka_bin(gl_S4, "DD012")   ! ДЛЯ ВОДОРОДА (Предыдущий расчёт)

        call Read_setka_bin(SS, "CC012")      ! ОСНОВНАЯ СЕТКА

		print*, "A1"
        call SUR_init(gl_surf1, SS)
		print*, "A2"
        call Int_Init(gl_S2, SS)
		print*, "A3"

        call Culc_f_pui(gl_S4, gl_S2)
        !call PUI_F_integr_Set(gl_S4)
        !call PUI_Culc_h0(gl_S4)
        !call PUI_F_integr_Culc(gl_S4)
        call PUI_n_T_culc(gl_S4)

        call Print_PUI_1D(gl_S4)

        print*, "END"

    end subroutine Print_PUI_algoritm

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
        logical :: tvd1, tvd2, null_un
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


            !$omp do private(null_un, KOBL, kdir, idgod, sosed_center, phi3, gran_center, Vr, Vphi, phi1, phi2, r, par1_TVD, wc, &
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
                    null_un = .False.
					POTOK2 = 0.0
                    KOBL = 0
                    kdir = 0
                    idgod = 0
                    gran = SS%gl_Cell_gran(gr, cell)
                    sosed = SS%gl_Cell_neighbour(gr, cell)
                    if(sosed == 0) CYCLE
                    gran_center = SS%gl_Gran_Center(:, gran, now)

                    normal = SS%gl_Gran_normal(:, gran, now)
                    Sqv = SS%gl_Gran_length(gran, now)
                    lenght = SS%gl_Cell_gran_dist(gr, cell, now)

                    if(SS%gl_Gran_neighbour(1, gran) /= cell) normal = -normal

                    ! Нужно вычислить скорость движения грани
                    wc = DOT_PRODUCT((SS%gl_Gran_Center(:, gran, now2) -  gran_center)/TT, normal)

                    call Get_gran_parameter(SS, gran, cell, par1_TVD, par2, now)  !! Основная функция получения параметров

                    ! if(center(1) > 25 .and. center(1) < 30 .and. center(2) < 4 .and. &
                    ! ( SS%gl_Gran_neighbour(2, gran) == -4 .or. SS%gl_Gran_neighbour(1, gran) == -4) ) then
                    !     print*, "EEEEEE"
                    !     print*, par1_TVD
                    !     print*, "EEEEEE"
                    !     print*, par2
                    !     print*, "EEEEEE"
                    !     print*, par1
                    !     print*, "EEEEEE"
                    !     pause
                    ! end if

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

                        if(SS%gl_Gran_type(gran) == 2 .and. gran_center(2) < 8) then  ! Контакт  5
                            wc = 0.0
                            qqq2 = qqq1
                            wc = DOT_PRODUCT(qqq1(2:3), normal)
                            qqq2(2:3) = qqq2(2:3) - 2.0 * wc * normal
                            wc = 0.0

                            call cgod3d(KOBL, 0, 0, 0, kdir, idgod, &
                            normal(1), normal(2), 0.0_8, 1.0_8, &
                            wc, qqq1(1:8), qqq2(1:8), &
                            dsl, dsp, dsc, 1.0_8, 1.66666666666666_8, &
                            POTOK2)

                            POTOK2(1) = 0.0
                            POTOK2(4) = 0.0
                            POTOK2(5) = 0.0
                            POTOK2(6) = 0.0
                            POTOK2(7) = 0.0
                            POTOK2(8) = 0.0
                            POTOK2(9) = 0.0

                            ! print*, "____________________"
                            ! print*, qqq1
                            ! print*, "____________________"
                            ! print*, qqq2
                            ! print*, "____________________"
                            ! print*, POTOK2
                            ! print*, "____________________"
                            ! pause

                            if (idgod == 2) STOP "ERROR okrfi9uhebrtomeevjoerhbbvecwwvertbhyrvgf"

                        else
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
                        end if

                    else
                        call chlld_Q(SS%gl_Gran_shem(gran), normal(1), normal(2), 0.0_8, &
                        wc, qqq1, qqq2, dsl, dsp, dsc, POTOK2, .False.)
                    end if

                    ! if(SS%gl_Cell_neighbour(gr, cell) == -4 .and. gran_center(1) > 35 .and. gran_center(1) < 60) then
                    !     print*, POTOK2
                    !     print*, "_______________"
                    !     pause
                    ! end if

                    ! if(gran == 77) then
                    !     print*, qqq1
                    !     print*, "__________________"
                    !     print*, qqq2
                    !     print*, "__________________"
                    !     print*, SS%gl_Gran_shem(gran), wc
                    !     print*, "__________________"
                    !     print*, POTOK2
                    !     print*, "__________________"
                    !     print*, dsl, dsp, dsc
                    !     print*, "__________________"
                    !     print*, normal
                    !     print*, "__________________"
                    !     pause
                    ! end if

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


                     !if(SS%gl_Gran_POTOK(gran) < -99000) then
                     !    SS%gl_Gran_POTOK(gran) = POTOK2(1)
                     !else if(dabs(dabs(POTOK2(1)) - dabs(SS%gl_Gran_POTOK(gran))) > 0.00001) then
                     !    print*, "________________"
                     !    print*, POTOK2(1)
                     !    print*, "________________"
                     !    print*, SS%gl_Gran_POTOK(gran)
                     !    print*, "________________"
                     !    print*, gran_center
                     !    print*, "________________"
                     !    print*, gran
                     !    print*, "________________"
						               !
                     !    pause
                     !end if

                end do

                Vol = SS%gl_Cell_square(cell, now)
                Vol2 = SS%gl_Cell_square(cell, now2)

                call Calc_sourse_MF(SS, cell, source, now, use_koeff_ = .True.)

                !if(center(1) < -215) source = 0.0  !! УБРАТЬ

                if(ieee_is_nan(source(2)) .or. ieee_is_nan(source(1)) .or. ieee_is_nan(source(4))) then
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

                ! if(SS%gl_Cell_type(cell) == 'A' .and. SS%gl_Cell_number(2, cell) == 1) then !! УБРАТЬ
                !     source(3) = 0.0
                ! end if

                ! Законы сохранения в ячейке
                ro2 = ro * Vol/Vol2 - TT * (POTOK(1) / Vol2 + ro * v/center(2) - source(1))
                if(ro2 <= 0.0) then
                    print*, "Ro < 0", ro2, ro, TT, Vol, cell
                    print*, "centr = ", center
                    ! print*, "_____________"
                    ! print*, "POTOK = ", POTOK
                    ! print*, "_____________"
                    ! print*, par1, "||||||||| ", par2
                    ! print*, "_____________"
                    ! stop
                    ro2 = 0.03
                end if
                Q2 = Q * Vol/Vol2 - TT * (POTOK(5) / Vol2 + Q * v/center(2) - (Q/ro) * source(1))
                u2 = (ro * u * Vol/Vol2 - TT * ( POTOK(3) / Vol2 + ro * v * u/center(2) - source(2) )) / ro2
                v2 = (ro * v * Vol/Vol2 - TT * ( POTOK(4) / Vol2 + ro * v * v/center(2) - source(3) )) / ro2

                pp = v * (SS%par_ggg * p / (SS%par_ggg - 1.0) + ro * (u * u + v * v) * 0.5) / center(2)
                p2 = ((  ( p / (SS%par_ggg - 1.0) + 0.5 * ro * (u**2 + v**2))  * Vol/Vol2   &
                        - TT * (POTOK(2)/ Vol2 + pp - source(4)) ) - 0.5 * ro2 * (u2**2 + v2**2) ) * (SS%par_ggg - 1.0)

                if(p2 <= 0.0) then
                    p2 = 0.1! 0.000001   !! УУУБРАТЬ
                end if

                SS%gd(1, cell, now2) = ro2
                SS%gd(2, cell, now2) = p2
                SS%gd(3, cell, now2) = u2
                SS%gd(4, cell, now2) = v2
                SS%gd(5, cell, now2) = Q2

			end do
			!$omp end do
            !$omp end parallel

     !         print*, "END"
			!pause
        end do


        

        ! print*, "ALL_TIME = ", ALL_TIME

    end subroutine Start_GD_algoritm

    subroutine Calc_move_velosity(SS, step)
        !! Вычисляем скорости движения граней
        TYPE (Setka), intent(in out) :: SS
        integer(4), intent(in) :: step
        integer(4) :: Num, i, s1, s2, gran, Num2, Num3, s3
        real(8) :: normal(2), qqq1(8), qqq2(8), POTOK(8), nat_HP
        real(8) :: dsl, dsc, dsp
        real(8) :: c1(2), c2(2), c3(2), wc
        integer(4) :: kdir, idgod, KOBL


        ! SS%par_koeff_TS = 0.002
        ! SS%par_koeff_HP = 0.1
        ! SS%par_koeff_BS = 0.1

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

            normal = normal * dsl * SS%par_koeff_TS
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

        !$omp do private(nat_HP, KOBL, wc, kdir, idgod, gran, normal, s1, s2, s3, qqq1, qqq2, POTOK, dsl, dsc, dsp, c1, c2, c3)
        do i = 1, Num2
            nat_HP = SS%par_nat_HP
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

            !! Попробуем убрать тангенциальные скорости (хотим убрать неустойчивость)
            ! qqq1(2:3) = DOT_PRODUCT(qqq1(2:3), normal) * normal
            ! qqq2(2:3) = DOT_PRODUCT(qqq2(2:3), normal) * normal

            call cgod3d(KOBL, 0, 0, 0, kdir, idgod, &
                normal(1), normal(2), 0.0_8, 1.0_8, &
                wc, qqq1, qqq2, &
                dsl, dsp, dsc, 1.0_8, 1.66666666666666_8, &
                POTOK)

            if (idgod == 2) then
                call chlld(2, normal(1), normal(2), 0.0_8, &
				0.0_8, qqq1, qqq2, dsl, dsp, dsc, POTOK)
            end if

            normal = normal * dsc * SS%par_koeff_HP
            s1 = SS%gl_all_Gran(1, gran)
            s2 = SS%gl_all_Gran(2, gran)

            c1 = SS%gl_yzel(:, s1, step)
            c2 = SS%gl_yzel(:, s2, step)


            SS%gl_yzel_Vel(:, s1) = SS%gl_yzel_Vel(:, s1) + normal
            SS%gl_yzel_Vel(:, s2) = SS%gl_yzel_Vel(:, s2) + normal

            SS%gl_Point_num(s1) = SS%gl_Point_num(s1) + 1
            SS%gl_Point_num(s2) = SS%gl_Point_num(s2) + 1

            !! Сюда же напишем поверхностное натяжение
            !! Это поверхностое натяжение работает только в случае, если скорость движения узла везде как среднее движения его граней
            ! s1 = SS%gl_all_Gran(1, gran)
            ! s2 = SS%gl_all_Gran(2, gran)

            ! c1 = SS%gl_yzel(:, s1, step)
            ! c2 = SS%gl_yzel(:, s2, step)

            ! if(c1(1) < -80.0) then
			! 	continue
            !     CYCLE
            ! end if

            !if(c1(1) > -50) nat_HP = nat_HP * 10

            !nat_HP = GET_HP_nat(polar_angle(c1(2), c1(1)))

            if(i > 1) then
                gran = SS%gl_HP(i - 1)
                s3 = SS%gl_all_Gran(2, gran)
                c3 = SS%gl_yzel(:, s3, step)
                c3 = (c2 + c3)/2.0
                SS%gl_yzel_Vel(:, s1) = SS%gl_yzel_Vel(:, s1) + 2.0 * (c3 - c1) * nat_HP
            end if

            if(i < Num2) then
                gran = SS%gl_HP(i + 1)
                s3 = SS%gl_all_Gran(1, gran)
                c3 = SS%gl_yzel(:, s3, step)
                c3 = (c1 + c3)/2.0
                SS%gl_yzel_Vel(:, s2) = SS%gl_yzel_Vel(:, s2) + 2.0 * (c3 - c2) * nat_HP
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

            normal = normal * dsp * SS%par_koeff_BS
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

        do j = 2, 2
            the = (j - 1) * par_pi/2.0/(N2 - 1)
            the2 = (j) * par_pi/2.0/(N2 - 1)

            node = SS%gl_RAY_A(SS%par_n_TS, 2) !! --------------------------------------------------
            coord = SS%gl_yzel(:, node, step2)
            norma = norm2(coord)
            R_TS = norma


            node = SS%gl_RAY_A(SS%par_n_HP, 2)
            coord2 = SS%gl_yzel(:, node, step2)

            node = SS%gl_RAY_A(SS%par_n_HP, 3)
            coord = SS%gl_yzel(:, node, step2)
            if(coord2(1) < coord(1)) then
                norma = norm2(coord)
                R_HP = norma * cos(the2)/cos(the)
            else
                norma = norm2(coord2)
                R_HP = norma
            end if

            node = SS%gl_RAY_A(SS%par_n_BS, 2)
            coord2 = SS%gl_yzel(:, node, step2)
            node = SS%gl_RAY_A(SS%par_n_BS, 3)
            coord = SS%gl_yzel(:, node, step2)
            if(coord2(1) < coord(1)) then
                norma = norm2(coord)
                R_BS = norma * cos(the2)/cos(the)
            else
                norma = norm2(coord2)
                R_BS = norma
            end if

            do i = 1, N1
                call Set_Ray_A(SS, i, j, R_TS, R_HP, R_BS, step2)
            end do
        end do

        do j = 1, 1
            the = (j - 1) * par_pi/2.0/(N2 - 1)
            the2  = (j) * par_pi/2.0/(N2 - 1)

            node = SS%gl_RAY_A(SS%par_n_TS, 2)
            coord = SS%gl_yzel(:, node, step2)
            norma = norm2(coord)
            R_TS = norma

            node = SS%gl_RAY_A(SS%par_n_HP, 2)
            coord = SS%gl_yzel(:, node, step2)
            norma = norm2(coord)
            R_HP = norma * cos(the2)/cos(the)

            node = SS%gl_RAY_A(SS%par_n_BS, 2)
            coord = SS%gl_yzel(:, node, step2)
            norma = norm2(coord)
            R_BS = norma * cos(the2)

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
            !print*, R_TS, the
            !pause
			
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

    subroutine Algoritm_Reinterpol(SS, XX, gd_)
        !! Переинтерполирует значения атомарного водорода из сетки интерполяции в центры ячеек новой сетки
        !! Нужно ли переинтерполировать газовую динамику, или только водород
	    ! Передвигает поверхности сетки согласно поверхностям в SURF
	    TYPE (Setka), intent(in out) :: SS
		TYPE (Inter_Setka), intent(in out) :: XX
        logical, intent(in), OPTIONAL :: gd_
        integer(4) :: N, i, num, j
        real(8) :: center(2)
        real(8) :: parH(5, 4)
        real(8) :: par(5)
        real(8) :: source(4)
        logical :: gd

        gd = .True.
        if(present(gd_)) gd = gd_

        if(size(SS%hydrogen(:, 1, 1, 2)) /= size(parH(:, 1))) STOP "ERROR 1 size Algoritm_Reinterpol 5y65gw4fervsgf "
        if(size(SS%hydrogen(1, :, 1, 2)) /= size(parH(1, :))) STOP "ERROR 2 size Algoritm_Reinterpol 98y8t4uhvtiewvgtssfvgs "
    

        N = size(SS%gl_Cell_Centr(1, :, 1))
        num = 1

        do i = 1, N
            center = SS%gl_Cell_Centr(:, i, 1)
            call Int_Get_Parameter(XX, center(1), center(2), num, PAR_hydrogen = parH, PAR_gd = par, PAR_atom_source = source)
            do j = 1, 4
                if(parH(1, j) <= 0.0) parH(1, j) = 0.0000001
                if(parH(2, j) <= 0.0) parH(2, j) = 0.0000001
            end do

            if(parH(1, 4) > 3.0) then
                print*, "ERROR  roh4"
                print*, center
                pause
            end if


            SS%hydrogen(:, :, i, 1) = parH
            SS%hydrogen(:, :, i, 2) = parH

            if(gd == .True.) then
                if(par(1) <= 0.0) par(1) = 0.0000001
                if(par(2) <= 0.0) par(2) = 0.000001
                SS%gd(:, i, 1) = par
                SS%gd(:, i, 2) = par
            end if

            !if(center(1) > 110) source = 0.0  !! УБРАТЬ
            SS%atom_source(1:4, i) = source
        end do
    end subroutine Algoritm_Reinterpol

    subroutine Algoritm_Reinterpol_S_pui(SS, XX)
        !! Переинтерполирует значения атомарного водорода из сетки интерполяции в центры ячеек новой сетки
        !! Нужно ли переинтерполировать газовую динамику, или только водород
	    ! Передвигает поверхности сетки согласно поверхностям в SURF
	    TYPE (Setka), intent(in out) :: SS
	    TYPE (Setka), intent(in) :: XX
        
        integer(4) :: N, i, num, j, cell, num2
        real(8) :: r(2)


        N = size(SS%gl_Cell_Centr(1, :, 1))
        num = 1
        cell = 1

        SS%pui_Sm = 0.0
        SS%pui_Sp = 0.0

        if(size(SS%pui_Sm(:, 1)) /= size(XX%pui_Sm(:, 1))) then
            print*, "ERROR Algoritm_Reinterpol_S_pui 9i847yt65vgwohe9p4ivu5t098ytc03r"
            print*, size(SS%pui_Sm(:, 1)), size(XX%pui_Sm(:, 1))
            pause
            STOP
        end if


        do i = 1, N
            if(SS%gl_all_Cell_zone(i) >= 3) CYCLE
            r = SS%gl_Cell_Centr(:, i, 1)
            call Geo_Find_Cell(XX, r(1), r(2), cell)
            num = XX%f_pui_num2(cell)
            if(num <= 0) CYCLE
            num2 = SS%f_pui_num2(i)
            SS%pui_Sm(:, num2) = XX%pui_Sm(:, num)
            SS%pui_Sp(:, num2) = XX%pui_Sp(:, num)
        end do
    end subroutine Algoritm_Reinterpol_S_pui

end module Algoritm