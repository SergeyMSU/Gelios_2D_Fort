module Algoritm
    USE STORAGE
    USE GEOMETRY
    USE Phys_parameter
	USE SURFACE
    USE Solvers
    implicit none 

    contains

    subroutine Gas_dynamic_algoritm(SS)
        TYPE (Setka), intent(in out) :: SS

        if(SS%init_geo == .False.) then
            STOP "Gas_dynamic_algoritm  error init_geo erty67543"  
        end if

        call Algoritm_Initial_condition(SS)  ! Зададим начальные условия для всех ячеек сетки
        call Algoritm_Bound_condition(SS)    ! Зададим граничные условия на внутренней сфере

        call Start_GD_algoritm(SS, 1000)
    end subroutine Gas_dynamic_algoritm

    subroutine Start_GD_algoritm(SS, all_step_)
        TYPE (Setka), intent(in out) :: SS
        integer(4), intent(in) :: all_step_

        integer(4) :: now, now2, step, cell, Ncell, gr, gran, sosed, all_step
        real(8) :: time, TT, center(2), par1(5), par2(5), ro, p, u, v, Q, normal(2), Sqv, Vol
        real(8) :: qqq1(9), qqq2(9), POTOK2(9), POTOK(5)
        real(8) :: dsl, dsc, dsp, loc_time, ALL_TIME
        logical :: null_bn

		all_step = all_step_
        if(mod(all_step, 2) == 1) all_step = all_step + 1   ! Делаем общее число шагов чётным

        now = 2                        ! Какие параметры сейчас будут считаться (1 или 2). Они меняются по очереди
        time = 0.000002_8               ! Начальная инициализация шага по времени 
        Ncell = size(SS%gl_all_Cell(1, :))
        qqq1 = 0.0
        qqq2 = 0.0
        null_bn = .False.
        ALL_TIME = 0.0

        do step = 1, all_step

            if (mod(step, 10) == 0) then
                print*, "Step = ", step, " dt = ", time
            end if

            now2 = now
            now = mod(now, 2) + 1

            TT = time
            time = 1000000.0
            ALL_TIME = ALL_TIME + TT

            ! Предлагаю пробегаться сразу по ячейкам, а не по граням

            do cell = 1, Ncell

                center = SS%gl_Cell_Centr(:, cell, now)

                if(norm2(center) < 10.0) CYCLE

                par1 = SS%gd(:, cell, now)    ! Получили газодинамические параметры
                POTOK = 0.0

                ! Пробегаемся по граням
                do gr = 1, 4
                    gran = SS%gl_Cell_gran(gr, cell)
                    sosed = SS%gl_Cell_neighbour(gr, cell)
                    ! Задаём параметры соседа для распада разрыва
                    if(sosed == 0) then
                        CYCLE
                    else if(sosed == -1) then
                        call Phys_input_flow(SS, par2)
                    else if(sosed == -2) then
                        par2 = par1
                    else if(sosed == -3) then
                        par2 = par1
                    else if(sosed == -4) then
                        par2 = par1
                        par2(4) = -par2(4)
                    else
                        par2 = SS%gd(:, sosed, now)    ! Получили газодинамические параметры
                    end if

                    normal = SS%gl_Gran_normal(:, gran, now)
                    Sqv = SS%gl_Gran_length(gran, now)

                    qqq1(1) = par1(1)
                    qqq1(5) = par1(2)
                    qqq1(2) = par1(3)
                    qqq1(3) = par1(4)
                    qqq1(4) = 0.0
                    qqq1(6:8) = 0.0
                    qqq1(9) = par1(5)

                    qqq2(1) = par2(1)
                    qqq2(5) = par2(2)
                    qqq2(2) = par2(3)
                    qqq2(3) = par2(4)
                    qqq2(4) = 0.0
                    qqq2(6:8) = 0.0
                    qqq2(9) = par2(5)

                    call chlld_Q(0, normal(1), normal(2), 0.0_8, &
                    0.0_8, qqq1, qqq2, dsl, dsp, dsc, POTOK2, .False.)

                    loc_time = 0.99 * Sqv/( max(dabs(dsl), dabs(dsp)) )

                    POTOK(5) = POTOK(5) + POTOK2(9) * Sqv
                    POTOK(1) = POTOK(1) + POTOK2(1) * Sqv
                    POTOK(3) = POTOK(3) + POTOK2(2) * Sqv
                    POTOK(4) = POTOK(4) + POTOK2(3) * Sqv
                    POTOK(2) = POTOK(2) + POTOK2(5) * Sqv

                    time = min(time, loc_time)

                end do

                Vol = SS%gl_Cell_square(cell, now)

                ! Законы сохранения в ячейке
                ro = par1(1) - TT * POTOK(1) / Vol
                Q = par1(5) - TT * POTOK(5) / Vol
                u = (par1(1) * par1(3) - TT * POTOK(3) / Vol) / ro
                v = (par1(1) * par1(4) - TT * POTOK(4) / Vol) / ro
                p = ((  ( par1(2) / (SS%par_ggg - 1.0) + 0.5 * par1(1) * norm2(par1(3:4))**2 )  &
                        - TT * POTOK(2)/ Vol) - 0.5 * ro * (u**2 + v**2) ) * (SS%par_ggg - 1.0)

                SS%gd(1, cell, now2) = ro
                SS%gd(2, cell, now2) = p
                SS%gd(3, cell, now2) = u
                SS%gd(4, cell, now2) = v
                SS%gd(5, cell, now2) = Q

			end do
			
		

        end do

        print*, "ALL_TIME = ", ALL_TIME
		pause

    end subroutine Start_GD_algoritm


    subroutine Algoritm_Bound_condition(SS)
        !! Задаём граничные условия в сетке (на внутренней сфере для газовой динамики)
        TYPE (Setka), intent(in out) :: SS
        integer(4) :: n1, i, cell
        real(8) :: x, y, par(5)

        n1 = size(SS%gl_Cell_A(1,:))
        do i = 1, n1
            cell = SS%gl_Cell_A(1,i)
            x = SS%gl_Cell_Centr(1, cell, 1)
            y = SS%gl_Cell_Centr(2, cell, 1)
            call Inner_Conditions(SS, x, y, par)
            SS%gd(:, cell, 1) = par
            SS%gd(:, cell, 2) = SS%gd(:, cell, 1)
        end do

        n1 = size(SS%gl_Cell_B(1,:))
        do i = 1, n1
            cell = SS%gl_Cell_B(1,i)
            x = SS%gl_Cell_Centr(1, cell, 1)
            y = SS%gl_Cell_Centr(2, cell, 1)
            call Inner_Conditions(SS, x, y, par)
            SS%gd(:, cell, 1) = par

            SS%gd(:, cell, 2) = SS%gd(:, cell, 1)
        end do

	end subroutine Algoritm_Bound_condition

    subroutine Algoritm_Initial_condition(SS)
        !! Задаём начальные условия в сетке
        TYPE (Setka), intent(in out) :: SS
        integer(4) :: N, i
        real(8) :: x, y, par(5)

        N = size(SS%gl_Cell_Centr(1, :, 1))

        do i = 1, N
            x = SS%gl_Cell_Centr(1, i, 1)
            y = SS%gl_Cell_Centr(2, i, 1)

            call Phys_Innitial_Conditions(SS, x, y, par)

            SS%gd(:, i, 1) = par

            SS%gd(:, i, 2) = SS%gd(:, i, 1)
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

        call Geo_Culc_normal(SS, 1)
        call Geo_Culc_normal(SS, 2)
        call Culc_Cell_Centr(SS, 1)
        call Culc_Cell_Centr(SS, 2)

        call Belong_Init(SS)

        call Geo_Culc_length_area(SS, 1)
        call Geo_Culc_length_area(SS, 2)
		
	end subroutine Algoritm_ReMove_Surface


end module Algoritm