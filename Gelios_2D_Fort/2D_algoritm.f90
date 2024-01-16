module Algoritm
    USE STORAGE
    USE GEOMETRY
    USE Phys_parameter
	USE SURFACE
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

        call Algoritm_Bound_condition(SS, par)  ! «ададим граничные услови€ на внутренней сфере

    end subroutine Gas_dynamic_algoritm


    subroutine Algoritm_Bound_condition(SS, par)
        !! «адаЄм граничные услови€ в сетке (на внутренней сфере дл€ газовой динамики)
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

	end subroutine Algoritm_Bound_condition
	
	
    subroutine Algoritm_ReMove_Surface(SS, SURF)
	    ! ѕередвигает поверхности сетки согласно поверхност€м в SURF
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
		
		do ii = 1, 5  !! Ётот цикл нужен, т.к. HP ищетс€ по старой x - координате, нужно делать итерации
            do j = 1, N2

                the = par_pi/2.0 + (j) * SS%par_triple_point/(N2)
		        the2 = par_pi/2.0 + (j) * SS%par_triple_point_2/(N2)

                R_TS = SUR_GET_TS(SURF, the)
                yz = SS%gl_RAY_B(SS%par_n_HP, j)
                xx = SS%gl_yzel(1, yz, 1)
                yy = SUR_GET_HP_x(SURF, xx)
                ! «десь дл€ HP нужна еЄ y - координата (т.е. yy - это просто высота)
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

        do ii = 1, 5  !! Ётот цикл нужен, т.к. HP ищетс€ по старой x - координате, нужно делать итерации
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
            the = par_pi/2.0 + SS%par_triple_point + (N2 - j + 1) * (par_pi/2.0 - SS%par_triple_point)/(N2)  !TODO небезопасно, лучше организовать функцию, котора€ по номеру луча будет выдавать его угол
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