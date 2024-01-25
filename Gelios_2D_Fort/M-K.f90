module Monte_Karlo
    use STORAGE 
    use GEOMETRY
    USE ieee_arithmetic
    USE My_func
    implicit none 


    contains

    subroutine M_K_start(SS, SS_int)
        USE OMP_LIB
        ! Variables
        TYPE (Setka), intent(in out) :: SS
        TYPE (Inter_Setka), intent(in) :: SS_int
        real(8) :: start_time, end_time
        integer(4) :: step, iter, potok
        real(8) :: mu_(par_n_zone + 1), Wt_(par_n_zone + 1), Wp_(par_n_zone + 1), Wr_(par_n_zone + 1), X_(par_n_zone + 1)
		logical :: bb
        integer(4) :: i, cell
        real(8) :: sin_, x, phi, y, z, ksi, Vx, Vy, Vz, r_peregel, no, ksi1, ksi2, ksi3, ksi4, ksi5

		integer(4) :: num, to_i, to_j, j, pp, k
		real(8) :: ll, rr, Vphi, Vr, pui_w2, pui_w1
		real(8), allocatable :: vol_sr(:)                                    ! Для осреднения в узлах
		real(8), allocatable :: M_K_Moment_print(:, :, :)
		real(8) :: PAR(9), MAS_PUI(2)
		integer :: nk

        call M_K_Set()    ! Создали массивы
        call M_K_init()   ! Инициализируем веса и т.д.

        end_time = 0.0
		start_time = 0.0
        call omp_set_num_threads(par_n_potok)
		start_time = omp_get_wtime()
        step = 1
        cell = 1

        call Get_sensor_sdvig()

        do iter = 1, par_n_potok * par_n_parallel

            potok = (omp_get_thread_num() + 1) 

            !$omp critical
                print*, "start potok = ", potok, " iter = ", iter, "   step = ", step, "from = ", par_n_potok * par_n_parallel
                step = step + 1
			!$omp end critical

            ! Запускаем частицы первого типа (с полусферы)
            do num = 1, MK_N1
                call MK_Init_Parametrs(SS, potok, mu_, Wt_, Wp_, Wr_, X_, bb)

                do i = 1, par_n_zone + 1
                    sin_ = sqrt(1.0 - (X_(i)**2))
					x = (par_Rmax) * X_(i)
					call M_K_rand(SS%sensor(1, 1, potok), SS%sensor(2, 1, potok), SS%sensor(3, 1, potok), ksi)
					phi = 2.0 * par_pi * ksi
					y = (par_Rmax) * sin_ * cos(phi)
					z = (par_Rmax) * sin_ * sin(phi)
					
                    call Geo_Find_Cell(SS, x, sqrt(y**2 + z**2), cell)
					call dekard_skorost(x, y, z, Wr_(i), Wp_(i), Wt_(i), Vx, Vy, Vz)
				
					if(cell < 1) then
						STOP "Error cell < 1 MK 0lhy976yihko"
					end if
					
					
					if(i /= par_n_zone + 1 .or. bb == .True.) then
						! Добавляем частицу в стек
						SS%stek(potok) = SS%stek(potok) + 1
						SS%M_K_particle(1:7, SS%stek(potok), potok) = (/ x, y, z, Vx, Vy, Vz, mu_(i) * SS%MK_mu1 * SS%MK_Mu_mult /)
						SS%M_K_particle_2(1, SS%stek(potok), potok) = cell       ! В какой ячейке находится
						SS%M_K_particle_2(2, SS%stek(potok), potok) = SS%gl_all_Cell_zone(cell) ! Сорт
						call MK_Distination(SS, SS%M_K_particle(1:3, SS%stek(potok), potok), SS%M_K_particle(4:6, SS%stek(potok), potok),&
							to_i, to_j, r_peregel)
						SS%M_K_particle(8, SS%stek(potok), potok) = r_peregel
						SS%M_K_particle_2(3, SS%stek(potok), potok) = to_i  ! Зона назначения
						SS%M_K_particle_2(4, SS%stek(potok), potok) = to_j  ! Зона назначения
					end if

                    call M_K_Fly(SS, SS_int, potok)

                end do


            end do


        end do


    end subroutine M_K_start

    subroutine M_K_Fly(SS, SS_int, potok)
        TYPE (Setka), intent(in out) :: SS
        TYPE (Inter_Setka), intent(in) :: SS_int
        integer(4), intent(in) :: potok

	end subroutine M_K_Fly


    subroutine M_K_Set()
	
	end subroutine M_K_Set


    subroutine M_K_init()
		
	end subroutine M_K_init


    subroutine Get_sensor_sdvig()
		
	end subroutine Get_sensor_sdvig

	subroutine M_K_Change_Velosity4(SS, potok, Ur, Uthe, Uphi, Vr, Vthe, Vphi, Wr_, Wthe_, Wphi_, mu_, cp, r, I_, x_ex, y_ex, z_ex, bb)
		! Как первая часть, но розыгрышь идёт по-частям
		! Запускаесть один основной и I_ дополнительных атомов
        TYPE (Setka), intent(in out) :: SS
		integer(4), intent(in) :: potok, I_
		real(8), intent(in) ::  Ur, Uthe, Uphi, Vr, Vthe, Vphi, cp, r, x_ex, y_ex, z_ex
		real(8), intent(out) ::   Wr_(I_ + 1), Wthe_(I_ + 1), Wphi_(I_ + 1), mu_(I_ + 1)
		logical, intent(out) :: bb
		
		real(8) :: X, uu
		real(8) :: gamma_(I_), Wa_(I_), Mho_(I_)
		real(8) :: ksi, gam1, gam2, Wr1, Wr2, Wr0, ksi1, ksi2, W1, W2, Wa, pp1, pp2, c, p, u
		real(8) :: p4, om1, om2, om3, lo, y1, y2, y3, v1, v2, v3, u1, u2, u3, uuu, yy, h, ksi3, ksi4, ksi5, ksi6, D, ko, gg
		integer(4) :: ii, met, k
		
		X = sqrt( (Vr - Ur)**2 + (Vthe - Uthe)**2 + (Vphi - Uphi)**2 )
		uu = exp(-(X**2)) / par_sqrtpi + (X + 1.0 / (2.0 * X)) * erf(X)
		Wr0 = -1.0
		met = 0

		do ii = 1, I_
			gamma_(ii) = 1.0 / ( (r / MK_R_zone(ii))**2 - 1.0)
		end do

		

		! Разыграем Wr
		
		do ii = 1, I_
			
			if (ii == 1) then
				gam1 = 0.0
				gam2 = gamma_(ii)
			else
				gam1 = gamma_(ii - 1)
				gam2 = gamma_(ii)
			end if
			
			
			
			pp1 = MK_for_Wr_1(0.0_8, gam2, Ur) - MK_for_Wr_1(0.0_8, gam1, Ur)
			pp2 = MK_for_Wr_2(0.0_8, gam2, Ur, Uthe) - MK_for_Wr_2(0.0_8, gam1, Ur, Uthe)
			
			
			call M_K_rand(SS%sensor(1, 2, potok), SS%sensor(2, 2, potok), SS%sensor(3, 2, potok), ksi)

			if (ksi < pp1 / (pp1 + pp2)) then
				call M_K_rand(SS%sensor(1, 2, potok), SS%sensor(2, 2, potok), SS%sensor(3, 2, potok), ksi)
				Wr1 = -3.0;
				Wr2 = 0.0;

				do while (MK_H_Wr_1(gam1, gam2, Wr1, Ur, pp1, ksi) >= 0.0)
					Wr1 = Wr1 - 1.0
				end do
					
				k = 0
				do while (dabs(Wr2 - Wr1) > 0.0001)     ! Деление пополам, иначе разваливается
					Wr0 = (Wr1 + Wr2) / 2.0
					if (MK_H_Wr_1(gam1, gam2, Wr0, Ur, pp1, ksi) < 0) then
						Wr1 = Wr0
					else
						Wr2 = Wr0
					end if
					k = k + 1
				end do
				Wr_(ii) = Wr1
			else
				call M_K_rand(SS%sensor(1, 2, potok), SS%sensor(2, 2, potok), SS%sensor(3, 2, potok), ksi)
				Wr1 = -3.0
				Wr2 = 0.0

				do while (MK_H_Wr_2(gam1, gam2, Wr1, Ur, Uthe, pp2, ksi) >= 0.0)
					Wr1 = Wr1 - 1.0
				end do
				
				k = 0
				do while (dabs(Wr2 - Wr1) > 0.0001)     ! Деление пополам, иначе разваливается
					Wr0 = (Wr1 + Wr2) / 2.0
					if (MK_H_Wr_2(gam1, gam2, Wr0, Ur, Uthe, pp2, ksi) < 0) then
						Wr1 = Wr0
					else
						Wr2 = Wr0
					end if
					k = k + 1
				end do
				Wr_(ii) = Wr1
			end if
		end do

		
		! Разыгрываем  Wa
		do ii = 1, I_
			
			if (ii == 1) then
				gam1 = 0.0
				gam2 = gamma_(ii)
			else
				gam1 = gamma_(ii - 1)
				gam2 = gamma_(ii)
			end if

			do while(.True.)
				call M_K_rand(SS%sensor(1, 2, potok), SS%sensor(2, 2, potok), SS%sensor(3, 2, potok), ksi1)
				call M_K_rand(SS%sensor(1, 2, potok), SS%sensor(2, 2, potok), SS%sensor(3, 2, potok), ksi2)
				W1 = sqrt(gam1 * (Wr_(ii)**2))
				W2 = sqrt(gam2 * (Wr_(ii)**2))
				Wa = sqrt(-log(exp(-(W1**2)) - ksi1 * (exp(-(W1**2)) - exp(-(W2**2)))))
				if ((1.0 + (Uthe * Wa)**2) / (1.0 + (Uthe * W2)**2) >= ksi2) EXIT
			end do
			Wa_(ii) = Wa
		end do

		
		! Считаем веса и Mho
		do ii = 1, I_
			c = Uthe * Wa_(ii)
			p = MK_norm_mho(c)
			
			Mho_(ii) = MK_play_mho(potok, c)
			

			Wthe_(ii) = Wa_(ii) * cos(Mho_(ii))
			Wphi_(ii) = Wa_(ii) * sin(Mho_(ii))
			
			if (ii == 1) then
				gam1 = 0.0
				gam2 = gamma_(ii)
			else
				gam1 = gamma_(ii - 1)
				gam2 = gamma_(ii)
			end if

			u = sqrt( (Vr - Wr_(ii))**2 + (Vthe - Wthe_(ii))**2 +(Vphi - Wphi_(ii))**2 )

			if (X > 7.0) then
				mu_(ii) = (u * MK_sigma2(u, cp) / (uu * MK_sigma2(uu, cp))) * (MK_f2(0.0_8, gam1, Ur, Uthe) - MK_f2(0.0_8, gam2, Ur, Uthe)) * exp(-(Uthe**2)) * (p) / (1.0 + (c**2))
			else
				mu_(ii) = (u * MK_sigma2(u, cp) / (MK_int_1(X * cp, cp) / cp)) * (MK_f2(0.0_8, gam1, Ur, Uthe) - MK_f2(0.0_8, gam2, Ur, Uthe)) * exp(-(Uthe**2)) * (p) / (1.0 + (c**2))
			end if
		end do


		! Розыгрыш основного атома
		p4 = 0.5 * par_sqrtpi * X / (1.0 + 0.5 * par_sqrtpi * X)
		

		if (I_ > 0) then
			gg = gamma_(I_)
		else
			gg = 0.0
		end if

		do while(.True.)
			call M_K_rand(SS%sensor(1, 2, potok), SS%sensor(2, 2, potok), SS%sensor(3, 2, potok), ksi1)
			call M_K_rand(SS%sensor(1, 2, potok), SS%sensor(2, 2, potok), SS%sensor(3, 2, potok), ksi2)
			call M_K_rand(SS%sensor(1, 2, potok), SS%sensor(2, 2, potok), SS%sensor(3, 2, potok), ksi3)
			call M_K_rand(SS%sensor(1, 2, potok), SS%sensor(2, 2, potok), SS%sensor(3, 2, potok), ksi4)
			call M_K_rand(SS%sensor(1, 2, potok), SS%sensor(2, 2, potok), SS%sensor(3, 2, potok), ksi5)
			call M_K_rand(SS%sensor(1, 2, potok), SS%sensor(2, 2, potok), SS%sensor(3, 2, potok), ksi6)
			
			if (p4 < ksi1) then
				om1 = 1.0 - 2.0 * ksi4
				om2 = sqrt(1.0 - (om1**2)) * cos(2.0 * par_pi * ksi5)
				om3 = sqrt(1.0 - (om1**2)) * sin(2.0 * par_pi * ksi5)

				lo = sqrt(-log(ksi2 * ksi3))
				y1 = lo * om1
				y2 = lo * om2
				y3 = lo * om3
			else
				y1 = sqrt(-log(ksi2)) * cos(par_pi * ksi3)
				y2 = sqrt(-log(ksi4)) * cos(2.0 * par_pi * ksi5)
				y3 = sqrt(-log(ksi4)) * sin(2.0 * par_pi * ksi5)
			end if
			
			v1 = y1 + Ur
			v2 = y2 + Uthe
			v3 = y3 + Uphi
			u1 = Vr - v1
			u2 = Vthe - v2
			u3 = Vphi - v3
			uuu = sqrt(kvv(u1, u2, u3))
			yy = sqrt(kvv(y1, y2, y3))
			h = ((uuu * MK_sigma2(uuu, cp)) / (MK_sigma2(X, cp) * (X + yy)))
			
			if (h >= ksi6) EXIT
		end do


		Wr_(I_ + 1) = v1
		Wthe_(I_ + 1) = v2
		Wphi_(I_ + 1) = v3


		if (Wr_(I_ + 1) >= 0.0 .or. (Wthe_(I_ + 1)**2) + (Wphi_(I_ + 1)**2) > gg * (Wr_(I_ + 1)**2)) then
			mu_(I_ + 1) = 1.0
			bb = .True.
			return
		else
			mu_(I_ + 1) = 0.0  ! Чтобы не запускать этот атом
			bb = .False.
			return
		end if

		bb = .True.
		return
	end subroutine M_K_Change_Velosity4

	subroutine MK_Velosity_initial2(SS, potok, Vx, Vy, Vz)
        TYPE (Setka), intent(in out) :: SS
		integer(4), intent(in) :: potok
		real(8), intent(out) :: Vx, Vy, Vz
		
		real(8) :: ksi1, ksi2, a, ksi3, ksi4, ksi5, ksi6, z, p1
		
		call M_K_rand(SS%sensor(1, 1, potok), SS%sensor(2, 1, potok), SS%sensor(3, 1, potok), ksi1)
		call M_K_rand(SS%sensor(1, 1, potok), SS%sensor(2, 1, potok), SS%sensor(3, 1, potok), ksi2)
		
		a = sqrt(-log(1.0 - ksi2))
		Vy = a * cos(2.0 * par_pi * ksi1)
		Vz = a * sin(2.0 * par_pi * ksi1)
		
		z = 0
		p1 = 0.5 * dabs(SS%par_Velosity_inf) * par_sqrtpi / (0.5 + 0.5 * dabs(SS%par_Velosity_inf) * par_sqrtpi);

		do
			call M_K_rand(SS%sensor(1, 1, potok), SS%sensor(2, 1, potok), SS%sensor(3, 1, potok), ksi3)
			call M_K_rand(SS%sensor(1, 1, potok), SS%sensor(2, 1, potok), SS%sensor(3, 1, potok), ksi4)
			call M_K_rand(SS%sensor(1, 1, potok), SS%sensor(2, 1, potok), SS%sensor(3, 1, potok), ksi5)
			call M_K_rand(SS%sensor(1, 1, potok), SS%sensor(2, 1, potok), SS%sensor(3, 1, potok), ksi6)

			if (p1 > ksi3) then
				z = cos(par_pi * ksi5) * sqrt(-log(ksi4))
			else
				z = sqrt(-log(1.0 - ksi4))
			end if
			
			if((dabs(z + SS%par_Velosity_inf) / (dabs(SS%par_Velosity_inf) + dabs(z)) > ksi6 .and. z >= -SS%par_Velosity_inf)) EXIT
		end do

		Vx = z + SS%par_Velosity_inf
		
		if (Vx <= 0.0) then
			print*, "Error iuygvbnmklo9890pljiouytrtyjhg"
		end if
		
		return
	end subroutine MK_Velosity_initial2
	
	subroutine MK_Velosity_initial(SS, potok, Vx, Vy, Vz)
        TYPE (Setka), intent(in out) :: SS
		integer(4), intent(in) :: potok
		real(8), intent(out) :: Vx, Vy, Vz
		
		real(8) :: ksi1, ksi2, a, ksi3, ksi4, ksi5, ksi6, z, p1
		
		call M_K_rand(SS%sensor(1, 1, potok), SS%sensor(2, 1, potok), SS%sensor(3, 1, potok), ksi1)
		call M_K_rand(SS%sensor(1, 1, potok), SS%sensor(2, 1, potok), SS%sensor(3, 1, potok), ksi2)
		
		a = sqrt(-log(ksi2))
		Vy = a * cos(2.0 * par_pi * ksi1)
		Vz = a * sin(2.0 * par_pi * ksi1)
		
		z = 0
		p1 = dabs(SS%par_Velosity_inf) * par_sqrtpi / (1.0 + dabs(SS%par_Velosity_inf) * par_sqrtpi)

		do
			call M_K_rand(SS%sensor(1, 1, potok), SS%sensor(2, 1, potok), SS%sensor(3, 1, potok), ksi3)
			call M_K_rand(SS%sensor(1, 1, potok), SS%sensor(2, 1, potok), SS%sensor(3, 1, potok), ksi4)
			call M_K_rand(SS%sensor(1, 1, potok), SS%sensor(2, 1, potok), SS%sensor(3, 1, potok), ksi5)
			call M_K_rand(SS%sensor(1, 1, potok), SS%sensor(2, 1, potok), SS%sensor(3, 1, potok), ksi6)

			if (p1 > ksi3) then
				z = cos(par_pi * ksi5) * sqrt(-log(ksi4))
			else
				if (ksi4 <= 0.5) then
					z = -sqrt(-log(2.0 * ksi4))
				else
					z = sqrt(-log(2.0 * (1.0 - ksi4)))
				end if
			end if
			
			if(dabs(z + SS%par_Velosity_inf) / (dabs(SS%par_Velosity_inf) + dabs(z)) >= ksi6 .and. z <= -SS%par_Velosity_inf) EXIT
		end do

		Vx = z + SS%par_Velosity_inf
		
		return
	end subroutine MK_Velosity_initial

    subroutine MK_ruletka(SS,n_potok, to_i, to_j, from_i, from_j, area, r, r_per, mu3, bb2)
        TYPE (Setka), intent(in out) :: SS
		real(8), intent(in out) :: mu3
		real(8), intent(in) :: r, r_per
		integer(4), intent(in) :: n_potok, to_i, to_j, from_i, from_j, area
		logical, intent(out) :: bb2
		
		real(8) :: mu, ksi
		
		mu = min( MK_Mu(to_i, to_j, area) * 0.3 * MK_SINKR(to_j) * (r_per/par_Rmax)**2, &
					MK_Mu(from_i, from_j, area) * 0.3 * MK_SINKR(from_j) * (r/par_Rmax)**2) 
		
		if (mu3 >= mu) then
			bb2 = .True.
			return
		else
			call M_K_rand(SS%sensor(1, 2, n_potok), SS%sensor(2, 2, n_potok), SS%sensor(3, 2, n_potok), ksi)
			if (mu3 >= ksi * mu) then
				mu3 = mu
				bb2 = .True.
				return
			else
				mu3 = 0.0
				bb2 = .False.
				return
			end if
		end if
		
		bb2 = .False.
		return
	
	end subroutine MK_ruletka

    subroutine MK_Init_Parametrs(SS, potok, mu_, Wt_, Wp_, Wr_, X_, bb)
		!! Розыгрышь атомов на начальной сфере (скорости и положения)
        TYPE (Setka), intent(in out) :: SS
		integer(4), intent(in) :: potok
		real(8), intent(out) :: mu_(par_n_zone + 1), Wt_(par_n_zone + 1), Wp_(par_n_zone + 1), Wr_(par_n_zone + 1), X_(par_n_zone + 1)
		logical, intent(out) :: bb
		
		real(8) :: Y, Wr1, Wr2, Wr0, W1, W2, Wa, Yt, c, p
		real(8) :: ksi, ksi1, ksi2
		real(8) :: X1
		real(8) :: X2
		real(8) :: X0                ! для метода хорд
		real(8) :: split, gam1, gam2
		real(8) :: Wa_(par_n_zone + 1)
		real(8) :: Mho_(par_n_zone + 1)
		integer(4) :: i
		real(8) :: p1, ksi3, ksi4, z, h, gg, p4_, ksi5, ksi6
		
		
		X0 = 1.0
		Y = dabs(SS%par_Velosity_inf)
		p1 = erf(Y) / (SS%MK_A1_ * (Y**2))
		
		

		! Разыгрываем  X
		do i = 1, par_n_zone
			
			call M_K_rand(SS%sensor(1, 1, potok), SS%sensor(2, 1, potok), SS%sensor(3, 1, potok), ksi)
			X1 = 0.0
			X2 = 1.0
			
			if (i == 1) then
				gam1 = 0.0
				gam2 = SS%MK_gam_zone(i)
			else
				gam1 = SS%MK_gam_zone(i - 1)
				gam2 = SS%MK_gam_zone(i)
			end if
			

			do while (dabs(X2 - X1) > 0.0000001)     ! Деление пополам, иначе разваливается
				X0 = (X1 + X2) / 2.0
				if (MK_Hx(gam1, gam2, X0, Y, ksi) < 0.0) then
					X1 = X0
				else
					X2 = X0
				end if
			end do
			X_(i) = X0
		end do

		! Разыгрываем  Wr
		do i = 1, par_n_zone
			
			call M_K_rand(SS%sensor(1, 1, potok), SS%sensor(2, 1, potok), SS%sensor(3, 1, potok), ksi)
			
			Wr1 = -5.0
			Wr2 = 0.0
			
			if (i == 1) then
				gam1 = 0.0
				gam2 = SS%MK_gam_zone(i)
			else
				gam1 = SS%MK_gam_zone(i - 1)
				gam2 = SS%MK_gam_zone(i)
			end if
			
			do while (MK_Hwr(gam1, gam2, Wr1, X_(i), Y, ksi) >= 0.0)
				Wr1 = Wr1 - 1.0
			end do
				
			do while (dabs(Wr2 - Wr1) > 0.000001)     ! Деление пополам, иначе разваливается
				Wr0 = (Wr1 + Wr2) / 2.0
				if (MK_Hwr(gam1, gam2, Wr0, X_(i), Y, ksi) < 0) then
					Wr1 = Wr0
				else
					Wr2 = Wr0
				end if
			end do

			Wr_(i) = Wr0
		end do

		
		! Разыгрываем  Wa
		
		do i = 1, par_n_zone
			Yt = Y * sqrt(1.0 - (X_(i)**2))
			
			if (i == 1) then
				gam1 = 0.0
				gam2 = SS%MK_gam_zone(i)
			else
				gam1 = SS%MK_gam_zone(i - 1)
				gam2 = SS%MK_gam_zone(i)
			end if

			do while(.True.)
				call M_K_rand(SS%sensor(1, 1, potok), SS%sensor(2, 1, potok), SS%sensor(3, 1, potok), ksi1)
				call M_K_rand(SS%sensor(1, 1, potok), SS%sensor(2, 1, potok), SS%sensor(3, 1, potok), ksi2)
				W1 = sqrt(gam1 * (Wr_(i)**2))
				W2 = sqrt(gam2 * (Wr_(i)**2))
				Wa = sqrt(-log(exp(-(W1**2)) - ksi1 * (exp(-(W1**2)) - exp(-(W2**2)))))
				if ((1.0 + (Yt * Wa)**2) / (1.0 + (Yt * W2)**2) >= ksi2) EXIT
			end do


			Wa_(i) = Wa
		end do

		! Считаем веса
		do i = 1, par_n_zone
			
			c = Y * sqrt(1.0 - (X_(i)**2)) * Wa_(i)
			p = MK_norm_mho(c)
			
			
			Mho_(i) = MK_play_mho(SS, potok, c)
			

			Wt_(i) = Wa_(i) * cos(Mho_(i))
			Wp_(i) = Wa_(i) * sin(Mho_(i))
			
			if (i == 1) then
				gam1 = 0.0
				gam2 = SS%MK_gam_zone(i)
			else
				gam1 = SS%MK_gam_zone(i - 1)
				gam2 = SS%MK_gam_zone(i)
			end if
			
			mu_(i) = (MK_F(1.0_8, gam2, Y) - MK_F(1.0_8, gam1, Y)) * (p) / (SS%MK_A0_ * (1.0 + (c**2)))
		end do

		! Разыгрываем основной атом

		! Розыгрыш X целиком из препринта

		
		call M_K_rand(SS%sensor(1, 1, potok), SS%sensor(2, 1, potok), SS%sensor(3, 1, potok), ksi1)
		
		
		if (p1 < ksi1) then
			
			do while(.True.)
				call M_K_rand(SS%sensor(1, 1, potok), SS%sensor(2, 1, potok), SS%sensor(3, 1, potok), ksi2)
				call M_K_rand(SS%sensor(1, 1, potok), SS%sensor(2, 1, potok), SS%sensor(3, 1, potok), ksi3)
				X2 = sqrt(ksi2)
				h = (1.0 + erf(X2 * Y)) / (1.0 + erf(Y))
				if (h > ksi3) EXIT
			end do
		else
			do
				call M_K_rand(SS%sensor(1, 1, potok), SS%sensor(2, 1, potok), SS%sensor(3, 1, potok), ksi2)
				call M_K_rand(SS%sensor(1, 1, potok), SS%sensor(2, 1, potok), SS%sensor(3, 1, potok), ksi3)
				X2 = (1.0 / Y) * sqrt(-log(ksi2)) * cos(par_pi * ksi3 / 2.0)
				if(X2 < 1.0) EXIT
			end do
		end if
		

		X_(par_n_zone + 1) = X2

		gg = 0.0
		
		if (par_n_zone > 0) gg = SS%MK_gam_zone(par_n_zone)
		
		p4_ = par_sqrtpi * (Y * X2) / (1.0 + par_sqrtpi * (Y * X2))                 ! Для розыгрыша основного атома на границе по препринту Маламы
		
		do while(.True.)
			call M_K_rand(SS%sensor(1, 1, potok), SS%sensor(2, 1, potok), SS%sensor(3, 1, potok), ksi1)
			call M_K_rand(SS%sensor(1, 1, potok), SS%sensor(2, 1, potok), SS%sensor(3, 1, potok), ksi2)
			call M_K_rand(SS%sensor(1, 1, potok), SS%sensor(2, 1, potok), SS%sensor(3, 1, potok), ksi3)
			if (p4_ > ksi1) then
				z = sqrt(-log(ksi2)) * cos(par_pi * ksi3)
			else
				if (ksi2 <= 0.5) then
					z = -sqrt(-log(2.0 * ksi2))
				else
					z = sqrt(-log(2.0 * (1.0 - ksi2)))
				end if
			end if

			Wr_(par_n_zone + 1) = z - (Y * X2)
			h = dabs(-(Y * X2) + z) / ((Y * X2) + dabs(z))
			call M_K_rand(SS%sensor(1, 1, potok), SS%sensor(2, 1, potok), SS%sensor(3, 1, potok), ksi4)
			if (h > ksi4 .and. z < (Y * X2)) EXIT
		end do

			
			call M_K_rand(SS%sensor(1, 1, potok), SS%sensor(2, 1, potok), SS%sensor(3, 1, potok), ksi5)
			call M_K_rand(SS%sensor(1, 1, potok), SS%sensor(2, 1, potok), SS%sensor(3, 1, potok), ksi6)

			Wt_(par_n_zone + 1) = Y * sqrt(1.0 - (X_(par_n_zone + 1)**2)) + sqrt(-log(ksi5)) * cos(2.0 * par_pi * ksi6)
			Wp_(par_n_zone + 1) = sqrt(-log(ksi5)) * sin(2.0 * par_pi * ksi6)


		if (par_n_zone == 0) then
			mu_(par_n_zone + 1) = 1.0
			bb = .True.
			return
		else
			if (Wr_(par_n_zone + 1) >= 0.0 .or. Wt_(par_n_zone + 1)**2 + Wp_(par_n_zone + 1)**2 > &
				SS%MK_gam_zone(par_n_zone) * Wr_(par_n_zone + 1)**2) then
				mu_(par_n_zone + 1) = 1.0
			else
				mu_(par_n_zone + 1) = 0.0  ! Чтобы не запускать этот атом
				bb = .False.
				return
			end if
		end if

		bb = .True.
		return
	end subroutine MK_Init_Parametrs

    subroutine M_K_rand(s1, s2, s3, b)
		integer(4), intent(in out) :: s1, s2, s3
		real(8), intent(out) :: b
		integer(4):: ic15 = 32768, ic10 = 1024
		integer(4):: mz = 710, my = 17784, mx = 11973
		real(8):: xi = 9.0949470177292824E-13, c = 1.073741824E9
		integer(4) :: i13, i12, i11, ii
		i13 = mz * s1 + my * s2 + mx * s3
		i12 = my * s1 + mx * s2
		i11 = mx * s1
		ii = i11 / ic15
		i12 = i12 + ii
		s1 = i11 - ic15 * ii
		ii = i12 / ic15
		i13 = i13 + ii
		s2 = i12 - ic15 * ii
		s3 = mod(i13,ic10)
		b = xi * (c * s3 + ic15 * s2 + s1)
	end subroutine M_K_rand

    subroutine MK_Distination(SS, r, V, to_i, to_j, rr)
		! Variables
        TYPE (Setka), intent(in) :: SS
		real(8), intent(in) :: r(3), V(3)
		real(8), intent(out) :: rr
		integer(4), intent(out) :: to_i, to_j
		
		real(8) :: time_do_peregel, rk(3)
		
		time_do_peregel = -DOT_PRODUCT(r, V) / kvv(V(1), V(2), V(3))
		rk = r + time_do_peregel * V
		rr = norm2(rk)
		
		to_i = MK_geo_zones(SS, rr, 1.0_8)
		to_j = MK_alpha_zones(SS, polar_angle( rk(1), sqrt(rk(2)**2 + rk(3)**2) ) )
	
	end subroutine MK_Distination

    integer(4) pure function MK_alpha_zones(SS, al)
        TYPE (Setka), intent(in) :: SS
		real(8), intent(in) :: al
		integer(4) :: i
		
		do i = 1, par_m_zone
			if(al < SS%MK_al_zone(i)) then
				MK_alpha_zones = i
				return
			end if
		end do
		
		MK_alpha_zones = par_m_zone + 1
		return
		
	end function MK_alpha_zones
	
	integer(4) pure function MK_geo_zones(SS, r, k)
		! В какой геометрической зоне сейчас находится атом (зоны считаются с нуля)
        TYPE (Setka), intent(in) :: SS
		real(8), intent (in) :: r, k
		integer(4) :: i

		do i = 1, par_n_zone
			if (r < k * SS%MK_R_zone(i)) then
				MK_geo_zones = i
				return
			end if
		end do
		
		MK_geo_zones = par_n_zone + 1
		
	end function MK_geo_zones

    real(8) pure function MK_F(X, gam, Y)
		real(8), intent (in) :: X, gam, Y
		real(8) :: bb, gg, Y2, Y4, X2, A, B, C, D
		
		bb = sqrt(1.0 + gam)
		gg = 1.0 + gam
		Y2 = Y**2
		Y4 = Y**4
		X2 = X**2
		A = exp(-Y2) / (2.0 * par_sqrtpi * Y * bb**7)
		B = 2 * X * Y * bb**3 * (2.0 + gam * (3.0 + (X2 - 1.0) * Y2 + gam))
		C = par_sqrtpi * gg * (gg * (gg * (4.0 + gam) + Y2 * (2.0 + 3.0 * gam)) + &
			exp(X2 * Y2 / gg) * (2.0 * X2 * (X2 - 1.0) * Y4 * gam - (gg**2) * (4.0 + gam) + &
				Y2 * gg * (-2.0 - 3.0 * gam + X2 * (2.0 + gam))))
		D = par_sqrtpi * gg * exp(X2 * Y2 / gg) * &
			(-2.0 * X2 * (X2 - 1.0) * Y4 * gam + gg**2 * (4.0 + gam) - Y2 * gg * (-2.0 - 3.0 * gam + X2 * (2.0 + gam))) * &
			erf(X * Y / bb)
			
		MK_F =  A * (B + C - D)
	end function MK_F
	
	real(8) pure function MK_FI(Z, X, gam, Y)
		real(8), intent (in) :: Z, X, gam, Y
		real(8) :: bb, gg, X2, Y2, Y4, Z2, A, B, C, D, ee
		
		bb = sqrt(1.0 + gam)
		gg = 1.0 + gam
		X2 = X**2
		Y2 = Y**2
		Y4 = Y**4
		Z2 = Z**2
		A = exp(-Y2 - 2.0 * X * Y * Z - Z2 * gg) / (par_sqrtpi * bb**7)
		ee = exp((X * Y + Z + Z * gam)**2 / gg)
		B = ee * par_sqrtpi * X * Y * (2.0 * X2 * (X2 - 1.0) * Y4 * gam - 2.0 * gg**2 + (X2 - 1.0) * Y2 * gg * (2.0 + 5.0 * gam))
		C = 2.0 * bb * (X2 * (X2 - 1.0) * Y4 * gam - X * (X2 - 1.0) * Y**3 * Z * gam * gg - gg**2 +  &
			(X2 - 1.0) * Y2 * gg * (1 + gam * (2.0 + Z2 * gg)))
		D = ee * par_sqrtpi * X * Y * (2.0 * X2 * (X2 - 1.0) * Y4 * gam - 2.0 * gg**2 + (X2 - 1.0) * Y2 * gg * (2.0 + 5.0 * gam)) * &
			erf((X * Y + Z + Z * gam) / bb)
			
		MK_FI =  A * (B + C + D)
	end function MK_FI
	
	real(8) pure function MK_Hx(gam1, gam2, X, Y, ksi)
		real(8), intent (in) :: gam1, gam2, X, Y, ksi
		MK_Hx = MK_F(X, gam2, Y) - MK_F(X, gam1, Y) - ksi * (MK_F(1.0_8, gam2, Y) - MK_F(1.0_8, gam1, Y))
	end function MK_Hx
	
	real(8) pure function MK_Hwr(gam1, gam2, Z, X, Y, ksi)
		real(8), intent (in) :: gam1, gam2, Z, X, Y, ksi
		MK_Hwr = MK_FI(Z, X, gam2, Y) - MK_FI(Z, X, gam1, Y) - ksi * (MK_FI(0.0_8, X, gam2, Y) - MK_FI(0.0_8, X, gam1, Y))
	end function MK_Hwr
	
	real(8) pure function MK_norm_mho(c)
		real(8), intent (in) :: c
		real(8) :: s, d, cc, nn
		integer(4) :: i
		
		s = 1.0
		d = 0.0
		cc = 1.0
		nn = 1.0
		
		do i = 1, 99
			cc = cc * c
			nn = nn * i
			d = (cc / nn)**2
			s = s + d
			if (d < 0.00001) then
				EXIT
			end if
		end do
		
		MK_norm_mho = s
	end function MK_norm_mho
	
	real(8) function MK_play_mho(SS, potok, c)
        TYPE (Setka), intent(in out) :: SS
		real(8), intent(in) :: c
		integer(4), intent(in) :: potok
		real(8) :: x, ksi
		
		do while(.True.)
			call M_K_rand(SS%sensor(1, 2, potok), SS%sensor(2, 2, potok), SS%sensor(3, 2, potok), ksi)
			x = 2.0 * ksi * par_pi
			call M_K_rand(SS%sensor(1, 2, potok), SS%sensor(2, 2, potok), SS%sensor(3, 2, potok), ksi)
			if(exp(2.0 * c * cos(x)) > ksi * exp(2.0 * dabs(c)) ) EXIT
		end do

		MK_play_mho = x
	end function MK_play_mho

	real(8) pure function MK_Velosity_1(u, cp)
		real(8), intent (in) :: u, cp
		
		if (u < 0.00001) then
			MK_Velosity_1 = 2.0 * cp / par_sqrtpi + 2.0 * u * u / (3.0 * cp * par_sqrtpi) - &
				u * u * u * u / (15.0 * cp * cp * cp * par_sqrtpi)
		else
			MK_Velosity_1 =  exp(-u * u / cp**2) * cp / par_sqrtpi + (u + (cp**2) / (2.0 * u)) * erf(u / cp)
		end if
		
	end function MK_Velosity_1
	
	real(8) pure function MK_Velosity_2(u, cp)
		real(8), intent (in) :: u, cp
		
		if (u < 0.00001) then
			MK_Velosity_2 = (8.0 / 3.0) * cp**4 * par_pi * u + (8.0 / 15.0) * cp**2 * par_pi * u * u * u - &
				(4.0 / 105.0) * par_pi * u**5
		else
			MK_Velosity_2 =  cp**3 * par_pi * (exp(-u * u / cp**2) * cp * u * 2.0 * (cp**2 + 2.0 * u**2) + &
			par_sqrtpi * (4.0 * u**4 + 4.0 * cp * cp * u**2 - cp**4) * erf(u / cp)) / (4.0 * u * u)
		end if
		
	end function MK_Velosity_2
	
	real(8) pure function MK_Velosity_3(u, cp)
		real(8), intent (in) :: u, cp
		
		if (u < 0.00001) then
			MK_Velosity_3 = 8.0 * cp / (3.0 * par_sqrtpi) + 8.0 * u**2 / (9.0 * cp * par_sqrtpi) - &
				44.0 * u**4 / (135.0 * cp * cp * cp * par_sqrtpi)
		else
			MK_Velosity_3 =  exp(-u**2 / cp**2) * cp * (5.0 * cp**2 + 2.0 * u**2) / (par_sqrtpi * (3.0 * cp**2 + 2.0 * u**2)) + &
			(4.0 * u**4 + 12.0 * cp**2 * u**2 + 3.0 * cp**4) * erf(u / cp) / (2.0 * u * (3.0 * cp**2 + 2.0 * u**2))
		end if
		
	end function MK_Velosity_3
	
	real(8) pure function MK_int_1_f1(x)

		real(8), intent (in) :: x
	
		if (x <= 1.0) then
			MK_int_1_f1 =  6.283185155644284 + 0.000024846677279866114 * x + &
				2.0934329078277405 * x * x + 0.008055998193903208 * x * x * x - &
				0.2355169235647438 * x * x * x * x + 0.03820480582423355 * x * x * x * x * x + &
				0.006992274370591744 * x * x * x * x * x * x
		else if (x <= 3.0) then
			MK_int_1_f1 =  6.437524091973454 - 0.6331520099380095 * x + &
				3.1348881317268997 * x * x - 0.8454201478027856 * x * x * x + &
				0.1004702004260311 * x * x * x * x + 0.0009895488638964746 * x * x * x * x * x - &
				0.000920750276197054 * x * x * x * x * x * x
		else if (x <= 5) then
			MK_int_1_f1 =  4.4920780630505135 + 2.5133093267020654 * x + &
				1.1327223176567935 * x * x - 0.24648691152318875 * x * x * x + &
				0.031326738629523766 * x * x * x * x - 0.0021366031960331384 * x * x * x * x * x + &
				0.00005954097505746697 * x * x * x * x * x * x
		else if (x <= 7) then
			MK_int_1_f1 =  1.9138683588136232 + 5.350374732905213 * x - &
				0.16380205801427633 * x * x + 0.06765898334856263 * x * x * x - &
				0.011071118267864083 * x * x * x * x + 0.0008673476933852199 * x * x * x * x * x - &
				0.00002691859374483661 * x * x * x * x * x * x
		else if (x <= 50.0) then
			MK_int_1_f1 =  1.3138472469154294 + 5.336877156136497 * x + &
				0.020286308991329983 * x * x - 0.9780973533544137 * (x / 10.0)**3 + &
				0.26354051936651874 * (x / 10.0)**4 - 0.03711733070841712 * (x / 10.0)**5 + &
				0.002120935433043921 * (x / 10.0)**6
		else
			MK_int_1_f1 =  1.3138472469154294 + 5.336877156136497 * x + &
				0.020286308991329983 * x * x - 0.9780973533544137 * (x / 10.0)**3 + &
				0.26354051936651874 * (x / 10.0)**4 - 0.03711733070841712 * (x / 10.0)**5 + &
				0.002120935433043921 * (x / 10.0)**6
		end if
				 
		return
	end function MK_int_1_f1
	
	real(8) pure function MK_int_1_f2(x)
		real(8), intent (in) :: x
		
		if (x <= 1.0) then
			MK_int_1_f2 =  1.328216167939543 - 0.000004545681954848391 * x + &
				2.537368073155103 * x * x - 0.0020584991728545624 * x * x * x - &
				0.03742568018912792 * x * x * x * x - 0.010312136385277346 * x * x * x * x * x + &
				0.002767736179209713 * x * x * x * x * x * x
		else if (x <= 3.0) then
			MK_int_1_f2 =  1.2959616295629246 + 0.1533684067037866 * x + &
				2.2354849981206106 * x * x + 0.3113395567715921 * x * x * x - &
				0.21656309882941488 * x * x * x * x + 0.041957500887605075 * x * x * x * x * x - &
				0.0029978773724628604 * x * x * x * x * x * x
		else if (x <= 5.0) then
			MK_int_1_f2 =  1.903643456971281 - 1.4801836911099535 * x + 3.973958664572268 * x * x - &
				0.6482729779428982 * x * x * x + 0.07665007314658864 * x * x * x * x - &
				0.005369758193338703 * x * x * x * x * x + 0.00016605531871992049 * x * x * x * x * x * x
		else if (x <= 7.0) then
			MK_int_1_f2 =  -4.484415105552316 + 5.3747429756690766 * x + &
				0.8892806582308143 * x * x + 0.09767316152573671 * x * x * x - &
				0.025704749778475783 * x * x * x * x + 0.0021937998296249206 * x * x * x * x * x - &
				0.00006928845984076111 * x * x * x * x * x * x
		end if
				 
		return
	end function MK_int_1_f2
	
	real(8) pure function MK_int_1_f3(x)
		real(8), intent (in) :: x
 
		if (x <= 1.0) then
			MK_int_1_f3 = 1.2938345594193854 - 0.000031719847351174835 * x + &
				1.3183710041280094 * x * x - 0.014150512069488197 * x * x * x + &
				0.4226114681928129 * x * x * x * x - 0.06985750969880078 * x * x * x * x * x - &
				0.015347864048406958 * x * x * x * x * x * x
		else if (x <= 3.0) then
			MK_int_1_f3 = 0.9667460440956788 + 1.336271810704016 * x - &
				0.8687355257991665 * x * x + 1.7676868273627229 * x * x * x - &
				0.2731222764016417 * x * x * x * x + 0.004801770033831665 * x * x * x * x * x + &
				0.001780776080720323 * x * x * x * x * x * x
		else if (x <= 5.0) then
			MK_int_1_f3 = 4.760566734123174 - 5.048204299463048 * x + 3.332342585228025 * x * x + &
				0.47584339615235993 * x * x * x - 0.12072786272726124 * x * x * x * x + &
				0.011870955604980658 * x * x * x * x * x - 0.0004580199652402304 * x * x * x * x * x * x
		else if (x <= 7.0) then
			MK_int_1_f3 = 9.370493362469261 - 10.848615619383615 * x + 6.423326878282571 * x * x - &
				0.4148977656870439 * x * x * x + 0.025300923044176957 * x * x * x * x - &
				0.0010108688120876522 * x * x * x * x * x + 0.00001864423130429156 * x * x * x * x * x * x
		end if
			 
		return
	end function MK_int_1_f3
	
	real(8) pure function MK_int_1(x, cp)
		real(8), intent (in) :: x, cp
		real(8) :: b
		
		b = 1.0 - par_a_2 * log(cp)
		MK_int_1 = (cp / (par_sqrtpi**3)) * (b * b * MK_int_1_f1(x / cp) - &
			2.0 * par_a_2 * b * MK_int_1_f2(x / cp) + par_a_2**2 * MK_int_1_f3(x / cp))
	end function MK_int_1
	
	real(8) pure function MK_int_2_f1(x)
		real(8), intent (in) :: x
	
		if (x <= 1.0) then
				MK_int_2_f1 =   8.377571213788123 * x + 0.00047608508679086725 * x * x + &
					1.6710478320575737 * x * x * x + 0.016857530811432916 * x * x * x * x - &
					0.15132474125724812 * x * x * x * x * x + 0.030723378194358945 * x * x * x * x * x * x
				return
		else if (x <= 3.0) then
				MK_int_2_f1 =   -0.11788367995598747 + 8.937936705157014 * x - &
						1.119886471634001 * x * x + 2.8831031948885917 * x * x * x - &
						0.735250146386514 * x * x * x * x + 0.10356311378423572 * x * x * x * x * x - &
						0.006231417172309398 * x * x * x * x * x * x
				return
		else if (x <= 5.0) then
				MK_int_2_f1 =   2.9044497739429858 + 2.757712415967557 * x + 4.239161941189675 * x * x + &
						0.36198838294786784 * x * x * x - 0.05737777787138304 * x * x * x * x + &
						0.004956250079677106 * x * x * x * x * x - 0.0001809238236975877 * x * x * x * x * x * x
				return
		else if (x <= 7.0) then
				MK_int_2_f1 =  41.6323689028892 - 38.118317864344135 * x + 22.211189528076645 * x * x -&
					3.8547348524931246 * x * x * x + 0.5000517174807501 * x * x * x * x -&
					0.03446294709493891 * x * x * x * x * x + 0.0009860204070962582 * x * x * x * x * x * x
				return
		end if
	
		MK_int_2_f1 =  0.0
	end  function MK_int_2_f1
	
	real(8) pure function MK_int_2_f2(x)

		real(8), intent (in) :: x
	
		if (x <= 1.0) then
			MK_int_2_f2 =  3.8653461103376667 * x + 0.0001975300512691014 * x * x + &
				2.4468141895384012 * x * x * x + 0.005987984681429616 * x * x * x * x - &
				0.06453987836713967 * x * x * x * x * x + 0.0066920981111229004 * x * x * x * x * x * x
			return
		else if (x <= 3.0) then
			MK_int_2_f2 =  -0.10983889480661446 + 4.321087890898017 * x - &
				0.7707850845797619 * x * x + 3.1237901158486583 * x * x * x - &
				0.31485222316123385 * x * x * x * x + 0.010270249760261791 * x * x * x * x * x + &
				0.0008259803934338584 * x * x * x * x * x * x
			return
		else if (x <= 5.0) then
			MK_int_2_f2 =  1.8468011847509729 + 1.1986396743254275 * x + &
				1.1421489029509448 * x * x + 2.606316149569781 * x * x * x - &
				0.2788783468089509 * x * x * x * x + 0.019815317035281846 * x * x * x * x * x - &
				0.0006379970557448899 * x * x * x * x * x * x
			return
		else if (x <= 7.0) then
			MK_int_2_f2 =  9.480707804348185 - 8.022988228952784 * x + 5.823555900242488 * x * x + &
				1.3277220473440972 * x * x * x - 0.08074921612981413 * x * x * x * x + &
				0.0033058587723954185 * x * x * x * x * x - 0.00006041810279926061 * x * x * x * x * x * x
			return
		end if
			 
		MK_int_2_f2 = 0.0
	end function MK_int_2_f2
		
	real(8) pure function MK_int_2_f3(x)
        real(8), intent (in) :: x
		if (x <= 1.0) then
			MK_int_2_f3 =  2.6106039258326 * x - 0.0008357997793049243 * x * x + &
				2.0764571907368174 * x * x * x - 0.03182275644841273 * x * x * x * x + &
				0.26310521962808975 * x * x * x * x * x - 0.06034325471643871 * x * x * x * x * x * x
			return
		else if (x <= 3.0) then
			MK_int_2_f3 =   0.20784760901369737 + 1.5920325291316857 * x + &
				2.0985329535259014 * x * x - 0.26286255221171206 * x * x * x + &
				1.4610329096120886 * x * x * x * x - 0.25626862029131897 * x * x * x * x * x + &
				0.01684969647300594 * x * x * x * x * x * x
			return
		else if (x <= 5.0) then
			MK_int_2_f3 =   -6.284115352064703 + 15.665162343948523 * x &
				- 10.766105772158252 * x**2 + 6.0821074614870465 * x**3 - &
				0.3181501403196319 * x**4 + 0.01232319194701587 * x**5 &
				- 0.00017890661550876597 * x**6
			return
		else if (x <= 7.0) then
			MK_int_2_f3 =   -4.355962170454177 + 13.835332665069274 * x &
				- 10.12766071646978 * x**2 + 5.999392227482686 * x**3 - &
				0.32171318647989955 * x**4 + 0.014181987261856027 * x**5 &
				- 0.00030579035551447497 * x**6
			return
		end if
				
		MK_int_2_f3 = 0.0
		
	end function MK_int_2_f3
	
	real(8) pure function MK_int_2(x, cp)
		real(8), intent (in) :: x, cp
		real(8) :: b
		
		b = 1.0 - par_a_2 * log(cp)
		MK_int_2 = -(cp * cp / (par_sqrtpi * par_sqrtpi * par_sqrtpi)) * (b * b * MK_int_2_f1(x / cp) - &
			2.0 * par_a_2 * b * MK_int_2_f2(x / cp) + (par_a_2**2) * MK_int_2_f3(x / cp))
	end function MK_int_2
	
	real(8) pure function MK_int_3(x, cp)
		real(8), intent (in) :: x, cp
		real(8) :: b
		
		b = 1.0 - par_a_2 * log(cp)
		MK_int_3 = (cp**3 / (par_sqrtpi**3)) * (b * b * MK_int_3_f1(x / cp) - 2.0 * par_a_2 * b * MK_int_3_f2(x / cp) + (par_a_2**2) * MK_int_3_f3(x / cp))
	end function MK_int_3

	real(8) pure function MK_int_3_f1(x)
		real(8), intent (in) :: x
		
		if (x <= 1.0) then
			MK_int_3_f1 = 12.566370586001975 - 0.00001944816384202852 * x + 12.567558607381049 * x**2 - &
				0.010507444068349692 * x**3 + 1.2911398125420694 * x**4 - 0.05048482363937502 * x**5 - &
				0.029947937607835744 * x**6
		else if (x <= 3.0) then
			MK_int_3_f1 = 12.451555448799724 + 0.40252674353016715 * x + 12.081033298182223 * x**2 + 0.12939193776331415 * x**3 + &
				1.478876561367302 * x**4 - 0.2237491583356496 * x**5 + 0.014474521138620033 * x**6
		else if (x <= 5.0) then
			MK_int_3_f1 = 8.281962852844913 + 9.783032429527339 * x + 3.165848344887614 * x**2 + 4.711462666119614 * x**3 + &
				0.13739453130827392 * x**4 - 0.012096015315889195 * x**5 + 0.0004514225555943018 * x**6
		else if (x <= 7.0) then
			MK_int_3_f1 = 17.29966669035025 + 2.1457227895928668 * x + 5.572082327920818 * x**2 + 4.4083748449004645 * x**3 + &
				0.13640200155890422 * x**4 - 0.00854302147917508 * x**5 + 0.00022205921430504255 * x**6
		end if
		
	end function MK_int_3_f1

	real(8) pure function MK_int_3_f2(x)
		real(8), intent (in) :: x
		
		if (x <= 1.0) then
			MK_int_3_f2 = 5.798024979296493 - 0.00001772671478406096 * x + 9.98769025405073 * x**2 - 0.0073593516014156535 * x**3 + &
				2.27820901023418 * x**4 - 0.03135086958655956 * x**5 - 0.030403716978821237 * x**6
		else if (x <= 3.0) then
			MK_int_3_f2 = 5.864728705834779 - 0.34799875480550213 * x + 10.760249318358127 * x**2 - 0.9493738978205943 * x**3 + &
				2.948810798239708 * x**4 - 0.29690823082625284 * x**5 + 0.015284639719532296 * x**6
		else if (x <= 5.0) then
			MK_int_3_f2 = -3.21810405152587 + 17.08054504204813 * x - 3.43427364926184 * x**2 + 5.342946243936558 * x**3 + &
				1.3459970589832742 * x**4 - 0.07448103623442687 * x**5 + 0.0021616888853670958 * x**6
		else if (x <= 7.0) then
			MK_int_3_f2 = -59.42860988375287 + 79.24709914283142 * x - 32.304295129821256 * x**2 + 12.558114389999929 * x**3 + &
				0.32138208180756495 * x**4 + 0.003976327853989966 * x**5 - 0.0003703205838879336 * x**6
		end if
		
	end function MK_int_3_f2

	real(8) pure function MK_int_3_f3(x)
		real(8), intent(in) :: x
		
		if (x <= 1.0) then
			MK_int_3_f3 = 3.915885322797866 + 0.000044284982651632276 * x + 7.778946280540595 * x**2 + 0.01990833982643192 * x**3 + &
				2.710531013102172 * x**4 + 0.09480498076917274 * x**5 + 0.026454483470905545 * x**6
		else if (x <= 3.0) then
			MK_int_3_f3 = 4.371710139648201 - 1.787430730965795 * x + 10.530386754306065 * x**2 - 1.9945949141813966 * x**3 + &
				3.310710912254515 * x**4 + 0.1356460090430464 * x**5 - 0.019853464614841342 * x**6
		else if (x <= 5.0) then
			MK_int_3_f3 = 4.6940208896644435 - 6.423081982183987 * x + 18.137713177954524 * x**2 - 7.298291873534797 * x**3 + &
				5.202095367169497 * x**4 - 0.20620743501245353 * x**5 + 0.005094092519200813 * x**6
		else if (x <= 7.0) then
			MK_int_3_f3 = 194.2826492092263 - 195.77327124311523 * x + 95.83253215419927 * x**2 - 23.99281397751645 * x**3 + &
				7.170549820159753 * x**4 - 0.32568903981020303 * x**5 + 0.007955090180048264 * x**6
		end if
		
	end function MK_int_3_f3

	real(8) pure function MK_for_Wr_1(Z, gam, ur)
		!! Для розыгрыша Wr в перезарядке по-частям (первая часть)
		real(8), intent (in) :: Z, gam, ur
		
		MK_for_Wr_1 = -exp(-gam * ur**2 / (gam + 1.0)) * (1.0 + erf((-ur + gam * Z + Z) / sqrt(gam + 1.0))) / (2.0 * sqrt(gam + 1.0))
	end function MK_for_Wr_1

	real(8) pure function MK_for_Wr_2(Z, gam, ur, ut)
		!! Для розыгрыша Wr в перезарядке по-частям 
		real(8), intent (in) :: Z, gam, ur, ut
		real(8) :: gam1, ur2
		
		gam1 = gam + 1.0  ! gam + 1
		ur2 = ur**2
		
		MK_for_Wr_2 =  1.0 / (4.0 * gam1 * gam1 * gam1 * par_sqrtpi) * ut**2 * exp((-1.0 - gam / gam1) * ur2 - gam1 * Z * Z) * &
			(2.0 * exp(ur * (gam * ur / gam1 + 2.0 * Z)) * gam * gam1 * (ur + Z + gam * Z) + &
				exp(ur2 + gam1 * Z * Z) * sqrt(gam1) * par_sqrtpi * (2.0 + gam * (5.0 + 3.0 * gam + 2.0 * ur2)) * &
				(-1.0 - erf((-ur + Z + gam * Z) / sqrt(gam1))))
	end function MK_for_Wr_2

	real(8) pure function MK_H_Wr_1(gam1, gam2,	V, ur, p, ksi)
		real(8), intent (in) :: gam1, gam2, V, ur, p, ksi
		MK_H_Wr_1 = MK_for_Wr_1(V, gam2, ur) - MK_for_Wr_1(V, gam1, ur) - ksi * p
	end function MK_H_Wr_1

	real(8) pure function MK_H_Wr_2(gam1, gam2, V, ur, ut, p, ksi)
		! Для розыгрыша Wr в перезарядке по-частям 
		real(8), intent (in) :: gam1, gam2, V, ur, ut, p, ksi
		MK_H_Wr_2 =  MK_for_Wr_2(V, gam2, ur, ut) - MK_for_Wr_2(V, gam1, ur, ut) - ksi * p
	end function MK_H_Wr_2

	real(8) pure function MK_f2(V, gam, ur, ut)
		real(8), intent (in) :: V, gam, ur, ut
		real(8) :: b, b4, b2, A
		
		b = sqrt(1.0 + gam)
		b4 = b**4
		b2 = b**2

		A = (exp(2.0 * V * ur - (ur**2) - (V**2) * b2) / (2.0 * par_sqrtpi * b4))
		
		if (A <= 0.0000000000001) then
			MK_f2 = 0.0
			return
		end if
			
		MK_f2 = A * (-gam * (ut**2) * (ur + gam * V + V) + &
			(par_sqrtpi / (2.0 * b)) * exp( (V + gam * V - ur)**2 / b2) * (2.0 * b4 + (ut**2) * (b2 * (3.0 * gam + 2.0) + 2.0 * gam * (ur**2))) * &
			(1.0 + erf((-ur + gam * V + V) / b)))
	end function MK_f2


end module Monte_Karlo