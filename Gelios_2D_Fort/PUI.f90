module PUI
    USE STORAGE
    USE GEOMETRY
    USE OMP_LIB
	USE Interpol
    implicit none 

    contains

	!! Функция инициализации находится в модуле STORAGE


    subroutine PUI_Add(SS, cell, wr, nu_ex, mu, time)
		! wr - скорость атома в СК, связанной со средней скоростью плазмы (модуль этой скорости)
        TYPE (Setka), intent(in out) :: SS
		integer, intent(in) :: cell
		real(8), intent(in) :: wr, nu_ex, mu, time
		integer i, j

		j = SS%f_pui_num2(cell)
		if(j > 0) then
			i = min(INT(wr/SS%pui_wR * SS%pui_nW) + 1, SS%pui_nW)

			call omp_set_lock(SS%pui_lock(j))
			SS%pui_Sm(i, j) = SS%pui_Sm(i, j) + mu * time
			SS%pui_Sp(i, j) = SS%pui_Sp(i, j) + mu * time * nu_ex
			call omp_unset_lock(SS%pui_lock(j))
		end if
	end subroutine PUI_Add

    subroutine PUI_calc_Sm(SS)
		! Расчёт Sm - он не считается "на лету" в Монте-карло и нужна постобработка
        TYPE (Setka), intent(in out) :: SS
		real(8) :: pui_Sm2(SS%pui_nW)           ! (pui_nW, :, potok)
		integer :: ij, i, j, k, num_all
		real(8) :: dthe, Vh, ff, the, d, w

		dthe = par_pi/40.0
		num_all = 0
		print*, "Start PUI_calc_Sm"
	 	!$omp parallel

	 	!$omp do private(pui_Sm2, ij, j, k, ff, Vh, the, d, w)
		do i = 1, size(SS%pui_Sm(1, :))

			!$omp critical
			num_all = num_all + 1
			if(mod(num_all, 5000) == 0) then
				print*, num_all, "from", size(SS%pui_Sm(1, :))
			end if
			!$omp end critical

			pui_Sm2 = 0.0
			do ij = 1, SS%pui_nW
				w = ((ij - 0.5) * SS%pui_wR / SS%pui_nW);
				do j = 1, SS%pui_nW
					Vh = ((j - 0.5) * SS%pui_wR / SS%pui_nW);
					ff = SS%pui_Sm(j, i)
					if (ff <= 0.0) CYCLE
					do k = 1, 40
						the = dthe * k
						d = sqrt(Vh**2 + w**2 - 2.0 * w * Vh * cos(the))
						if (d > 0.000000001) then
							pui_Sm2(ij) = pui_Sm2(ij) + ff * d * MK_sigma(SS, d) * sin(the) * dthe * 2.0 * par_pi
						end if
					end do
				end do
				pui_Sm2(ij) = pui_Sm2(ij) / (4.0 * par_pi)
			end do

			SS%pui_Sm(:, i) = pui_Sm2/SS%par_Kn !/2.0  !TODO НУЖНО БУДЕТ УБРАТЬ ДЕЛЕНИЕ НА ДВА В СЛЕДУЮЩИЙ РАЗ !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		end do
		!$omp end do
		!$omp end parallel

		print*, "End PUI_calc_Sm"

	end subroutine PUI_calc_Sm

    subroutine PUI_print_S(SS, x, y)
        TYPE (Setka), intent(in) :: SS
        real(8), intent(in) :: x, y
        integer(4) :: num, n2, i
        real(8) :: w
        LOGICAL :: inzone
        character(len=5) :: name

        num = 1
        call Geo_Find_Cell(SS, x, y, num, inzone)

        if(SS%gl_all_Cell_zone(num) > 2) return
        write(unit=name,fmt='(i5.5)') num

        n2 = SS%f_pui_num2(num)
        open(3, file = "S_PUI_" // name // ".txt")
		if(n2 > 0) then
			do i = 1, SS%pui_nW
				w = (i-0.5) * SS%pui_wR/SS%pui_nW
				write(3,*) w, SS%pui_Sm(i, n2), SS%pui_Sp(i, n2)
			end do
		end if
		close(3)


    end subroutine PUI_print_S

	subroutine PUI_print_pui(SS, x, y)
        TYPE (Setka), intent(in) :: SS
        real(8), intent(in) :: x, y
        integer(4) :: num, n2, i
        real(8) :: w
        LOGICAL :: inzone
        character(len=5) :: name

        num = 1
        call Geo_Find_Cell(SS, x, y, num, inzone)

        if(SS%gl_all_Cell_zone(num) > 2) return
        write(unit=name,fmt='(i5.5)') num

        n2 = SS%f_pui_num2(num)
        open(3, file = "f_PUI_" // name // ".txt")
		if(n2 > 0) then
			do i = 1, SS%pui_nW
				w = (i-0.5) * SS%pui_wR/SS%pui_nW
				write(3,*) w, SS%f_pui(i, n2)
			end do
		end if
		close(3)


    end subroutine PUI_print_pui

	subroutine Culc_f_pui(SS, SS_int)
		TYPE (Setka), intent(in out) :: SS
		TYPE (Inter_Setka), intent(in) :: SS_int
		integer :: i, k, num, yzel, tetraedron, iw, numw, num_do, i1, num2, cell_n, iter, cell
		real(8) :: r(2), PAR(SS%n_par), dt, Sm, w0, q1, rho, rho0, w, qInt, Sp, rho_do, rr
		real(8) :: aa(3), bb(3), cc(3)
		real(8) :: s, cospsi, C, A, B
		real(8) :: PAR_k(4), normal(3)
		real(8) :: f0_pui(SS%pui_nW)
		real(8) :: mas_w(SS%pui_nW)
		real(8) :: mas_w0(SS%pui_nW)
		real(8) :: mas_Sm(SS%pui_nW)
		real(8) :: mas_Sm2(SS%pui_nW)
		logical :: find_n, inzone

		SS%pui_Sm = SS%pui_Sm * SS%par_n_H_LISM
		SS%pui_Sp = SS%pui_Sp * SS%par_n_H_LISM


		print*, "START Culc_f_pui"

		dt = 0.01/50.0

		! Находим функцию распределения для ячеек перед ударной волной
		print*, "Do TS"
		!$omp parallel
	 	!$omp do private(cell, iter, k, num, yzel, tetraedron, iw, numw, num_do, i1, num2, r, PAR, Sm, w0, q1, rho, rho0, w, qInt, Sp, rho_do, s, cospsi, C, A, B, PAR_k, normal, f0_pui, mas_w, mas_w0, mas_Sm, mas_Sm2, find_n)
		do i = 1, size(SS%f_pui_num)
			if(mod(i, 300) == 0) print*, "i = ", i, " from", size(SS%f_pui_num)
			k = SS%f_pui_num(i)      ! Номер ячейки
			if(SS%gl_all_Cell_zone(k) == 2) CYCLE  ! Пока обрабатываем точки внутри TS

			do iw = 1, SS%pui_nW
				mas_w0(iw) = ((iw - 0.5) * SS%pui_wR / SS%pui_nW)
			end do

			rho0 = SS%gd(1, k, 1)
			f0_pui = 0.0
			mas_Sm = 0.0

			r = SS%gl_Cell_Centr(:, k, 1)  ! Координаты этого узла
			mas_w = mas_w0

			num = 3
			cell = k
			qInt = 0.0            ! Интеграл от источника массы при ионизации

			! Бежим до Солнца
			iter = 1
			do while (.TRUE.)
				iter = iter + 1
				if(iter > 100000) then
					continue
				end if
				!call Int2_Get_par_fast(r(1), r(2), r(3), num, PAR, PAR_k = PAR_k)
				call Int_Get_Parameter(SS_int, r(1), r(2), num, PAR_gd = PAR, PAR_atom_source = PAR_k)
				call Geo_Find_Cell(SS, r(1), r(2), cell, inzone)
				if(inzone == .False.) then
					print*, "ERROR wvretbyunim67574321c4tv5y"
					pause
					STOP
				end if
				q1 = PAR_k(4)
				rho = PAR(1)
				qInt = qInt + dt * q1/rho
				r = r - PAR(3:4) * dt
				if(dt * norm2(PAR(3:4)) > 0.01) then
					print*, "ERROR wr4tv5byvwr4c3qe"
					pause
					STOP
				end if
				if(r(2) <= 0.0) r(2) = 0.00001
				tetraedron = SS%f_pui_num2(cell) ! Номер тетраэдра в массиве источников

				do iw = 1, SS%pui_nW
					numw = min(INT(mas_w(iw)/SS%pui_wR * SS%pui_nW) + 1, SS%pui_nW)
					if(mas_w(iw) < SS%pui_wR .and. mas_w(iw) > 0) then
						mas_Sm(iw) = mas_Sm(iw) + SS%pui_Sm(numw, tetraedron) * dt
						f0_pui(iw) = f0_pui(iw) + SS%pui_Sp(numw, tetraedron) * dt * exp(-mas_Sm(iw))  ! Это S+, просто сразу накапливаем в функцию распределения
					end if
					mas_w(iw) = mas_w0(iw) / ( (rho0/rho)**(1.0/3.0) * exp(-1.0/3.0 * qInt) )
				end do

				if(norm2(r) <= SS%par_R0) EXIT
			end do

			SS%f_pui(:, i) = f0_pui(:)
		end do
		!$omp end do
		!$omp end parallel



		! Находим функцию распределения для ячеек за ударной волной
		print*, "Posle TS"
		!$omp parallel
	 	!$omp do private(cell, iter, aa, bb, cc, k, num, yzel, tetraedron, iw, numw, num_do, i1, num2, r, PAR, Sm, w0, q1, rho, rho0, w, qInt, Sp, rho_do, s, cospsi, C, A, B, PAR_k, normal, f0_pui, mas_w, mas_w0, mas_Sm, mas_Sm2, find_n)
		do i = 1, size(SS%f_pui_num)

			do iw = 1, SS%pui_nW
				mas_w0(iw) = ((iw - 0.5) * SS%pui_wR / SS%pui_nW)
			end do

			if(mod(i, 300) == 0) print*, "i = ", i, " from", size(SS%f_pui_num)
			k = SS%f_pui_num(i)      ! Номер узла сетки интерполяции, в которой считаем PUI
			if(SS%gl_all_Cell_zone(k) == 1) CYCLE  ! Пропускаем ячейки до TS
			rho0 = SS%gd(1, k, 1)
			rho = rho0
			f0_pui = 0.0
			mas_Sm = 0.0

			r = SS%gl_Cell_Centr(:, k, 1)  ! Координаты этого узла
			mas_w = mas_w0

			num = 1
			cell = k
			qInt = 0.0            ! Интеграл от источника массы при ионизации

			!print*, r
			!pause

			! Бежим до TS
			do while (.TRUE.)
				num_do = cell
				call Int_Get_Parameter(SS_int, r(1), r(2), num, PAR_gd = PAR, PAR_atom_source = PAR_k)
				call Geo_Find_Cell(SS, r(1), r(2), cell, inzone)
				if(inzone == .False.) then
					print*, "ERROR wvretbyunim67574321c4tv5y"
					pause
					STOP
				end if

				q1 = PAR_k(4)
				rho_do = rho
				rho = PAR(1)
				
				if(SS%gl_all_Cell_zone(cell) == 1) EXIT   ! Если попали в ячейку из области 1 - область до TS

				if(SS%gl_all_Cell_zone(cell) == 3) then   ! Если попали в ячейку из области 3 - область за HP
					!print*, "Popal v zonu 3"
					r = r * 0.999
					CYCLE
				end if

				qInt = qInt + dt * q1/rho
				r = r - PAR(3:4) * dt
				if(r(2) <= 0.0) r(2) = 0.00001
				
				tetraedron = SS%f_pui_num2(cell) ! Номер тетраэдра в массиве источников

				do iw = 1, SS%pui_nW
					numw = min(INT(mas_w(iw)/SS%pui_wR * SS%pui_nW) + 1, SS%pui_nW)
					if(mas_w(iw) < SS%pui_wR .and. mas_w(iw) > 0) then
						mas_Sm(iw) = mas_Sm(iw) + SS%pui_Sm(numw, tetraedron) * dt
						f0_pui(iw) = f0_pui(iw) + SS%pui_Sp(numw, tetraedron) * dt * exp(-mas_Sm(iw))  ! Это S+, просто сразу накапливаем в функцию распределения
					end if
					mas_w(iw) = mas_w0(iw) / ( (rho0/rho)**(1.0/3.0) * exp(-1.0/3.0 * qInt) )
				end do

			end do

			s = rho_do/rho
			C = s

			rho0 = rho
			qInt = 0.0
			mas_w0 = mas_w/sqrt(C)
			mas_w = mas_w0

			! Бежим до Солнца
			do while (.TRUE.)
				call Int_Get_Parameter(SS_int, r(1), r(2), num, PAR_gd = PAR, PAR_atom_source = PAR_k)
				call Geo_Find_Cell(SS, r(1), r(2), cell, inzone)
				q1 = PAR_k(4)
				rho = PAR(1)
				qInt = qInt + dt * q1/rho
				r = r - PAR(3:4) * dt
				if(r(2) <= 0.0) r(2) = 0.00001
				!print*, r
				!pause
				tetraedron = SS%f_pui_num2(cell) ! Номер тетраэдра в массиве источников

				do iw = 1, SS%pui_nW
					numw = min(INT(mas_w(iw)/SS%pui_wR * SS%pui_nW) + 1, SS%pui_nW)
					if(mas_w(iw) < SS%pui_wR .and. mas_w(iw) > 0) then
						mas_Sm(iw) = mas_Sm(iw) + SS%pui_Sm(numw, tetraedron) * dt
						f0_pui(iw) = f0_pui(iw) + SS%pui_Sp(numw, tetraedron) * dt * exp(-mas_Sm(iw)) * s/C**(1.5)  ! Это S+, просто сразу накапливаем в функцию распределения
					end if
					mas_w(iw) = mas_w0(iw) / ( (rho0/rho)**(1.0/3.0) * exp(-1.0/3.0 * qInt) )
				end do

				if(norm2(r) <= SS%par_R0) EXIT
			end do

			SS%f_pui(:, i) = f0_pui(:)

		end do
		!$omp end do
		!$omp end parallel

		print*, "END Culc_f_pui"
	end subroutine Culc_f_pui

	subroutine PUI_Culc_h0(SS)
		TYPE (Setka), intent(in out) :: SS
		integer :: i, j, k
		real(8) :: the, UH, w, u, S

		!$omp parallel
	 	!$omp do private(S, UH, u, the, w, i, j)
		do k = 1, SS%pui_h0_n
			S = 0.0
			UH = (k - 0.5) * SS%pui_wR/SS%pui_h0_n
			do i = 1, 1000
				w = (i - 0.5) * SS%pui_wR / 1000
				do j = 1, 180
					the = j * par_pi / 180.0
					u = sqrt((w * sin(the))**2 + (w * cos(the) - UH)**2)
					S = max(S, u * MK_sigma(SS, u)/((w + SS%pui_h0_wc) * MK_sigma(SS, w + SS%pui_h0_wc)) )
				end do
			end do
			SS%h0_pui(k) = S
		end do
		!$omp end do
		!$omp end parallel

		open(1, file = "h0_PUI.txt")
		do k = 1, SS%pui_h0_n
			UH = (k - 0.5) * SS%pui_wR/SS%pui_h0_n
			write(1, *) UH, SS%h0_pui(k)
		end do
		close(1)

		open(1, file = "n0_PUI_2d_UH=40.txt")
		UH = 30.0
		k = min(INT(UH/SS%pui_wR * SS%pui_h0_n) + 1, SS%pui_h0_n)
		do i = 1, 1000
			w = (i - 0.5) * SS%pui_wR / 1000
			do j = 1, 180
				the = j * par_pi / 180.0
				u = sqrt((w * sin(the))**2 + (w * cos(the) - UH)**2)
				S = u * MK_sigma(SS, u)/((w + SS%pui_h0_wc) * MK_sigma(SS, w + SS%pui_h0_wc))/SS%h0_pui(k)
				write(1, *) w, the, S
			end do
		end do
		close(1)

	end subroutine PUI_Culc_h0

	real(8) pure function PUI_get_f(SS, n, w)
		! Получить значение функции распределения PUI с номером n, со скоростью w (используется линейная интерполяция)
		TYPE (Setka), intent(in) :: SS
		real(8), intent (in) :: w
		integer, intent (in) :: n
		integer :: k1, k2
		real(8) :: x1, x2,  f1, f2
		k1 = max(min(INT(w/SS%pui_wR * SS%pui_nW + 0.5), SS%pui_nW), 1)
		k2 = min(k1 + 1, SS%pui_nW)
		x1 = (k1 - 0.5) * SS%pui_wR/SS%pui_nW
		x2 = (k2 - 0.5) * SS%pui_wR/SS%pui_nW
		f1 = SS%f_pui(k1, n)
		f2 = SS%f_pui(k2, n)
		!PUI_get_f = (f1 - f2)/(x1 - x2) * x + (f1 * x2 - f2 * x1)/(x2 - x1)

		! if(w > 5.0 .and.(w < x1 .or. w > x2)) then
		! 	print*, "ERROR PUI_get_f 390", w, x1, x2
		! end if

		if(w < x1 .or. dabs(x1 - x2) <= 0.0000001) then
			PUI_get_f = f1
		else if (w > x2) then
			PUI_get_f = f2
		else
			PUI_get_f = f1 * (w - x2)/(x1 - x2) + f2 * (w - x1)/(x2 - x1)
		end if
    end function PUI_get_f
	
	subroutine PUI_F_integr_Set(SS)
		! Создаём массивы для интегрирования функции распределения для последующего розыгрыша
		! и сразу вычисляет их
		! Также считаем концентрацию PUI
		TYPE (Setka), intent(in out) :: SS
		integer(4) :: N

		print*, "PUI_F_integr_Set  START"
		N = size(SS%f_pui(1, :))
		allocate(SS%F_integr_pui(SS%pui_F_n, N))
		allocate(SS%nu_integr_pui(SS%pui_F_n, N))
		allocate(SS%Mz_integr_pui(SS%pui_F_n, N))
		allocate(SS%E_integr_pui(SS%pui_F_n, N))
		print*, "PUI_F_integr_Set  END"
	end subroutine PUI_F_integr_Set

	subroutine PUI_F_integr_Culc(SS)
		! Создаём массивы для интегрирования функции распределения для последующего розыгрыша
		! и сразу вычисляет их
		! Также считаем концентрацию PUI
		TYPE (Setka), intent(in out) :: SS
		integer :: N, n2, i, step, j, k
		real(8) :: S, SSS, w, S1, S2, ff, UH, u, the

		print*, "PUI_F_integr_Culc  START"
		N = size(SS%f_pui(1, :))

		step = 0
		! Посчитаем этот интеграл (см. документацию PUI \rho_w)
		!$omp parallel
	 	!$omp do private(SSS, S, w, S1, S2, i, ff, UH, j, u, the, k)
		do n2 = 1, N
			!$omp critical
				step = step + 1
				if(mod(step, 500) == 0) then
					print*, step, "  from = ", N
				end if
			!$omp end critical
			SSS = 0.0

			do i = 1, 10000
				w = (i-0.5) * SS%pui_wR/10000
				ff = PUI_get_f(SS, n2, w)
				SSS = SSS + ff * w**2 * (w + SS%pui_h0_wc) * MK_sigma(SS, w + SS%pui_h0_wc) * (SS%pui_wR/10000)
			end do

			S1 = 0.0
			w = 0.0
			!! Далее алгоритм посложнее, там шаг по w должен быть как можно меньше
			do i = 1, SS%pui_F_n
				S2 = (i - 0.5) * 1.0/SS%pui_F_n
				!print*, "S2 = ", S2
				do while (S1 < S2)
					w = w + (SS%pui_wR/10000)
					ff = PUI_get_f(SS, n2, w)
					S1 = S1 + ff * w**2 * (w + SS%pui_h0_wc) * MK_sigma(SS, w + SS%pui_h0_wc) * (SS%pui_wR/10000)/SSS
				end do
				SS%F_integr_pui(i, n2) = w
			end do


			! Считаем частоту перезарядки и источники импульса и энергии
			do i = 1, SS%pui_F_n
				S = 0.0
				SSS = 0.0
				S1 = 0.0
				UH = (i - 0.5) * SS%pui_wR/SS%pui_F_n
				do j = 1, 1000
					w = (j - 0.5) * SS%pui_wR/1000
					ff = PUI_get_f(SS, n2, w)
					do k = 1, 180
						the = par_pi * k/180	
						u = sqrt(w**2 * sin(the)**2 + (w * cos(the) - UH)**2)
						u = max(u, 0.000001_8)
						S = S + 2.0 * par_pi * ff * w**2 * u * MK_sigma(SS, u) * sin(the) * (par_pi/180) * (SS%pui_wR/1000)
						SSS = SSS + 2.0 * par_pi * ff * w**2 * u * MK_sigma(SS, u) * sin(the) * w * cos(the) * (par_pi/180) * (SS%pui_wR/1000)
						S1 = S1 + par_pi * ff * w**2 * u * MK_sigma(SS, u) * sin(the) * w**2 * (par_pi/180) * (SS%pui_wR/1000)   !TODO 
					end do
				end do
				SS%nu_integr_pui(i, n2) = S
				SS%Mz_integr_pui(i, n2) = SSS
				SS%E_integr_pui(i, n2) = S1
			end do


		end do
		!$omp end do
		!$omp end parallel
		print*, "PUI_F_integr_Culc  END"
	end subroutine PUI_F_integr_Culc

	subroutine PUI_n_T_culc(SS)
		TYPE (Setka), intent(in out) :: SS
		integer :: n2, N, i
		real(8) :: w, S, S2

		N = size(SS%f_pui(1, :))

		!$omp parallel
	 	!$omp do private(i, w, S, S2)
		do n2 = 1, N
			S = 0.0
			S2 = 0.0
			do i = 1, SS%pui_nW
				w = (i-0.5) * SS%pui_wR/SS%pui_nW
				S = S + SS%f_pui(i, n2) * 4 * par_pi * w**2 * (SS%pui_wR/SS%pui_nW)
				S2 = S2 + SS%f_pui(i, n2) * 4 * par_pi * w**4 * (SS%pui_wR/SS%pui_nW)
			end do
			S2 = S2/(S * 3.0)
			
			SS%par_pui(1, n2) = S
			SS%par_pui(2, n2) = S2
		end do
		!$omp end do
		!$omp end parallel

	end subroutine PUI_n_T_culc

	real(8) pure function PUI_get_nu_integr(SS, n, w)
		TYPE (Setka), intent(in) :: SS
		!real(8) function PUI_get_f(n, w)
		! Получить значение функции распределения PUI с номером n, со скоростью w (используется линейная интерполяция)
		real(8), intent (in) :: w
		integer, intent (in) :: n
		integer :: k1, k2
		real(8) :: x1, x2,  f1, f2
		k1 = max(min(INT(w/SS%pui_wR * SS%pui_F_n + 0.5), SS%pui_F_n), 1)
		k2 = min(k1 + 1, SS%pui_F_n)
		x1 = (k1 - 0.5) * SS%pui_wR/SS%pui_F_n
		x2 = (k2 - 0.5) * SS%pui_wR/SS%pui_F_n
		f1 = SS%nu_integr_pui(k1, n)
		f2 = SS%nu_integr_pui(k2, n)
		!PUI_get_f = (f1 - f2)/(x1 - x2) * x + (f1 * x2 - f2 * x1)/(x2 - x1)

		! if(w > 5.0 .and.(w < x1 .or. w > x2)) then
		! 	print*, "ERROR PUI_get_f 390", w, x1, x2
		! end if

		if(w < x1 .or. dabs(x1 - x2) <= 0.0000001) then
			PUI_get_nu_integr = f1
		else if (w > x2) then
			PUI_get_nu_integr = f2
		else
			PUI_get_nu_integr = f1 * (w - x2)/(x1 - x2) + f2 * (w - x1)/(x2 - x1)
		end if
    end function PUI_get_nu_integr

	real(8) pure function PUI_get_Mz_integr(SS, n, w)
		TYPE (Setka), intent(in) :: SS
		!real(8) function PUI_get_f(n, w)
		! Получить значение функции распределения PUI с номером n, со скоростью w (используется линейная интерполяция)
		real(8), intent (in) :: w
		integer, intent (in) :: n
		integer :: k1, k2
		real(8) :: x1, x2,  f1, f2
		k1 = max(min(INT(w/SS%pui_wR * SS%pui_F_n + 0.5), SS%pui_F_n), 1)
		k2 = min(k1 + 1, SS%pui_F_n)
		x1 = (k1 - 0.5) * SS%pui_wR/SS%pui_F_n
		x2 = (k2 - 0.5) * SS%pui_wR/SS%pui_F_n
		f1 = SS%Mz_integr_pui(k1, n)
		f2 = SS%Mz_integr_pui(k2, n)

		if(w < x1 .or. dabs(x1 - x2) <= 0.0000001) then
			PUI_get_Mz_integr = f1
		else if (w > x2) then
			PUI_get_Mz_integr = f2
		else
			PUI_get_Mz_integr = f1 * (w - x2)/(x1 - x2) + f2 * (w - x1)/(x2 - x1)
		end if
    end function PUI_get_Mz_integr

	real(8) pure function PUI_get_E_integr(SS, n, w)
		TYPE (Setka), intent(in) :: SS
		!real(8) function PUI_get_f(n, w)
		! Получить значение функции распределения PUI с номером n, со скоростью w (используется линейная интерполяция)
		real(8), intent (in) :: w
		integer, intent (in) :: n
		integer :: k1, k2
		real(8) :: x1, x2,  f1, f2
		k1 = max(min(INT(w/SS%pui_wR * SS%pui_F_n + 0.5), SS%pui_F_n), 1)
		k2 = min(k1 + 1, SS%pui_F_n)
		x1 = (k1 - 0.5) * SS%pui_wR/SS%pui_F_n
		x2 = (k2 - 0.5) * SS%pui_wR/SS%pui_F_n
		f1 = SS%E_integr_pui(k1, n)
		f2 = SS%E_integr_pui(k2, n)
		!PUI_get_f = (f1 - f2)/(x1 - x2) * x + (f1 * x2 - f2 * x1)/(x2 - x1)

		! if(w > 5.0 .and.(w < x1 .or. w > x2)) then
		! 	print*, "ERROR PUI_get_f 390", w, x1, x2
		! end if

		if(w < x1 .or. dabs(x1 - x2) <= 0.0000001) then
			PUI_get_E_integr = f1
		else if (w > x2) then
			PUI_get_E_integr = f2
		else
			PUI_get_E_integr = f1 * (w - x2)/(x1 - x2) + f2 * (w - x1)/(x2 - x1)
		end if
    end function PUI_get_E_integr

end module PUI