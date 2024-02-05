module Interpol
    use STORAGE 
    use My_func
    USE ieee_arithmetic
    implicit none 

    contains

    subroutine Int_Init(SS_int, SS)
        TYPE (Setka), intent(in) :: SS
        TYPE (Inter_Setka), intent(in out) :: SS_int

        integer(4) :: N1, N2, i, cell, N3, M1, M2, M3, j, it, nn, yzel, a1, a2, cell2, n_yzel
        real(8) :: p1(2)

        if(SS%init_geo == .False.) STOP "ERROR  Int_Init  11 987r4yg3hjfkgf34f43"
        if(SS_int%init == .True.) STOP "ERROR  Int_Init  12 12edrc4tveyrbunimo98,0p78i6n7ub6y"

        SS_int%init = .True.

        SS_int%n_Hidrogen = SS%n_Hidrogen
        SS_int%n_par = SS%n_par

        N1 = size(SS%gl_Cell_A(:, 1))
        N2 = size(SS%gl_Cell_A(1, :))
        allocate(SS_int%gl_Cell_A(N1, N2 + 1))

        N1 = size(SS%gl_Cell_B(:, 1))
        N2 = size(SS%gl_Cell_B(1, :))
        allocate(SS_int%gl_Cell_B(N1, N2 + 1))

        N1 = size(SS%gl_Cell_C(:, 1))
        N2 = size(SS%gl_Cell_C(1, :))
        allocate(SS_int%gl_Cell_C(N1 - 1, N2 - 1))

        SS_int%gl_Cell_A = 1
        SS_int%gl_Cell_B = 1
        SS_int%gl_Cell_C = 1

        N1 = size(SS_int%gl_Cell_A) + size(SS_int%gl_Cell_B) + size(SS_int%gl_Cell_C) - SS%par_n_TS + 1
        allocate(SS_int%gl_all_Cell(4, N1))
        allocate(SS_int%gl_all_triangle(N1))
        allocate(SS_int%gl_Cell_center(2, N1))
        allocate(SS_int%gl_Cell_neighbour(4, N1))
        allocate(SS_int%gl_Cell_Belong(3, 4, N1))
        allocate(SS_int%gl_Cell_interpol_matrix(4, 4, N1))
		
		SS_int%gl_all_Cell = 0
        SS_int%gl_Cell_Belong = 0.0
        SS_int%gl_Cell_neighbour = 0
        SS_int%gl_Cell_center = 0.0
        SS_int%gl_Cell_interpol_matrix = 0.0
        SS_int%gl_all_triangle = .False.

        N1 = size(SS%gl_all_Cell(1, :))
        allocate(SS_int%gl_yzel(2, N1 + SS%par_n_END - 1 + SS%par_n_TS + SS%par_m_O - 1 + 1))  ! Добавляем центр координат
		SS_int%gl_yzel = 0.0
        n_yzel = N1 + SS%par_n_END - 1 + SS%par_n_TS + SS%par_m_O - 1 + 1
        allocate(SS_int%gd(SS_int%n_par, n_yzel))
        allocate(SS_int%hydrogen(5, SS_int%n_Hidrogen, n_yzel))

        allocate(SS_int%atom_all_source(4, SS_int%n_Hidrogen, n_yzel))
        allocate(SS_int%atom_source(7, n_yzel))


		
        ! Заполняем координаты узлов --------------------------------------------------------------------------
		do i = 1, N1
			SS_int%gl_yzel(:, i) = SS%gl_Cell_Centr(:, i, 1)
			SS_int%gd(:, i) = SS%gd(:, i, 1)
			SS_int%hydrogen(:, :, i) = SS%hydrogen(:, :, i, 1)
			SS_int%atom_all_source(:, :, i) = SS%atom_all_source(:, :, i)
			SS_int%atom_source(:, i) = SS%atom_source(:, i)
		end do

        N2 = size(SS%gl_Cell_A(:, 1))

        do i = 1, N2
            cell = SS%gl_Cell_A(i, 1)
            p1 = SS%gl_Cell_Centr(:, cell, 1)
            p1(1) = norm2(p1)
            p1(2) = 0.0  ! Проекция узла на ось симметрии
            SS_int%gl_yzel(:, N1 + i) = p1

            SS_int%gd(:, N1 + i) = SS%gd(:, cell, 1)
			SS_int%hydrogen(:, :, N1 + i) = SS%hydrogen(:, :, cell, 1)
			SS_int%atom_all_source(:, :, N1 + i) = SS%atom_all_source(:, :, cell)
			SS_int%atom_source(:, N1 + i) = SS%atom_source(:, cell)

            SS_int%gd(4, N1 + i) = 0.0
            SS_int%hydrogen(4, :, N1 + i) = 0.0
            SS_int%atom_source(2, N1 + i) = 0.0
		end do
		
		N3 = size(SS%gl_Cell_B(:, 1))
		
		do i = 1, N3
            cell = SS%gl_Cell_B(i, 1)
            p1 = SS%gl_Cell_Centr(:, cell, 1)
            p1(1) = -norm2(p1)
            p1(2) = 0.0  ! Проекция узла на ось симметрии
            SS_int%gl_yzel(:, N1 + N2 + i) = p1

            SS_int%gd(:,N1 + N2 + i) = SS%gd(:, cell, 1)
			SS_int%hydrogen(:, :, N1 + N2 + i) = SS%hydrogen(:, :, cell, 1)
            SS_int%atom_all_source(:, :, N1 + N2 + i) = SS%atom_all_source(:, :, cell)
			SS_int%atom_source(:, N1 + N2 + i) = SS%atom_source(:, cell)

            SS_int%gd(4, N1 + N2 + i) = 0.0
            SS_int%hydrogen(4, :, N1 + N2 + i) = 0.0
            SS_int%atom_source(2, N1 + N2 + i) = 0.0
		end do
		
		! Заполняем A - ячейки --------------------------------------------------------------------------
		
		M1 = size(SS_int%gl_Cell_A(:, 1))
        M2 = size(SS_int%gl_Cell_A(1, :))
		
		it = 1
		
		do j = 1, M2
			do i = 1, M1
				SS_int%gl_Cell_A(i, j) = it
				
				if(j == M2 .and. i < SS%par_n_TS)then
					nn = size(SS_int%gl_Cell_B(1, :))
					SS_int%gl_Cell_B(i, nn) = it
				end if

                if(j == 1) then

                    if( i /= 1) then
                        cell = N1 + i - 1
                    else
                        cell = N1 + N2 + N3 + 1
                    end if
                    SS_int%gl_all_Cell(1, it) = cell

                    cell = N1 + i
                    SS_int%gl_all_Cell(2, it) = cell

                    cell = SS%gl_Cell_A(i, j)
                    SS_int%gl_all_Cell(3, it) = cell

                    if( i /= 1) then
                        cell = SS%gl_Cell_A(i - 1, j)
                    else
                        cell = N1 + N2 + N3 + 1
                    end if
                    SS_int%gl_all_Cell(4, it) = cell
                else if (j < M2) then

                    if( i /= 1) then
                        cell = SS%gl_Cell_A(i - 1, j - 1)
                    else
                        cell = N1 + N2 + N3 + 1
                    end if
                    SS_int%gl_all_Cell(1, it) = cell

                    cell = SS%gl_Cell_A(i, j - 1)
                    SS_int%gl_all_Cell(2, it) = cell

                    cell = SS%gl_Cell_A(i, j)
                    SS_int%gl_all_Cell(3, it) = cell

                    if( i /= 1) then
                        cell = SS%gl_Cell_A(i - 1, j)
                    else
                        cell = N1 + N2 + N3 + 1
                    end if
                    SS_int%gl_all_Cell(4, it) = cell
				else
					if( i /= 1) then
                        cell = SS%gl_Cell_A(i - 1, j - 1)
                    else
                        cell = N1 + N2 + N3 + 1
                    end if
                    SS_int%gl_all_Cell(1, it) = cell

                    cell = SS%gl_Cell_A(i, j - 1)
                    SS_int%gl_all_Cell(2, it) = cell

					if(i < SS%par_n_TS) then
                        nn = size(SS%gl_Cell_B(1, :))
                        cell = SS%gl_Cell_B(i, nn)
                    else
                        cell = SS%gl_Cell_C(i - SS%par_n_TS + 1, 1)
                    end if
                    
                    SS_int%gl_all_Cell(3, it) = cell

                    if( i /= 1) then
                        if(i < SS%par_n_TS) then
                            nn = size(SS%gl_Cell_B(1, :))
                            cell = SS%gl_Cell_B(i - 1, nn)
                        else
                            cell = SS%gl_Cell_C( max(i - SS%par_n_TS, 1), 1)
                        end if
                    else
                        cell = N1 + N2 + N3 + 1
                    end if

                    SS_int%gl_all_Cell(4, it) = cell
                end if

				it = it + 1
			end do
		end do

        ! Заполняем B - ячейки --------------------------------------------------------------------------
		
		M1 = size(SS_int%gl_Cell_B(:, 1))
        M2 = size(SS_int%gl_Cell_B(1, :))
		
		do j = 1, M2
			do i = 1, M1
				if(j == M2 .and. i < SS%par_n_TS) then
					CYCLE
				end if
				
				SS_int%gl_Cell_B(i, j) = it

                if(j == 1) then

                    if( i /= 1) then
                        cell = N1 + N2 +  i - 1
                    else
                        cell = N1 + N2 + N3 + 1
                    end if
                    SS_int%gl_all_Cell(1, it) = cell

                    cell = N1 + N2 + i
                    SS_int%gl_all_Cell(2, it) = cell

                    cell = SS%gl_Cell_B(i, j)
                    SS_int%gl_all_Cell(3, it) = cell

                    if( i /= 1) then
                        cell = SS%gl_Cell_B(i - 1, j)
                    else
                        cell = N1 + N2 + N3 + 1
                    end if
                    SS_int%gl_all_Cell(4, it) = cell
				else if (j < M2) then
					
                    if( i /= 1) then
                        cell = SS%gl_Cell_B(i - 1, j - 1)
                    else
                        cell = N1 + N2 + N3 + 1
                    end if
                    SS_int%gl_all_Cell(1, it) = cell

                    cell = SS%gl_Cell_B(i, j - 1)
                    SS_int%gl_all_Cell(2, it) = cell

                    cell = SS%gl_Cell_B(i, j)
                    SS_int%gl_all_Cell(3, it) = cell

                    if( i /= 1) then
                        cell = SS%gl_Cell_B(i - 1, j)
                    else
                        cell = N1 + N2 + N3 + 1
                    end if
                    SS_int%gl_all_Cell(4, it) = cell
				else
					if(i == SS%par_n_TS) then
						cell = SS%gl_Cell_B(i - 1, j - 1)
                        SS_int%gl_all_Cell(1, it) = cell

                        cell = SS%gl_Cell_B(i, j - 1)
                        SS_int%gl_all_Cell(2, it) = cell

                        cell = SS%gl_Cell_C(1, i - SS%par_n_TS + 1)
                        SS_int%gl_all_Cell(3, it) = cell
                        
						nn = size(SS%gl_Cell_A(SS%par_n_TS, :))
					    cell = SS%gl_Cell_A(SS%par_n_TS - 1, nn)
                        SS_int%gl_all_Cell(4, it) = cell
					else
					    cell = SS%gl_Cell_B(i - 1, j - 1)
                        SS_int%gl_all_Cell(1, it) = cell

                        cell = SS%gl_Cell_B(i, j - 1)
                        SS_int%gl_all_Cell(2, it) = cell

                        cell = SS%gl_Cell_C(1, i - SS%par_n_TS + 1)
                        SS_int%gl_all_Cell(3, it) = cell

					    cell = SS%gl_Cell_C(1, i - SS%par_n_TS)
                        SS_int%gl_all_Cell(4, it) = cell
					end if
                end if

				it = it + 1
			end do
		end do
		
        ! Заполняем C - ячейки --------------------------------------------------------------------------
		
		M1 = size(SS_int%gl_Cell_C(:, 1))
        M2 = size(SS_int%gl_Cell_C(1, :))
		
		do j = 1, M2
			do i = 1, M1
				SS_int%gl_Cell_C(i, j) = it

                cell = SS%gl_Cell_C(i, j)
                SS_int%gl_all_Cell(1, it) = cell

                cell = SS%gl_Cell_C(i + 1, j)
                SS_int%gl_all_Cell(2, it) = cell

                cell = SS%gl_Cell_C(i + 1, j + 1)
                SS_int%gl_all_Cell(3, it) = cell

                cell = SS%gl_Cell_C(i, j + 1)
                SS_int%gl_all_Cell(4, it) = cell
               

				it = it + 1
			end do
		end do

        ! Посчитаем центры ячеек
        M1 = size(SS_int%gl_Cell_center(1, :))
        do i = 1, M1
            p1 = 0.0
            do j = 1, 4
                yzel = SS_int%gl_all_Cell(j, i)
                p1 = p1 + SS_int%gl_yzel(:, yzel)
            end do
            SS_int%gl_Cell_center(:, i) = p1/4.0
		end do

		! Теперь необходимо заполнить соседей (для быстрого поиска по сетке)

		! А - ячейки
        M1 = size(SS_int%gl_Cell_A(:, 1))
        M2 = size(SS_int%gl_Cell_A(1, :))

        do j = 1, M2
			do i = 1, M1 - 1
                cell = SS_int%gl_Cell_A(i, j)
                a1 = SS_int%gl_all_Cell(2, cell)
                a2 = SS_int%gl_all_Cell(3, cell)

                cell2 = SS_int%gl_Cell_A(i + 1, j)
                if(cell2 < 1) CYCLE

                if ( ANY(SS_int%gl_all_Cell(:, cell2) == a1) .and. ANY(SS_int%gl_all_Cell(:, cell2) == a2) ) then
                    SS_int%gl_Cell_neighbour(2, cell) = cell2
                end if
            end do
        end do

        do j = 2, M2
			do i = 1, M1
                cell = SS_int%gl_Cell_A(i, j)
                a1 = SS_int%gl_all_Cell(1, cell)
                a2 = SS_int%gl_all_Cell(2, cell)

                cell2 = SS_int%gl_Cell_A(i, j - 1)
                if(cell2 < 1) CYCLE

                if ( ANY(SS_int%gl_all_Cell(:, cell2) == a1) .and. ANY(SS_int%gl_all_Cell(:, cell2) == a2) ) then
                    SS_int%gl_Cell_neighbour(1, cell) = cell2
                end if
            end do
		end do
		
		do j = 1, M2 - 1
			do i = 1, M1
                cell = SS_int%gl_Cell_A(i, j)
                a1 = SS_int%gl_all_Cell(3, cell)
                a2 = SS_int%gl_all_Cell(4, cell)

                cell2 = SS_int%gl_Cell_A(i, j + 1)
                if(cell2 < 1) CYCLE

                if ( ANY(SS_int%gl_all_Cell(:, cell2) == a1) .and. ANY(SS_int%gl_all_Cell(:, cell2) == a2) ) then
                    SS_int%gl_Cell_neighbour(3, cell) = cell2
                end if
            end do
		end do
        
		do j = 1, M2
			do i = 2, M1
                cell = SS_int%gl_Cell_A(i, j)
                a1 = SS_int%gl_all_Cell(4, cell)
                a2 = SS_int%gl_all_Cell(1, cell)

                cell2 = SS_int%gl_Cell_A(i - 1, j)
                if(cell2 < 1) CYCLE

                if ( ANY(SS_int%gl_all_Cell(:, cell2) == a1) .and. ANY(SS_int%gl_all_Cell(:, cell2) == a2) ) then
                    SS_int%gl_Cell_neighbour(4, cell) = cell2
                end if
            end do
		end do

        ! B - ячейки
        M1 = size(SS_int%gl_Cell_B(:, 1))
        M2 = size(SS_int%gl_Cell_B(1, :))

        do j = 1, M2
			do i = 1, M1 - 1
                cell = SS_int%gl_Cell_B(i, j)
                a1 = SS_int%gl_all_Cell(2, cell)
                a2 = SS_int%gl_all_Cell(3, cell)

                cell2 = SS_int%gl_Cell_B(i + 1, j)
                if(cell2 < 1) CYCLE

                if ( ANY(SS_int%gl_all_Cell(:, cell2) == a1) .and. ANY(SS_int%gl_all_Cell(:, cell2) == a2) ) then
                    SS_int%gl_Cell_neighbour(2, cell) = cell2
                end if
            end do
        end do

        do j = 2, M2
			do i = 1, M1
                cell = SS_int%gl_Cell_B(i, j)
                a1 = SS_int%gl_all_Cell(1, cell)
                a2 = SS_int%gl_all_Cell(2, cell)

                cell2 = SS_int%gl_Cell_B(i, j - 1)
                if(cell2 < 1) CYCLE

                if ( ANY(SS_int%gl_all_Cell(:, cell2) == a1) .and. ANY(SS_int%gl_all_Cell(:, cell2) == a2) ) then
                    SS_int%gl_Cell_neighbour(1, cell) = cell2
                end if
            end do
		end do
		
		do j = 1, M2 - 1
			do i = 1, M1
                cell = SS_int%gl_Cell_B(i, j)
                a1 = SS_int%gl_all_Cell(3, cell)
                a2 = SS_int%gl_all_Cell(4, cell)

                cell2 = SS_int%gl_Cell_B(i, j + 1)
                if(cell2 < 1) CYCLE

                if ( ANY(SS_int%gl_all_Cell(:, cell2) == a1) .and. ANY(SS_int%gl_all_Cell(:, cell2) == a2) ) then
                    SS_int%gl_Cell_neighbour(3, cell) = cell2
                end if
            end do
		end do
        
		do j = 1, M2
			do i = 2, M1
                cell = SS_int%gl_Cell_B(i, j)
                a1 = SS_int%gl_all_Cell(4, cell)
                a2 = SS_int%gl_all_Cell(1, cell)

                cell2 = SS_int%gl_Cell_B(i - 1, j)
                if(cell2 < 1) CYCLE

                if ( ANY(SS_int%gl_all_Cell(:, cell2) == a1) .and. ANY(SS_int%gl_all_Cell(:, cell2) == a2) ) then
                    SS_int%gl_Cell_neighbour(4, cell) = cell2
                end if
            end do
		end do
		
		! C - ячейки
		M1 = size(SS_int%gl_Cell_C(:, 1))
        M2 = size(SS_int%gl_Cell_C(1, :))
		
		do j = 2, M2
			do i = 1, M1
                cell = SS_int%gl_Cell_C(i, j)
                a1 = SS_int%gl_all_Cell(1, cell)
                a2 = SS_int%gl_all_Cell(2, cell)

                cell2 = SS_int%gl_Cell_C(i, j - 1)
                if(cell2 < 1) CYCLE

                if ( ANY(SS_int%gl_all_Cell(:, cell2) == a1) .and. ANY(SS_int%gl_all_Cell(:, cell2) == a2) ) then
                    SS_int%gl_Cell_neighbour(1, cell) = cell2
                end if
            end do
		end do
		
		do j = 1, M2
			do i = 1, M1 - 1
                cell = SS_int%gl_Cell_C(i, j)
                a1 = SS_int%gl_all_Cell(2, cell)
                a2 = SS_int%gl_all_Cell(3, cell)

                cell2 = SS_int%gl_Cell_C(i + 1, j)
                if(cell2 < 1) CYCLE

                if ( ANY(SS_int%gl_all_Cell(:, cell2) == a1) .and. ANY(SS_int%gl_all_Cell(:, cell2) == a2) ) then
                    SS_int%gl_Cell_neighbour(2, cell) = cell2
                end if
            end do
		end do
		
		do j = 1, M2 - 1
			do i = 1, M1
                cell = SS_int%gl_Cell_C(i, j)
                a1 = SS_int%gl_all_Cell(3, cell)
                a2 = SS_int%gl_all_Cell(4, cell)

                cell2 = SS_int%gl_Cell_C(i, j + 1)
                if(cell2 < 1) CYCLE

                if ( ANY(SS_int%gl_all_Cell(:, cell2) == a1) .and. ANY(SS_int%gl_all_Cell(:, cell2) == a2) ) then
                    SS_int%gl_Cell_neighbour(3, cell) = cell2
                end if
            end do
		end do
		
		do j = 1, M2
			do i = 2, M1
                cell = SS_int%gl_Cell_C(i, j)
                a1 = SS_int%gl_all_Cell(4, cell)
                a2 = SS_int%gl_all_Cell(1, cell)

                cell2 = SS_int%gl_Cell_C(i - 1, j)
                if(cell2 < 1) CYCLE

                if ( ANY(SS_int%gl_all_Cell(:, cell2) == a1) .and. ANY(SS_int%gl_all_Cell(:, cell2) == a2) ) then
                    SS_int%gl_Cell_neighbour(4, cell) = cell2
                end if
            end do
		end do

        ! Связи мкжду A-C ячейками

        M1 = size(SS_int%gl_Cell_A(:, 1))
        M2 = size(SS_int%gl_Cell_A(1, :))

        do j = M2, M2
			do i = SS%par_n_TS + 1, M1
                cell = SS_int%gl_Cell_A(i, j)
                a1 = SS_int%gl_all_Cell(3, cell)
                a2 = SS_int%gl_all_Cell(4, cell)

                cell2 = SS_int%gl_Cell_C(i - SS%par_n_TS, 1)
                if(cell2 < 1) CYCLE

                if ( ANY(SS_int%gl_all_Cell(:, cell2) == a1) .and. ANY(SS_int%gl_all_Cell(:, cell2) == a2) ) then
                    SS_int%gl_Cell_neighbour(3, cell) = cell2
                end if
            end do
		end do

        ! Связи мкжду A-B ячейками

        do j = M2, M2
			do i = SS%par_n_TS, SS%par_n_TS
                cell = SS_int%gl_Cell_A(i, j)
                a1 = SS_int%gl_all_Cell(3, cell)
                a2 = SS_int%gl_all_Cell(4, cell)

                cell2 = SS_int%gl_Cell_B(i, size(SS_int%gl_Cell_B(1, :)))
                if(cell2 < 1) CYCLE

                if ( ANY(SS_int%gl_all_Cell(:, cell2) == a1) .and. ANY(SS_int%gl_all_Cell(:, cell2) == a2) ) then
                    SS_int%gl_Cell_neighbour(3, cell) = cell2
                end if
            end do
		end do

        do j = M2, M2
			do i = 1, SS%par_n_TS
                cell = SS_int%gl_Cell_A(i, j)
                a1 = SS_int%gl_all_Cell(3, cell)
                a2 = SS_int%gl_all_Cell(4, cell)

                cell2 = SS_int%gl_Cell_B(i, size(SS_int%gl_Cell_B(1, :)) - 1)
                if(cell2 < 1) CYCLE

                if ( ANY(SS_int%gl_all_Cell(:, cell2) == a1) .and. ANY(SS_int%gl_all_Cell(:, cell2) == a2) ) then
                    SS_int%gl_Cell_neighbour(3, cell) = cell2
                end if
            end do
		end do

        ! Связи мкжду C-A ячейками

        M1 = size(SS_int%gl_Cell_C(:, 1))
        M2 = size(SS_int%gl_Cell_C(1, :))
		
		do j = 1, 1
			do i = 1, M1
                cell = SS_int%gl_Cell_C(i, j)
                a1 = SS_int%gl_all_Cell(1, cell)
                a2 = SS_int%gl_all_Cell(2, cell)

                cell2 = SS_int%gl_Cell_A(i + SS%par_n_TS, size(SS_int%gl_Cell_A(1, :)))
                if(cell2 < 1) CYCLE

                if ( ANY(SS_int%gl_all_Cell(:, cell2) == a1) .and. ANY(SS_int%gl_all_Cell(:, cell2) == a2) ) then
                    SS_int%gl_Cell_neighbour(1, cell) = cell2
                end if
            end do
		end do

        ! Связи мкжду C-B ячейками

        M1 = size(SS_int%gl_Cell_C(:, 1))
        M2 = size(SS_int%gl_Cell_C(1, :))

        do j = 1, M2
			do i = 1, 1
                cell = SS_int%gl_Cell_C(i, j)
                a1 = SS_int%gl_all_Cell(4, cell)
                a2 = SS_int%gl_all_Cell(1, cell)

                cell2 = SS_int%gl_Cell_B(j + SS%par_n_TS, size(SS_int%gl_Cell_B(1, :)))
                if(cell2 < 1) CYCLE

                if ( ANY(SS_int%gl_all_Cell(:, cell2) == a1) .and. ANY(SS_int%gl_all_Cell(:, cell2) == a2) ) then
                    SS_int%gl_Cell_neighbour(4, cell) = cell2
                end if
            end do
		end do

        ! Связи мкжду B-C ячейками

        M1 = size(SS_int%gl_Cell_B(:, 1))
        M2 = size(SS_int%gl_Cell_B(1, :))

        do j = M2, M2
			do i = SS%par_n_TS + 1, M1
                cell = SS_int%gl_Cell_B(i, j)
                a1 = SS_int%gl_all_Cell(3, cell)
                a2 = SS_int%gl_all_Cell(4, cell)

                cell2 = SS_int%gl_Cell_C(1, i - SS%par_n_TS)
                if(cell2 < 1) CYCLE

                if ( ANY(SS_int%gl_all_Cell(:, cell2) == a1) .and. ANY(SS_int%gl_all_Cell(:, cell2) == a2) ) then
                    SS_int%gl_Cell_neighbour(3, cell) = cell2
                end if
            end do
		end do

        ! B-A тройная точка

        M1 = size(SS_int%gl_Cell_B(:, 1))
        M2 = size(SS_int%gl_Cell_B(1, :))

        do j = M2, M2
			do i = SS%par_n_TS, SS%par_n_TS
                cell = SS_int%gl_Cell_B(i, j)
                a1 = SS_int%gl_all_Cell(3, cell)
                a2 = SS_int%gl_all_Cell(4, cell)

                cell2 = SS_int%gl_Cell_A(SS%par_n_TS, size(SS_int%gl_Cell_A(1, :)))
                if(cell2 < 1) CYCLE

                if ( ANY(SS_int%gl_all_Cell(:, cell2) == a1) .and. ANY(SS_int%gl_all_Cell(:, cell2) == a2) ) then
                    SS_int%gl_Cell_neighbour(3, cell) = cell2
                end if
            end do
		end do

        ! Считаем Ax + By + C для каждой грани
        call Int_Belong_init(SS_int)
        call Int_Matrix_init(SS_int)

    end subroutine Int_Init

    subroutine Int_Matrix_init(SS)
	    USE ieee_arithmetic
        TYPE (Inter_Setka), intent(in out) :: SS
        integer(4) :: i, N, a1, a2, a3, a4
        real(8) :: M(4, 4), p1(2), p2(2), p3(2), p4(2), MM(3, 3), S

        N = size(SS%gl_all_Cell(1, :))

        do i = 1, N
            a1 = SS%gl_all_Cell(1, i)
            a2 = SS%gl_all_Cell(2, i)
            a3 = SS%gl_all_Cell(3, i)
            a4 = SS%gl_all_Cell(4, i)

            p1 = SS%gl_yzel(:, a1)
            p2 = SS%gl_yzel(:, a2)
            p3 = SS%gl_yzel(:, a3)
            p4 = SS%gl_yzel(:, a4)

            if(a1 == a4 .or. a3 == a4) then
                ! В этом случае у нас треугольник
                SS%gl_all_triangle(i) = .True.
                MM(1, 1) = p1(1)
                MM(1, 2) = p1(2)
                MM(1, 3) = 1.0

                MM(2, 1) = p2(1)
                MM(2, 2) = p2(2)
                MM(2, 3) = 1.0

                MM(3, 1) = p3(1)
                MM(3, 2) = p3(2)
                MM(3, 3) = 1.0


                SS%gl_Cell_interpol_matrix(1:3, 1:3, i) = matinv3(MM)

                S = SS%gl_Cell_interpol_matrix(1, 1, i)

                if(ieee_is_normal(S) == .False.) STOP "Error 690 Int_Matrix_init vdbfhf46397kjhgefjk99"
            else
                SS%gl_all_triangle(i) = .False.
                M(1, 1) = p1(1)
                M(1, 2) = p1(2)
                M(1, 3) = p1(1) * p1(2)
                M(1, 4) = 1.0

                M(2, 1) = p2(1)
                M(2, 2) = p2(2)
                M(2, 3) = p2(1) * p2(2)
                M(2, 4) = 1.0

                M(3, 1) = p3(1)
                M(3, 2) = p3(2)
                M(3, 3) = p3(1) * p3(2)
                M(3, 4) = 1.0

                M(4, 1) = p4(1)
                M(4, 2) = p4(2)
                M(4, 3) = p4(1) * p4(2)
                M(4, 4) = 1.0

                SS%gl_Cell_interpol_matrix(:, :, i) = matinv4(M)
                S = SS%gl_Cell_interpol_matrix(1, 1, i)
                if(ieee_is_normal(S) == .False.) STOP "Error 690 Int_Matrix_init vdbfhf46397kjhgefjk99"
            end if
        end do

    end subroutine Int_Matrix_init

    subroutine Int_Belong_init(SS)
        TYPE (Inter_Setka), intent(in out) :: SS
        integer(4) :: N1, i, j, gran, yz1, yz2, jj
        real(8) :: c, p1(2), p2(2), n(2), centr(2)

        N1 = size(SS%gl_Cell_Belong(1, 1, :))

        do i = 1, N1  ! По всем ячейкам

            do j = 1, 4  ! По всем граням
                yz1 = SS%gl_all_Cell(j, i)
                jj = j + 1
                if(jj > 4) jj = 1
                yz2 = SS%gl_all_Cell(jj, i)

                if(yz1 == yz2) then
                    SS%gl_Cell_Belong(:, j, i) = 0.0
                    CYCLE
                end if

                p1 = SS%gl_yzel(:, yz1)
                p2 = SS%gl_yzel(:, yz2)


                n = p2 - p1
                c = n(1)
                n(1) = -n(2)
                n(2) = c

                c = -DOT_PRODUCT(p1, n)

                SS%gl_Cell_Belong(1, j, i) = n(1)
                SS%gl_Cell_Belong(2, j, i) = n(2)
                SS%gl_Cell_Belong(3, j, i) = c
                centr = SS%gl_Cell_center(:, i)

                if(centr(1) * n(1) + centr(2) * n(2) + c > 0) then
                    SS%gl_Cell_Belong(:, j, i) = -SS%gl_Cell_Belong(:, j, i)
                end if
            end do
        end do
    end subroutine Int_Belong_init

    subroutine Int_Find_Cell(SS, xx, yy, num, outer)
        ! Поиск номера ячейки по её координатам
        ! num = предположительный изначальный номер (если не знаем пусть будет равен 1)
        ! Если outer == .True. то точка находится за пределами расчётной области, но найдена максимально близкая к ней ячейка
        ! В этом случае просто не надо интерполировать, а надо просто взять значения в блихжайшём узле
        TYPE (Inter_Setka), intent(in) :: SS
        real(8), intent(in) :: xx, yy
        integer(4), intent(in out) :: num
        integer(4) :: j, gran, sosed, max_num
        LOGICAL, intent(out) :: outer
        real(8) :: x, y

        x = xx
        y = yy

        max_num = 0
        outer = .False.

        loop1:do while(.TRUE.)
            max_num = max_num + 1
            outer = .False.

            if(max_num > 3000) then
                x = x + 0.000001
                y = y + 0.000001
            end if

            if(max_num > 100000) then
                print*, x, y
                print*, num
                STOP "ERROR Int_Find_Cell 683 iuyhgw4it0pflkjhrfe"
            end if

            !print*, "A", num
            
            loop2:do j = 1, 4
                sosed = SS%gl_Cell_neighbour(j, num)

                !print*, "B", j, sosed
                !print*, "C", SS%gl_Cell_Belong(:, j, num)

                if(SS%gl_Cell_Belong(1, j, num) * x + SS%gl_Cell_Belong(2, j, num) * y + SS%gl_Cell_Belong(3, j, num) > 0) then
                    if(sosed == 0) then
                        outer = .True.
                        CYCLE loop2
                    end if
                    num = sosed
                    cycle loop1
                end if
            end do loop2

            !print*, "Stop", outer

            ! Если outer == .True. то точка находится за пределами расчётной области, но найдена максимально близкая к ней ячейка
            ! В этом случае просто не надо интерполировать, а надо просто взять значения в центре ячейки

            EXIT loop1
        end do loop1
    end subroutine Int_Find_Cell

    subroutine Int_Belong_Cell(SS, xx, yy, num, outer)
        TYPE (Inter_Setka), intent(in) :: SS
        real(8), intent(in) :: xx, yy
        integer(4), intent(in out) :: num
        integer(4) :: j, gran, sosed
        LOGICAL, intent(out) :: outer
        real(8) :: x, y

        x = xx
        y = yy
        outer = .True.

        loop2:do j = 1, 4
            if(SS%gl_Cell_Belong(1, j, num) * x + SS%gl_Cell_Belong(2, j, num) * y + SS%gl_Cell_Belong(3, j, num) > 0) then
                outer = .False.
                EXIT loop2
            end if
        end do loop2

    end subroutine Int_Belong_Cell

    subroutine Int_Get_Parameter(SS, x, y, num, PAR_gd, PAR_hydrogen, PAR_atom_source)
        TYPE (Inter_Setka), intent(in) :: SS
        real(8), intent(in) :: x, y
        integer(4), intent(in out) :: num
        real(8), intent(out), optional :: PAR_gd(:)
        real(8), intent(out), optional :: PAR_hydrogen(:, :)
        real(8), intent(out), optional :: PAR_atom_source(:)

        LOGICAL :: outer, outer2
        real(8) :: dist, r(2), dist_min, yzel_min, M(4, 4), v(4), MM(4), ccc, kk(4), kkk, p1(2), p2(2), p3(2), p4(2), n1(2), n2(2)
        real(8) :: a, b, c, d, c1, c2, x0, y0, d1(2), d2(2), t, e1, e2, e3, e4, l1, l2, f1, f2, a2, b2
        integer(4) :: i, yzel, j, yz1, yz2, yz3, yz4

        if(present(PAR_atom_source)) then
            if(size(PAR_atom_source) /= 4) STOP "ERROR Int 851 y2w0S3KFz4v7npnDl9rgziYnZNfNSc"
        end if

        call Int_Find_Cell(SS, x, y, num, outer)  ! Находим ячейку, которой принадлежит точка

        if(outer == .True.) then  ! Точка за пределами области не надо интерполировать, надо найти близжайший узел
            dist_min = 100000.0
            do i = 1, 4
                yzel = SS%gl_all_Cell(i, num)
                r = SS%gl_yzel(:, yzel)
                dist = sqrt((x - r(1))**2 + (y - r(2))**2)
                if(dist < dist_min) then
                    dist_min = dist
                    yzel_min = yzel
                end if
            end do

            if(present(PAR_gd)) then
                if(size(PAR_gd) /= SS%n_par) STOP "Error Int_Get_Parameter size(PAR_gd) /= SS%n_par   09876tyuinewufhugferafrgnmiutr"
                PAR_gd = SS%gd(:, yzel_min)
            end if

            if(present(PAR_atom_source)) then
                PAR_atom_source = SS%atom_source(1:4, yzel_min)
            end if

            if(present(PAR_hydrogen)) then
                if(size(PAR_hydrogen(1, :)) /= SS%n_Hidrogen) STOP "Error Int_Get_Parameter size(PAR_hydrogen) /= SS%n_par   oehwjfiehurgegvve"
                PAR_hydrogen = SS%hydrogen(:, :, yzel_min)
            end if

        else
            if(SS%gl_all_triangle(num) == .True.) then ! ТРЕУГОЛЬНИК
				
                M = SS%gl_Cell_interpol_matrix(:, :, num)

                MM = 0.0

                v(1) = x
                v(2) = y
                v(3) = 1.0

                do i = 1, 3
                    do j = 1, 3
                        MM(i) = MM(i) + v(j) * M(j, i)
                    end do
                end do

                if(present(PAR_gd)) then
                    if(size(PAR_gd) /= SS%n_par) STOP "Error Int_Get_Parameter size(PAR_gd) /= SS%n_par   09876tyuinewufhugferafrgnmiutr"
                    PAR_gd = 0.0
                    do i = 1, 3
                        yzel = SS%gl_all_Cell(i, num)
                        PAR_gd(:) = PAR_gd(:) + MM(i) * SS%gd(:, yzel)
                    end do

                    ! if(PAR_gd(1) <= 0.0) then
                    !     Print*, "ERROR 887 <0 ouw7n8yvn987rnt97wct5gewv"
                    !     pause
                    !     STOP
                    ! end if
                end if

                if(present(PAR_atom_source)) then
                    PAR_atom_source = 0.0
                    do i = 1, 3
                        yzel = SS%gl_all_Cell(i, num)
                        PAR_atom_source(:) = PAR_atom_source(:) + MM(i) * SS%atom_source(1:4, yzel)
                    end do
                end if

                if(present(PAR_hydrogen)) then
                    if(size(PAR_hydrogen(1, :)) /= SS%n_Hidrogen) STOP "Error Int_Get_Parameter size(PAR_hydrogen) /= SS%n_par   oehwjfiehurgegvve"
                    PAR_hydrogen = 0.0
                    do i = 1, 3
                        yzel = SS%gl_all_Cell(i, num)
                        PAR_hydrogen = PAR_hydrogen + MM(i) * SS%hydrogen(:, :, yzel)
                    end do

                    ! if(PAR_hydrogen(1, 4) > 3.0) then
                    !     Print*, "ERROR 887 <0 y54b6unewvbyerev5e"
                    !     pause
                    !     STOP
                    ! end if

                end if


			else   ! ПРЯМОУГОЛЬНИК
				
                do i = 1, 4
                    yzel = SS%gl_all_Cell(i, num)
                    r = SS%gl_yzel(:, yzel)
                    kk(i) = 1.0/(max(0.01, sqrt( (x - r(1))**2 + (y - r(2))**2)))
                end do
                kkk = sum(kk)

                yzel = SS%gl_all_Cell(1, num)
                p1 = SS%gl_yzel(:, yzel)

                yzel = SS%gl_all_Cell(2, num)
                p2 = SS%gl_yzel(:, yzel)

                yzel = SS%gl_all_Cell(3, num)
                p3 = SS%gl_yzel(:, yzel)

                yzel = SS%gl_all_Cell(4, num)
                p4 = SS%gl_yzel(:, yzel)

                n1 = SS%gl_Cell_Belong(1:2, 1, num)
                n2 = SS%gl_Cell_Belong(1:2, 3, num)

                n1 = n1/norm2(n1)
                n2 = n2/norm2(n2)

                if( dabs(DOT_PRODUCT(n1, n2)) < 0.999) then  ! Значит можно находить точку пересечения граней
                    a = SS%gl_Cell_Belong(1, 1, num)
                    b = SS%gl_Cell_Belong(2, 1, num)
                    c1 = SS%gl_Cell_Belong(3, 1, num)
                    c = SS%gl_Cell_Belong(1, 3, num)
                    d = SS%gl_Cell_Belong(2, 3, num)
                    c2 = SS%gl_Cell_Belong(3, 3, num)
                    ccc = a * d - b * c
                    x0 = (-c1 * d + c2 * b)/ccc
                    y0 = (c * c1 - a * c2)/ccc

                    a = SS%gl_Cell_Belong(1, 4, num)
                    b = SS%gl_Cell_Belong(2, 4, num)
                    c = SS%gl_Cell_Belong(3, 4, num)
                    t = (-c - a * x0 - b * y0)/(a * (x - x0) + b * (y - y0))
                    d1(1) = x0 + t * (x - x0)
                    d1(2) = y0 + t * (y - y0)

                    a = SS%gl_Cell_Belong(1, 2, num)
                    b = SS%gl_Cell_Belong(2, 2, num)
                    c = SS%gl_Cell_Belong(3, 2, num)
                    t = (-c - a * x0 - b * y0)/(a * (x - x0) + b * (y - y0))
                    d2(1) = x0 + t * (x - x0)
                    d2(2) = y0 + t * (y - y0)
                else
                    a = SS%gl_Cell_Belong(1, 1, num)
                    b = SS%gl_Cell_Belong(2, 1, num)
                    c = -a * x - b * y

                    a2 = SS%gl_Cell_Belong(1, 4, num)
                    b2 = SS%gl_Cell_Belong(2, 4, num)
                    c2 = SS%gl_Cell_Belong(3, 4, num)

                    x0 = (b * c2 - c * b2)/(a * b2 - a2 * b)
                    y0 = (a * c2 - c * a2)/(b * a2 - a * b2)

                    d1(1) = x0
                    d1(2) = y0


                    a2 = SS%gl_Cell_Belong(1, 2, num)
                    b2 = SS%gl_Cell_Belong(2, 2, num)
                    c2 = SS%gl_Cell_Belong(3, 2, num)

                    x0 = (b * c2 - c * b2)/(a * b2 - a2 * b)
                    y0 = (a * c2 - c * a2)/(b * a2 - a * b2)

                    d2(1) = x0
                    d2(2) = y0
                end if

                e1 = norm2(d1 - p4)
                e2 = norm2(d1 - p1)
                e3 = norm2(d2 - p3)
                e4 = norm2(d2 - p2)

                l1 = sqrt((x - d1(1))**2 + (y - d1(2))**2)
                l2 = sqrt((x - d2(1))**2 + (y - d2(2))**2)

                yz1 = SS%gl_all_Cell(1, num)
                yz2 = SS%gl_all_Cell(2, num)
                yz3 = SS%gl_all_Cell(3, num)
                yz4 = SS%gl_all_Cell(4, num)

                if(present(PAR_gd)) then
                    if(size(PAR_gd) /= SS%n_par) STOP "Error Int_Get_Parameter size(PAR_gd) /= SS%n_par   09876tyuinewufhugferafrgnmiutr"
                    PAR_gd = 0.0
                    do i = 1, size(PAR_gd)
                        f1 = (SS%gd(i, yz1) * e1 + SS%gd(i, yz4) * e2)/(e1 + e2)
                        f2 = (SS%gd(i, yz2) * e3 + SS%gd(i, yz3) * e4)/(e3 + e4)
                        PAR_gd(i) = (f1 * l2 + f2 * l1)/(l1 + l2)
                    end do
                end if

                if(present(PAR_atom_source)) then
                    PAR_atom_source = 0.0
                    do i = 1, size(PAR_atom_source)
                        f1 = (SS%atom_source(i, yz1) * e1 + SS%atom_source(i, yz4) * e2)/(e1 + e2)
                        f2 = (SS%atom_source(i, yz2) * e3 + SS%atom_source(i, yz3) * e4)/(e3 + e4)
                        PAR_atom_source(i) = (f1 * l2 + f2 * l1)/(l1 + l2)
                    end do
                end if

                if(present(PAR_hydrogen)) then
                    if(size(PAR_hydrogen(1, :)) /= SS%n_Hidrogen) STOP "Error Int_Get_Parameter size(PAR_hydrogen) /= SS%n_par   oehwjfiehurgegvve"
                    PAR_hydrogen = 0.0
                    do j = 1, size(PAR_hydrogen(1, :))
                        do i = 1, size(PAR_hydrogen(:, 1))
                            f1 = (SS%hydrogen(i, j, yz1) * e1 + SS%hydrogen(i, j, yz4) * e2)/(e1 + e2)
                            f2 = (SS%hydrogen(i, j, yz2) * e3 + SS%hydrogen(i, j, yz3) * e4)/(e3 + e4)
                            PAR_hydrogen(i, j) = (f1 * l2 + f2 * l1)/(l1 + l2)
                        end do
                    end do
                end if
                

                return

                ! end if

                if(present(PAR_gd)) then
                    if(size(PAR_gd) /= SS%n_par) STOP "Error Int_Get_Parameter size(PAR_gd) /= SS%n_par   09876tyuinewufhugferafrgnmiutr"
                    PAR_gd = 0.0
                    do i = 1, 4
                        yzel = SS%gl_all_Cell(i, num)
                        PAR_gd(:) = PAR_gd(:) + kk(i) * SS%gd(:, yzel)
                        !PAR_gd(:) = PAR_gd(:) + SS%gd(:, yzel)
                    end do
                    PAR_gd = PAR_gd/kkk
                    !PAR_gd = PAR_gd/4.0
                end if

                if(present(PAR_atom_source)) then
                    PAR_atom_source = 0.0
                    do i = 1, 4
                        yzel = SS%gl_all_Cell(i, num)
                        PAR_atom_source(:) = PAR_atom_source(:) + kk(i) * SS%atom_source(1:4, yzel)
                        !PAR_atom_source(:) = PAR_atom_source(:) + SS%atom_source(1:4, yzel)
                    end do
                    PAR_atom_source = PAR_atom_source/kkk
                    !PAR_atom_source = PAR_atom_source/4.0
                end if

                if(present(PAR_hydrogen)) then
                    if(size(PAR_hydrogen(1, :)) /= SS%n_Hidrogen) STOP "Error Int_Get_Parameter size(PAR_hydrogen) /= SS%n_par   oehwjfiehurgegvve"
                    PAR_hydrogen = 0.0
                    do i = 1, 4
                        yzel = SS%gl_all_Cell(i, num)
                        PAR_hydrogen = PAR_hydrogen + kk(i) * SS%hydrogen(:, :, yzel)
                        !PAR_hydrogen = PAR_hydrogen + SS%hydrogen(:, :, yzel)
                    end do
                    PAR_hydrogen = PAR_hydrogen/kkk
                    !PAR_hydrogen = PAR_hydrogen/4.0
                end if


                return

                M = SS%gl_Cell_interpol_matrix(:, :, num)

                MM = 0.0

                v(1) = x
                v(2) = y
                v(3) = x * y
                v(4) = 1.0

                do i = 1, 4
                    do j = 1, 4
						ccc = v(i) * M(j, i)
                        MM(i) = MM(i) + v(j) * M(j, i)
                    end do
                end do

                if(present(PAR_gd)) then
                    if(size(PAR_gd) /= SS%n_par) STOP "Error Int_Get_Parameter size(PAR_gd) /= SS%n_par   09876tyuinewufhugferafrgnmiutr"
                    PAR_gd = 0.0
                    do i = 1, 4
                        yzel = SS%gl_all_Cell(i, num)
                        PAR_gd(:) = PAR_gd(:) + MM(i) * SS%gd(:, yzel)
                    end do

                    ! if(PAR_gd(1) <= 0.0) then
                    !     Print*, "ERROR 887 <0 ve5bt4vv3v676543565e74yvbyr"
                    !     pause
                    !     STOP
                    ! end if

                end if

                if(present(PAR_atom_source)) then
                    PAR_atom_source = 0.0
                    do i = 1, 4
                        yzel = SS%gl_all_Cell(i, num)
                        PAR_atom_source(:) = PAR_atom_source(:) + MM(i) * SS%atom_source(1:4, yzel)
                    end do
                end if

                if(present(PAR_hydrogen)) then
                    if(size(PAR_hydrogen(1, :)) /= SS%n_Hidrogen) STOP "Error Int_Get_Parameter size(PAR_hydrogen) /= SS%n_par   oehwjfiehurgegvve"
                    PAR_hydrogen = 0.0
                    do i = 1, 4
                        yzel = SS%gl_all_Cell(i, num)
                        PAR_hydrogen = PAR_hydrogen + MM(i) * SS%hydrogen(:, :, yzel)
                    end do

                    ! if(PAR_hydrogen(1, 4) > 3.0) then
                    !     Print*, "ERROR 887 <0 y54b6unewvbyerev5e"
					! 	print*, x, y
                    !     do i = 1, 4
                    !         yzel = SS%gl_all_Cell(i, num)
					! 		print*, "------"
                    !         print*, SS%hydrogen(1, 4, yzel)
					! 		print*, SS%gl_yzel(:, yzel)
					! 		print*, "------"
					! 	end do
					! 	call Print_matrix_real(M)
                    !     pause
                    !     STOP
                    ! end if

                    ! if(PAR_hydrogen(1, 3) <= 0.0 .or. PAR_hydrogen(2, 3) <= 0.0) then
                    !     Print*, "ERROR 887 <0 veete4vt545vt43t"
                    !     print*, "__________________"
                    !     call Print_matrix_real(SS%hydrogen(:, :, yzel))
                    !     print*, "__________________"
                    !     call Print_matrix_real(PAR_hydrogen)
                    !     pause
                    !     STOP
                    ! end if

                end if

            end if
        end if

    end subroutine Int_Get_Parameter

    subroutine Int_Print_all_point(SS)
        TYPE (Inter_Setka), intent(in out) :: SS
        integer(4) :: N, i

        open(1, file = 'Interpol_all_point.txt')
        write(1,*) "TITLE = 'HP'  VARIABLES = 'X', 'Y'"

        N = size(SS%gl_yzel(1, :))
        do i = 1, N
            write(1,*) SS%gl_yzel(:, i)
        end do

        close(1)

    end subroutine Int_Print_all_point

    subroutine Int_Print_Cell(SS)
        ! Печатаем все ячейки (есть дублирование), каждая ячейка печатается отдельно
        TYPE (Inter_Setka), intent(in) :: SS
        integer :: N1, N2, i, j, N, node, yz

        N = size(SS%gl_Cell_A) + size(SS%gl_Cell_B) + size(SS%gl_Cell_C)
        open(1, file = 'Interpol_Print_Cell.txt')
        write(1,*) "TITLE = 'HP'  VARIABLES = 'X', 'Y'  ZONE T= 'HP', N= ", 4 * N, ", E =  ", 4 * N , ", F=FEPOINT, ET=LINESEG"

        do j = 1, size(SS%gl_Cell_A(1, :))
            do i = 1, size(SS%gl_Cell_A(:, 1))
                node = SS%gl_Cell_A(i, j)
                yz = SS%gl_all_Cell(1, node)
                write(1,*) SS%gl_yzel(:, yz)
                yz = SS%gl_all_Cell(2, node)
                write(1,*) SS%gl_yzel(:, yz)
                yz = SS%gl_all_Cell(3, node)
                write(1,*) SS%gl_yzel(:, yz)
                yz = SS%gl_all_Cell(4, node)
                write(1,*) SS%gl_yzel(:, yz)
            end do
        end do

        do j = 1, size(SS%gl_Cell_B(1, :))
            do i = 1, size(SS%gl_Cell_B(:, 1))
                node = SS%gl_Cell_B(i, j)
                yz = SS%gl_all_Cell(1, node)
                write(1,*) SS%gl_yzel(:, yz)
                yz = SS%gl_all_Cell(2, node)
                write(1,*) SS%gl_yzel(:, yz)
                yz = SS%gl_all_Cell(3, node)
                write(1,*) SS%gl_yzel(:, yz)
                yz = SS%gl_all_Cell(4, node)
                write(1,*) SS%gl_yzel(:, yz)
            end do
        end do

        do j = 1, size(SS%gl_Cell_C(1, :))
            do i = 1, size(SS%gl_Cell_C(:, 1))
                node = SS%gl_Cell_C(i, j)
                yz = SS%gl_all_Cell(1, node)
                write(1,*) SS%gl_yzel(:, yz)
                yz = SS%gl_all_Cell(2, node)
                write(1,*) SS%gl_yzel(:, yz)
                yz = SS%gl_all_Cell(3, node)
                write(1,*) SS%gl_yzel(:, yz)
                yz = SS%gl_all_Cell(4, node)
                write(1,*) SS%gl_yzel(:, yz)
            end do
        end do
        
        do j = 0, N
            write(1,*) 4 * j + 1, 4 * j + 2
            write(1,*) 4 * j + 2, 4 * j + 3
            write(1,*) 4 * j + 3, 4 * j + 4
            write(1,*) 4 * j + 4, 4 * j + 1
        end do

        close(1)

    end subroutine Int_Print_Cell

    subroutine Int_Print_center(SS)
        ! Печатаем все ячейки (есть дублирование), каждая ячейка печатается отдельно
        TYPE (Inter_Setka), intent(in) :: SS
        integer(4) :: N, i

        N = size(SS%gl_Cell_center(1, :)) 
        open(1, file = 'Interpol_Print_Cell_center.txt')
        write(1,*) "TITLE = 'HP'  VARIABLES = 'X', 'Y'"

        do i = 1, N
            write(1,*) SS%gl_Cell_center(:, i)
        end do

        close(1)
        
    end subroutine Int_Print_center

    subroutine Int_Print_connect(SS)
        ! Печатаем все ячейки (есть дублирование), каждая ячейка печатается отдельно
        TYPE (Inter_Setka), intent(in) :: SS
        integer(4) :: N, i, cell, j

        N = size(SS%gl_Cell_center(1, :)) 
        open(1, file = 'Interpol_Print_connect.txt')
        write(1,*) "TITLE = 'HP'  VARIABLES = 'X', 'Y'  ZONE T= 'HP', N= ", 5 * N, ", E =  ", 4 * N , ", F=FEPOINT, ET=LINESEG"

        do i = 1, N
            write(1,*) SS%gl_Cell_center(:, i)
            do j = 1, 4
                cell = SS%gl_Cell_neighbour(j, i)
                if(cell == 0) then
                    write(1,*) SS%gl_Cell_center(:, i)
                else
                    write(1,*) (SS%gl_Cell_center(:, cell) + SS%gl_Cell_center(:, i))/2.0
                end if
            end do
        end do

        do j = 0, N - 1
            write(1,*) 5 * j + 1, 5 * j + 2
            write(1,*) 5 * j + 1, 5 * j + 3
            write(1,*) 5 * j + 1, 5 * j + 4
            write(1,*) 5 * j + 1, 5 * j + 5
        end do

        close(1)
    end subroutine Int_Print_connect

end module Interpol