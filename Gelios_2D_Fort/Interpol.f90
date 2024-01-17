module Interpol
    use STORAGE 
    implicit none 

    contains

    subroutine Int_Init(SS_int, SS)
        TYPE (Setka), intent(in) :: SS
        TYPE (Inter_Setka), intent(in out) :: SS_int

        integer(4) :: N1, N2, i, cell, N3, M1, M2, M3, j, it, nn, yzel, a1, a2, cell2
        real(8) :: p1(2)

        if(SS%init_geo == .False.) STOP "ERROR  Int_Init  11 987r4yg3hjfkgf34f43"
        if(SS_int%init == .True.) STOP "ERROR  Int_Init  12 12edrc4tveyrbunimo98,0p78i6n7ub6y"

        SS_int%init = .True.

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
        allocate(SS_int%gl_Cell_center(2, N1))
        allocate(SS_int%gl_Cell_neighbour(4, N1))
        allocate(SS_int%gl_Cell_Belong(3, 4, N1))
		
		SS_int%gl_all_Cell = 0
        SS_int%gl_Cell_Belong = 0.0
        SS_int%gl_Cell_neighbour = 0
        SS_int%gl_Cell_center = 0.0

        N1 = size(SS%gl_all_Cell(1, :))
        allocate(SS_int%gl_yzel(2, N1 + SS%par_n_END - 1 + SS%par_n_TS + SS%par_m_O - 1 + 1))  ! Добавляем центр координат
		SS_int%gl_yzel = 0.0
		
        ! Заполняем координаты узлов --------------------------------------------------------------------------
		do i = 1, N1
			SS_int%gl_yzel(:, i) = SS%gl_Cell_Centr(:, i, 1)
		end do

        N2 = size(SS%gl_Cell_A(:, 1))

        do i = 1, N2
            cell = SS%gl_Cell_A(i, 1)
            p1 = SS%gl_Cell_Centr(:, cell, 1)
            p1(1) = norm2(p1)
            p1(2) = 0.0  ! Проекция узла на ось симметрии
            SS_int%gl_yzel(:, N1 + i) = p1
		end do
		
		N3 = size(SS%gl_Cell_B(:, 1))
		
		do i = 1, N3
            cell = SS%gl_Cell_B(i, 1)
            p1 = SS%gl_Cell_Centr(:, cell, 1)
            p1(1) = -norm2(p1)
            p1(2) = 0.0  ! Проекция узла на ось симметрии
            SS_int%gl_yzel(:, N1 + N2 + i) = p1
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

    end subroutine Int_Init

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

    subroutine Int_Find_Cell(SS, x, y, num)
        ! Поиск номера ячейки по её координатам
        ! num = предположительный изначальный номер (если не знаем пусть будет равен 1)
        TYPE (Inter_Setka), intent(in) :: SS
        real(8), intent(in) :: x, y
        integer(4), intent(in out) :: num
        integer(4) :: j, gran, sosed, max_num
        LOGICAL :: outer

        max_num = 0
        outer = .False.

        loop1:do while(.TRUE.)
            max_num = max_num + 1
            outer = .False.
            if(max_num > 1000000) then
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