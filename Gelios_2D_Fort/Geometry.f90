module GEOMETRY
    use STORAGE 
    use My_func
    USE ieee_arithmetic
    implicit none 

    contains

    subroutine Init_Setka(SS)        ! Выделение памяти под все массивы сетки
        ! Предполагается, что все параметры сетки определены
        TYPE (Setka), intent(in out) :: SS
        integer(4) :: N1, n2, n, i, cell, j

        if(SS%init_geo == .True.) then
            STOP "ERROR  Init_Setka  987ytuikjhyt"  
        end if
        
        SS%init_geo = .True.

        allocate(SS%gl_RAY_A(SS%par_n_END, SS%par_m_A))
        allocate(SS%gl_RAY_B(SS%par_n_HP, SS%par_m_BC))
        allocate(SS%gl_RAY_C(SS%par_n_END - SS%par_n_HP + 1, SS%par_m_BC))
        allocate(SS%gl_RAY_O(SS%par_n_END - SS%par_n_HP + 1, SS%par_m_O))
        allocate(SS%gl_RAY_K(SS%par_n_TS, SS%par_m_K))
        allocate(SS%gl_RAY_D(SS%par_m_O + 1, SS%par_m_K + 1))
        allocate(SS%gl_RAY_E(SS%par_n_HP - SS%par_n_TS + 1, SS%par_m_O))
        
        allocate(SS%gl_Cell_A(SS%par_n_END - 1, SS%par_m_A + SS%par_m_BC - 1))
        allocate(SS%gl_Cell_B( (SS%par_n_TS - 1) + SS%par_m_O, SS%par_m_K) )
        allocate(SS%gl_Cell_C( SS%par_n_END - SS%par_n_TS, SS%par_m_O))

        N1 = size(SS%gl_Cell_A(:,:)) + size(SS%gl_Cell_B(:,:)) + size(SS%gl_Cell_C(:,:))  ! Число ячеек в сетке
        
        allocate(SS%gl_all_Cell(4, N1) )
        allocate(SS%gl_Cell_neighbour(4, N1))
        allocate(SS%gl_Cell_Centr(2, N1, 2))
        allocate(SS%gl_Cell_type(N1))
        allocate(SS%gl_Cell_number(2, N1))
        allocate(SS%gl_all_Cell_zone(N1))
        SS%gl_all_Cell_zone = 0

        ALLOCATE(SS%gd(SS%n_par, N1 ,2))
        ALLOCATE(SS%hydrogen(5, SS%n_Hidrogen, N1 ,2))
        ALLOCATE(SS%atom_source(SS%n_atom_source, N1))
        ALLOCATE(SS%atom_all_source(4, SS%n_Hidrogen, N1))
        SS%gd = 0.0
        SS%atom_source = 0.0
        SS%atom_all_source = 0.0

        allocate(SS%gl_Cell_gran(4, N1))
        allocate(SS%gl_Cell_gran_dist(4, N1, 2))
        allocate(SS%gl_Cell_belong(3, 4, N1))
        allocate(SS%gl_Cell_square(N1, 2))

        ! Посчитаем число узлов в сетке
        SS%par_n_points = SS%par_n_END * (SS%par_m_A + SS%par_m_BC) + SS%par_m_K * (SS%par_n_TS + SS%par_m_O) + (SS%par_n_END - SS%par_n_TS + 1) * SS%par_m_O - &
                (SS%par_m_A + SS%par_m_BC + SS%par_m_K - 1)  ! Всего точек в сетке
        
        allocate(SS%gl_yzel(2, SS%par_n_points, 2))
        allocate(SS%gl_yzel_Vel(2, SS%par_n_points))
        allocate(SS%gl_Point_num(SS%par_n_points))

        N2 = size(SS%gl_Cell_A(1, :))
        N1 = size(SS%gl_Cell_A(:, 1))
        n = 2.0 * N1 * N2! - N1

        N2 = size(SS%gl_Cell_B(1, :))
        N1 = size(SS%gl_Cell_B(:, 1))
        n = n + 2.0 * N1 * N2 + N1

        N2 = size(SS%gl_Cell_C(1, :))
        N1 = size(SS%gl_Cell_C(:, 1))
        n = n + 2.0 * N1 * N2 + N1

    
        !n1 = 2 * (SS%par_n_END - 1) * (SS%par_m_A + SS%par_m_BC - 1) - (SS%par_n_END - 1) + 2 * (SS%par_n_TS - 1 + SS%par_m_O) * (SS%par_m_K) + &
        !    2 * (SS%par_n_END - SS%par_n_TS) * SS%par_m_O + (SS%par_n_END - SS%par_n_TS)
        !n1 = (n1 + n2)

        allocate(SS%gl_all_Gran(2, n))
        allocate(SS%gl_Gran_neighbour(2, n))
        allocate(SS%gl_Gran_type(n))
        allocate(SS%gl_Gran_shem(n))
        allocate(SS%gl_Gran_normal(2, n, 2))
        allocate(SS%gl_Gran_length(n, 2))
        allocate(SS%gl_Gran_Center(2, n, 2))
        allocate(SS%gl_Gran_neighbour_TVD(2, n))

        allocate(SS%gl_HP( (SS%par_m_O + SS%par_m_A + SS%par_m_BC - 1)  ))   ! Выделяем память под контакт
        allocate(SS%gl_TS( (SS%par_m_A + SS%par_m_BC + SS%par_m_K - 1) ))   ! Выделяем память под TS
        allocate(SS%gl_BS( (SS%par_m_A - 1) ))   ! Выделяем память под BS

        SS%gl_HP = 0
        SS%gl_TS = 0
        SS%gl_BS = 0
        SS%gl_Gran_type = 0
        SS%gl_Gran_shem = 1          !  HLL по умолчанию

        SS%gl_Gran_neighbour_TVD = 0
        SS%gl_all_Gran = 0
        SS%gl_Gran_neighbour = 0
        SS%gl_Cell_gran = 0
        SS%gl_Cell_gran_dist = 0.0
        SS%gl_Cell_belong = 0.0
        SS%gl_Gran_normal = 0.0
        SS%gl_Gran_length = 0.0
        SS%gl_Gran_Center = 0.0
        SS%gl_Cell_square = 0.0

        SS%gl_Cell_type = "-"
        SS%gl_Cell_number = 0

        SS%gl_RAY_A = 0
        SS%gl_RAY_B = 0
        SS%gl_RAY_C = 0
        SS%gl_RAY_O = 0
        SS%gl_RAY_K = 0
        SS%gl_RAY_D = 0
        SS%gl_RAY_E = 0

        SS%gl_Cell_A = 0
        SS%gl_Cell_B = 0
        SS%gl_Cell_C = 0

        SS%gl_all_Cell = 0

        SS%gl_Cell_neighbour = 0

        SS%gl_yzel = 0.0_8
        SS%gl_Cell_Centr = 0.0

        SS%gl_yzel_Vel = 0.0

    end subroutine Init_Setka

    subroutine Dell_Setka(SS)
        TYPE (Setka), intent(in out) :: SS

        deallocate(SS%gl_RAY_A)
        deallocate(SS%gl_RAY_B)
        deallocate(SS%gl_RAY_C)
        deallocate(SS%gl_RAY_O)
        deallocate(SS%gl_RAY_K)
        deallocate(SS%gl_RAY_D)
        deallocate(SS%gl_RAY_E)

        deallocate(SS%gl_Cell_A)
        deallocate(SS%gl_Cell_B)
        deallocate(SS%gl_Cell_C)

        deallocate(SS%gl_all_Cell)
        deallocate(SS%gl_Cell_neighbour)
        deallocate(SS%gl_Cell_Centr)
        deallocate(SS%gl_Cell_type)
        deallocate(SS%gl_Cell_number)

        deallocate(SS%gd)
        deallocate(SS%hydrogen)

        deallocate(SS%gl_Cell_gran)
        deallocate(SS%gl_Cell_gran_dist)
        deallocate(SS%gl_Cell_belong)
        deallocate(SS%gl_Cell_square)

        deallocate(SS%gl_yzel)
        deallocate(SS%gl_yzel_Vel)
        deallocate(SS%gl_Point_num)

        deallocate(SS%gl_all_Gran)
        deallocate(SS%gl_Gran_neighbour)
        deallocate(SS%gl_Gran_type)
        deallocate(SS%gl_Gran_normal)
        deallocate(SS%gl_Gran_length)
        deallocate(SS%gl_Gran_Center)

        deallocate(SS%gl_HP)
        deallocate(SS%gl_TS)
        deallocate(SS%gl_BS)
    end subroutine Dell_Setka

    subroutine Build_Setka_start(SS)        ! Начальное построение сетки
        TYPE (Setka), intent(in out) :: SS

        integer(4) :: i, j, N1, N2, i1, node, kk2, ni, num1, yz1, yz2, gran
        real(8) :: r, phi, the, xx, x, y, z, rr, x2, y2, z2, R_HP, the2

        node = 1
        ! Задаём первый узел
        SS%gl_yzel(:, node, 1) = (/ 0.0_8, 0.0_8 /)
        node = node + 1

        ! **************************************************************** Начнём с узлов сетки

        ! Цикл генерации точек на лучах А и их связывание с этими лучами ************************************************************
        N2 = size(SS%gl_RAY_A(1, :))
        N1 = size(SS%gl_RAY_A(:, 1))
        do j = 1, N2
            do i = 1, N1
                if (i == 1) then
                    SS%gl_RAY_A(i, j) = 1
                    CYCLE
                end if

                SS%gl_RAY_A(i, j) = node
                call Set_Ray_A(SS, i, j, SS%par_R_character, SS%par_R_character * 1.3, SS%par_R_character * 2.0, 1)
                node = node + 1

            end do
        end do

        ! Цикл генерации точек на лучах B и их связывание с этими лучами ************************************************************
        N2 = size(SS%gl_RAY_B(1, :))
        N1 = size(SS%gl_RAY_B(:, 1))
        do j = 1, N2
            do i = 1, N1
                if (i == 1) then
                    SS%gl_RAY_B(i, j) = 1
                    CYCLE
                end if

                the = par_pi/2.0 + (j) * SS%par_triple_point/(N2)
		        the2 = par_pi/2.0 + (j) * SS%par_triple_point_2/(N2)


                R_HP = (SS%par_R_character * 1.3 - SS%par_R_character * sin(the))/sin(the2)

                SS%gl_RAY_B(i, j) = node
                call Set_Ray_B(SS, i, j, SS%par_R_character, R_HP, 1)
                node = node + 1

            end do
        end do

        ! Цикл генерации точек на лучах C и их связывание с этими лучами ************************************************************
        N2 = size(SS%gl_RAY_C(1, :))
        N1 = size(SS%gl_RAY_C(:, 1))
        do j = 1, N2
            do i = 1, N1

                if (i == 1) then
                    SS%gl_RAY_C(i, j) = SS%gl_RAY_B(size(SS%gl_RAY_B(:, 1)), j)
                    CYCLE
                end if

                SS%gl_RAY_C(i, j) = node
                call Set_Ray_C(SS, i, j, 1)
                node = node + 1

            end do
        end do

        ! Цикл генерации точек на лучах O и их связывание с этими лучами ************************************************************
        N2 = size(SS%gl_RAY_O(1, :))
        N1 = size(SS%gl_RAY_O(:, 1))
        do j = 1, N2
            do i = 1, N1

                SS%gl_RAY_O(i, j) = node
                call Set_Ray_O(SS, i, j, SS%par_R_character * 1.3, 1)
                node = node + 1

            end do
        end do

        ! Цикл генерации точек на лучах K и их связывание с этими лучами ************************************************************
        N2 = size(SS%gl_RAY_K(1, :))
        N1 = size(SS%gl_RAY_K(:, 1))
        do j = 1, N2
            do i = 1, N1
                if (i == 1) then
                    SS%gl_RAY_K(i, j) = 1
                    CYCLE
                end if
                SS%gl_RAY_K(i, j) = node
                call Set_Ray_K(SS, i, j, SS%par_R_character, 1)
                node = node + 1

            end do
        end do

        ! Цикл генерации точек на лучах D и их связывание с этими лучами ************************************************************
        N2 = size(SS%gl_RAY_D(1, :))
        N1 = size(SS%gl_RAY_D(:, 1))
        do j = 1, N2
            do i = 1, N1
                if(i == 1) then
					if(j <= size(SS%gl_RAY_K(1, :))) then
                        SS%gl_RAY_D(i, j) = SS%gl_RAY_K(size(SS%gl_RAY_K(:, 1)), j)
					else
						SS%gl_RAY_D(i, j) = SS%gl_RAY_B(SS%par_n_TS, size(SS%gl_RAY_B(1, :)))
					end if
					CYCLE
                end if

                SS%gl_RAY_D(i, j) = node
                call Set_Ray_D(SS, i, j, 1)
                node = node + 1
            end do
        end do

        ! Цикл генерации точек на лучах E и их связывание с этими лучами ************************************************************
        N2 = size(SS%gl_RAY_E(1, :))
        N1 = size(SS%gl_RAY_E(:, 1))
        do j = 1, N2
            do i = 1, N1
                if(i == 1) then
                    SS%gl_RAY_E(i, j) = SS%gl_RAY_D(j + 1, size(SS%gl_RAY_D(1, :)))
                    CYCLE
                end if

                if(i == N1) then
                    SS%gl_RAY_E(i, j) = SS%gl_RAY_O(1, j)
                    CYCLE
                end if

                SS%gl_RAY_E(i, j) = node
                call Set_Ray_E(SS, i, j, 1)
                node = node + 1
            end do
        end do

        SS%gl_yzel(:, :, 2) = SS%gl_yzel(:, :, 1)

        ! Строим сами ячейки и связываем их с точками

        ! А - группа ячеек ************************************************************************************************************************
        N2 = size(SS%gl_Cell_A(1, :))
        N1 = size(SS%gl_Cell_A(:, 1))

        i1 = 1

        do j = 1, N2
            do i = 1, N1

                if (j < SS%par_m_A) then
                    SS%gl_all_Cell(1, i1) = SS%gl_RAY_A(i, j)
                    SS%gl_all_Cell(2, i1) = SS%gl_RAY_A(i + 1, j)
                    SS%gl_all_Cell(3, i1) = SS%gl_RAY_A(i + 1, j + 1)
                    SS%gl_all_Cell(4, i1) = SS%gl_RAY_A(i, j + 1)
                else if (j > SS%par_m_A) then
                    if (i < SS%par_n_HP) then
                        SS%gl_all_Cell(1, i1) = SS%gl_RAY_B(i, j - SS%par_m_A)
                        SS%gl_all_Cell(2, i1) = SS%gl_RAY_B(i + 1, j - SS%par_m_A)
                        SS%gl_all_Cell(3, i1) = SS%gl_RAY_B(i + 1, j - SS%par_m_A + 1)
                        SS%gl_all_Cell(4, i1) = SS%gl_RAY_B(i, j - SS%par_m_A + 1)
                    else if (i == SS%par_n_HP) then
                        SS%gl_all_Cell(1, i1) = SS%gl_RAY_B(i, j - SS%par_m_A)
                        SS%gl_all_Cell(2, i1) = SS%gl_RAY_C(2, j - SS%par_m_A)
                        SS%gl_all_Cell(3, i1) = SS%gl_RAY_C(2, j - SS%par_m_A + 1)
                        SS%gl_all_Cell(4, i1) = SS%gl_RAY_B(i, j - SS%par_m_A + 1)
                    else 
                        SS%gl_all_Cell(1, i1) = SS%gl_RAY_C(i - SS%par_n_HP + 1, j - SS%par_m_A)
                        SS%gl_all_Cell(2, i1) = SS%gl_RAY_C(i - SS%par_n_HP + 2, j - SS%par_m_A)
                        SS%gl_all_Cell(3, i1) = SS%gl_RAY_C(i - SS%par_n_HP + 2, j - SS%par_m_A + 1)
                        SS%gl_all_Cell(4, i1) = SS%gl_RAY_C(i - SS%par_n_HP + 1, j - SS%par_m_A + 1)
                    end if
                else
                    if (i < SS%par_n_HP) then
                        SS%gl_all_Cell(1, i1) = SS%gl_RAY_A(i, j)
                        SS%gl_all_Cell(2, i1) = SS%gl_RAY_A(i + 1, j)
                        SS%gl_all_Cell(3, i1) = SS%gl_RAY_B(i + 1, 1)
                        SS%gl_all_Cell(4, i1) = SS%gl_RAY_B(i, 1)
                    else if (i == SS%par_n_HP) then
                        SS%gl_all_Cell(1, i1) = SS%gl_RAY_A(i, j)
                        SS%gl_all_Cell(2, i1) = SS%gl_RAY_A(i + 1, j)
                        SS%gl_all_Cell(3, i1) = SS%gl_RAY_C(2, 1)
                        SS%gl_all_Cell(4, i1) = SS%gl_RAY_B(i, 1)
                    else ! i > par_n_HP
                        SS%gl_all_Cell(1, i1) = SS%gl_RAY_A(i, j)
                        SS%gl_all_Cell(2, i1) = SS%gl_RAY_A(i + 1, j)
                        SS%gl_all_Cell(3, i1) = SS%gl_RAY_C(i - SS%par_n_HP + 2, 1)
                        SS%gl_all_Cell(4, i1) = SS%gl_RAY_C(i - SS%par_n_HP + 1, 1)
                    end if
                end if

                SS%gl_Cell_A(i, j) = i1
                i1 = i1 + 1

            end do
        end do

        ! B - группа ячеек ***********************************************************************************************************************
        N2 = size(SS%gl_Cell_B(1, :))
        N1 = size(SS%gl_Cell_B(:, 1))

        do j = 1, N2
            do i = 1, N1

                if (j < SS%par_m_K) then
                    if (i < SS%par_n_TS) then
                        SS%gl_all_Cell(1, i1) = SS%gl_RAY_K(i, j)
                        SS%gl_all_Cell(2, i1) = SS%gl_RAY_K(i, j + 1)
                        SS%gl_all_Cell(3, i1) = SS%gl_RAY_K(i + 1, j + 1)
                        SS%gl_all_Cell(4, i1) = SS%gl_RAY_K(i + 1, j)
                    else if (i == SS%par_n_TS) then
                        SS%gl_all_Cell(1, i1) = SS%gl_RAY_K(i, j)
                        SS%gl_all_Cell(2, i1) = SS%gl_RAY_K(i, j + 1)
                        SS%gl_all_Cell(3, i1) = SS%gl_RAY_D(2, j + 1)
                        SS%gl_all_Cell(4, i1) = SS%gl_RAY_D(2, j)
                    else  ! i > par_n_TS
                        SS%gl_all_Cell(1, i1) = SS%gl_RAY_D(i - SS%par_n_TS + 1, j)
                        SS%gl_all_Cell(2, i1) = SS%gl_RAY_D(i - SS%par_n_TS + 1, j + 1)
                        SS%gl_all_Cell(3, i1) = SS%gl_RAY_D(i - SS%par_n_TS + 1 + 1, j + 1)
                        SS%gl_all_Cell(4, i1) = SS%gl_RAY_D(i - SS%par_n_TS + 1 + 1, j)
                    end if
                else  ! j == par_m_K
                    if (i < SS%par_n_TS) then
                        SS%gl_all_Cell(1, i1) = SS%gl_RAY_K(i, j)
                        SS%gl_all_Cell(2, i1) = SS%gl_RAY_B(i, size(SS%gl_RAY_B(i, :)))
                        SS%gl_all_Cell(3, i1) = SS%gl_RAY_B(i + 1, size(SS%gl_RAY_B(i, :)))
                        SS%gl_all_Cell(4, i1) = SS%gl_RAY_K(i + 1, j)
                    else if (i == SS%par_n_TS) then
                        SS%gl_all_Cell(1, i1) = SS%gl_RAY_K(i, j)
                        SS%gl_all_Cell(2, i1) = SS%gl_RAY_B(i, size(SS%gl_RAY_B(i, :)))
                        SS%gl_all_Cell(3, i1) = SS%gl_RAY_D(2, size(SS%gl_RAY_D(2, :)))
                        SS%gl_all_Cell(4, i1) = SS%gl_RAY_D(2, size(SS%gl_RAY_D(2, :)) - 1)
                    else   ! i > par_n_TS
                        SS%gl_all_Cell(1, i1) = SS%gl_RAY_D(i - SS%par_n_TS + 1, j)
                        SS%gl_all_Cell(2, i1) = SS%gl_RAY_D(i - SS%par_n_TS + 1, j + 1)
                        SS%gl_all_Cell(3, i1) = SS%gl_RAY_D(i - SS%par_n_TS + 1 + 1, j + 1)
                        SS%gl_all_Cell(4, i1) = SS%gl_RAY_D(i - SS%par_n_TS + 1 + 1, j)
                    end if
                end if


                SS%gl_Cell_B(i, j) = i1
                i1 = i1 + 1

            end do
        end do

        ! C - группа ячеек ************************************************************************************************************************
        N2 = size(SS%gl_Cell_C(1, :))
        N1 = size(SS%gl_Cell_C(:, 1))

        do j = 1, N2
            do i = 1, N1

                if (j > 1) then
                    if (i < SS%par_n_HP - SS%par_n_TS + 1) then
                        SS%gl_all_Cell(1, i1) = SS%gl_RAY_E(i, j - 1)
                        SS%gl_all_Cell(2, i1) = SS%gl_RAY_E(i + 1, j - 1)
                        SS%gl_all_Cell(3, i1) = SS%gl_RAY_E(i + 1, j - 1 + 1)
                        SS%gl_all_Cell(4, i1) = SS%gl_RAY_E(i, j - 1 + 1)
                    else if (i == SS%par_n_HP - SS%par_n_TS + 1) then
                        SS%gl_all_Cell(1, i1) = SS%gl_RAY_E(i, j - 1)
                        SS%gl_all_Cell(2, i1) = SS%gl_RAY_O(2, j - 1)
                        SS%gl_all_Cell(3, i1) = SS%gl_RAY_O(2, j - 1 + 1)
                        SS%gl_all_Cell(4, i1) = SS%gl_RAY_E(i, j - 1 + 1)
                    else  ! i > par_n_HP - par_n_TS + 1
                        SS%gl_all_Cell(1, i1) = SS%gl_RAY_O(i - (SS%par_n_HP - SS%par_n_TS + 1) + 1, j - 1)
                        SS%gl_all_Cell(2, i1) = SS%gl_RAY_O(i - (SS%par_n_HP - SS%par_n_TS + 1) + 2, j - 1)
                        SS%gl_all_Cell(3, i1) = SS%gl_RAY_O(i - (SS%par_n_HP - SS%par_n_TS + 1) + 2, j - 1 + 1)
                        SS%gl_all_Cell(4, i1) = SS%gl_RAY_O(i - (SS%par_n_HP - SS%par_n_TS + 1) + 1, j - 1 + 1)
                    end if
                else  ! j == 1
                    if (i < SS%par_n_HP - SS%par_n_TS + 1) then
                        SS%gl_all_Cell(1, i1) = SS%gl_RAY_B(SS%par_n_TS + i - 1, size(SS%gl_RAY_B(SS%par_n_TS + i - 1, :)))
                        SS%gl_all_Cell(2, i1) = SS%gl_RAY_B(SS%par_n_TS + i, size(SS%gl_RAY_B(SS%par_n_TS + i - 1, :)))
                        SS%gl_all_Cell(3, i1) = SS%gl_RAY_E(i + 1, j)
                        SS%gl_all_Cell(4, i1) = SS%gl_RAY_E(i, j)
                    else if (i == SS%par_n_HP - SS%par_n_TS + 1) then
                        SS%gl_all_Cell(1, i1) = SS%gl_RAY_B(SS%par_n_TS + i - 1, size(SS%gl_RAY_B(SS%par_n_TS + i - 1, :)))
                        SS%gl_all_Cell(2, i1) = SS%gl_RAY_C(2, size(SS%gl_RAY_C(2, :)))
                        SS%gl_all_Cell(3, i1) = SS%gl_RAY_O(2, j)
                        SS%gl_all_Cell(4, i1) = SS%gl_RAY_E(i, j)
                    else  ! i > par_n_HP - par_n_TS + 1
                        SS%gl_all_Cell(1, i1) = SS%gl_RAY_C(i - (SS%par_n_HP - SS%par_n_TS + 1) + 1, size(SS%gl_RAY_C(1, :)))
                        SS%gl_all_Cell(2, i1) = SS%gl_RAY_C(i - (SS%par_n_HP - SS%par_n_TS + 1) + 2, size(SS%gl_RAY_C(1, :)))
                        SS%gl_all_Cell(3, i1) = SS%gl_RAY_O(i - (SS%par_n_HP - SS%par_n_TS + 1) + 2, j)
                        SS%gl_all_Cell(4, i1) = SS%gl_RAY_O(i - (SS%par_n_HP - SS%par_n_TS + 1) + 1, j)
                    end if
                end if


                SS%gl_Cell_C(i, j) = i1
                i1 = i1 + 1

            end do
        end do

        ! Свяжем ячейки с их соседями (имеется строгий порядок)

        ! Начнём с ячеек группы А
        N2 = size(SS%gl_Cell_A(1, :))
        N1 = size(SS%gl_Cell_A(:, 1))

        do j = 1, N2
            do i = 1, N1

                ! Заполняем первого соседа (по лучу вверх)
                if (i == N1) then
                    if (j < SS%par_m_A) then
                        SS%gl_Cell_neighbour(1, SS%gl_Cell_A(i, j)) = -1   !! Граница (набегающий поток)
                    else
                        SS%gl_Cell_neighbour(1, SS%gl_Cell_A(i, j)) = -3   !! Граница верхний цилиндр
                    end if
                else
                    SS%gl_Cell_neighbour(1, SS%gl_Cell_A(i, j)) = SS%gl_Cell_A(i + 1, j)
                end if

                ! Заполняем второго соседа (по лучу вниз)
                if (i == 1) then
                    continue
                else
                    SS%gl_Cell_neighbour(2, SS%gl_Cell_A(i, j)) = SS%gl_Cell_A(i - 1, j)
                end if

                ! Заполняем третьего соседа (по углу в плоскости, в сторону увеличения)
                if (j == N2) then
                    if (i < SS%par_n_TS) then
                        SS%gl_Cell_neighbour(3, SS%gl_Cell_A(i, j)) = SS%gl_Cell_B(i, SS%par_m_K)
                    else
                        SS%gl_Cell_neighbour(3, SS%gl_Cell_A(i, j)) = SS%gl_Cell_C(i - SS%par_n_TS + 1, 1)
                    end if
                else
                    SS%gl_Cell_neighbour(3, SS%gl_Cell_A(i, j)) = SS%gl_Cell_A(i, j + 1)
                end if

                ! Заполняем четвёртого соседа (по углу в плоскости, в сторону уменьшения)
                if (j == 1) then
                    SS%gl_Cell_neighbour(4, SS%gl_Cell_A(i, j)) = -4  !! Ось симметрии
                    !TODO continue
                else
                    SS%gl_Cell_neighbour(4, SS%gl_Cell_A(i, j)) = SS%gl_Cell_A(i, j - 1)
                end if
            end do
        end do

        ! Ячееки группы B
        N2 = size(SS%gl_Cell_B(1, :))
        N1 = size(SS%gl_Cell_B(:, 1))

        do j = 1, N2
            do i = 1, N1

                ! Заполняем первого соседа (по лучу вверх)
                if (i == N1) then
                    SS%gl_Cell_neighbour(1, SS%gl_Cell_B(i, j)) = -2   ! Выходная граница
                else
                    SS%gl_Cell_neighbour(1, SS%gl_Cell_B(i, j)) = SS%gl_Cell_B(i + 1, j)
                end if

                ! Заполняем второго соседа (по лучу вниз)
                if (i == 1) then
                    continue
                else
                    SS%gl_Cell_neighbour(2, SS%gl_Cell_B(i, j)) = SS%gl_Cell_B(i - 1, j)
                end if

                ! Заполняем третьего соседа (по углу в плоскости наверх на схеме)
                if (j == N2) then
                    if(i < SS%par_n_TS) then
                        SS%gl_Cell_neighbour(3, SS%gl_Cell_B(i, j)) = SS%gl_Cell_A(i, size(SS%gl_Cell_A(i,:)))
                    else
                        SS%gl_Cell_neighbour(3, SS%gl_Cell_B(i, j)) = SS%gl_Cell_C(1, i - SS%par_n_TS + 1)
                    end if
                else
                    SS%gl_Cell_neighbour(3, SS%gl_Cell_B(i, j)) = SS%gl_Cell_B(i, j + 1)
                end if

                ! Заполняем четвёртого соседа (по углу в плоскости)
                if (j == 1) then
                    SS%gl_Cell_neighbour(4, SS%gl_Cell_B(i, j)) = -4
                else
                    SS%gl_Cell_neighbour(4, SS%gl_Cell_B(i, j)) = SS%gl_Cell_B(i, j - 1)
                end if

            end do
        end do

        ! Ячееки группы C
        N2 = size(SS%gl_Cell_C(1, :))
        N1 = size(SS%gl_Cell_C(:, 1))

        do j = 1, N2
            do i = 1, N1

                ! Заполняем первого соседа (по лучу вверх)
                if (i == N1) then
                    SS%gl_Cell_neighbour(1, SS%gl_Cell_C(i, j)) = -3   ! Верхний цилиндр
                else
                    SS%gl_Cell_neighbour(1, SS%gl_Cell_C(i, j)) = SS%gl_Cell_C(i + 1, j)
                end if

                ! Заполняем второго соседа (по лучу вниз)
                if (i == 1) then
                    SS%gl_Cell_neighbour(2, SS%gl_Cell_C(i, j)) = SS%gl_Cell_B(j + SS%par_n_TS - 1, size(SS%gl_Cell_B(j + SS%par_n_TS - 1, :)))
                else
                    SS%gl_Cell_neighbour(2, SS%gl_Cell_C(i, j)) = SS%gl_Cell_C(i - 1, j)
                end if

                ! Заполняем третьего соседа (по углу в плоскости, в сторону увеличения угла!)
                if (j == N2) then
                    SS%gl_Cell_neighbour(3, SS%gl_Cell_C(i, j)) = -2  ! Выходная граница
                else
                    SS%gl_Cell_neighbour(3, SS%gl_Cell_C(i, j)) = SS%gl_Cell_C(i, j + 1)
                end if

                ! Заполняем четвёртого соседа (по углу в плоскости, в сторону уменьшения)
                if (j == 1) then
                    SS%gl_Cell_neighbour(4, SS%gl_Cell_C(i, j)) = SS%gl_Cell_A(i + SS%par_n_TS - 1, size(SS%gl_Cell_A(i + SS%par_n_TS - 1, :)))
                else
                    SS%gl_Cell_neighbour(4, SS%gl_Cell_C(i, j)) = SS%gl_Cell_C(i, j - 1)
                end if

            end do
        end do

        ! Создадим грани

        node = 1

        ! Начнём с ячеек группы А
        N2 = size(SS%gl_Cell_A(1, :))
        N1 = size(SS%gl_Cell_A(:, 1))

        do j = 1, N2
            do i = 1, N1
                ni = SS%gl_Cell_A(i, j)

                ! Заполняем 1-ую грань для этой ячейки
                if (i < N1) then
                    SS%gl_all_Gran(1, node) = SS%gl_all_Cell(2, ni)
                    SS%gl_all_Gran(2, node) = SS%gl_all_Cell(3, ni)

                    SS%gl_Gran_neighbour(1, node) = ni
                    SS%gl_Gran_neighbour(2, node) = SS%gl_Cell_A(i + 1, j)

                    SS%gl_Cell_gran(1, ni) = node
                    SS%gl_Cell_gran(2, SS%gl_Cell_A(i + 1, j)) = node

                    node = node + 1
                else
                    SS%gl_all_Gran(1, node) = SS%gl_all_Cell(2, ni)
                    SS%gl_all_Gran(2, node) = SS%gl_all_Cell(3, ni)

                    SS%gl_Gran_neighbour(1, node) = ni

                    if (j < SS%par_m_A) then
                        SS%gl_Gran_neighbour(2, node) = -1
                    else
                        SS%gl_Gran_neighbour(2, node) = -3
                    end if

                    SS%gl_Cell_gran(1, ni) = node
                    node = node + 1
                end if

                ! Заполняем 4-ую грань для этой ячейки
                if (j /= 1) then
                    SS%gl_all_Gran(1, node) = SS%gl_all_Cell(1, ni)
                    SS%gl_all_Gran(2, node) = SS%gl_all_Cell(2, ni)

                    SS%gl_Gran_neighbour(1, node) = ni
                    SS%gl_Gran_neighbour(2, node) = SS%gl_Cell_A(i, j - 1)

                    SS%gl_Cell_gran(4, ni) = node
                    SS%gl_Cell_gran(3, SS%gl_Cell_A(i, j - 1)) = node

                    node = node + 1
                else !TODO новое (отличие от сетки 3Д)
                    SS%gl_all_Gran(1, node) = SS%gl_all_Cell(1, ni)
                    SS%gl_all_Gran(2, node) = SS%gl_all_Cell(2, ni)

                    SS%gl_Gran_neighbour(1, node) = ni
                    SS%gl_Gran_neighbour(2, node) = -4

                    SS%gl_Cell_gran(4, ni) = node

                    node = node + 1
                end if
                
            end do
        end do

        ! Ячееки группы B
        N2 = size(SS%gl_Cell_B(1, :))
        N1 = size(SS%gl_Cell_B(:, 1))

        do j = 1, N2
            do i = 1, N1

                ni = SS%gl_Cell_B(i, j)

                ! Заполняем 1-ую грань для этой ячейки
                if (i < N1) then
                    SS%gl_all_Gran(1, node) = SS%gl_all_Cell(3, ni)
                    SS%gl_all_Gran(2, node) = SS%gl_all_Cell(4, ni)

                    SS%gl_Gran_neighbour(1, node) = ni
                    SS%gl_Gran_neighbour(2, node) = SS%gl_Cell_B(i + 1, j)

                    SS%gl_Cell_gran(1, ni) = node
                    SS%gl_Cell_gran(2, SS%gl_Cell_B(i + 1, j)) = node

                    node = node + 1
                else
                    SS%gl_all_Gran(1, node) = SS%gl_all_Cell(3, ni)
                    SS%gl_all_Gran(2, node) = SS%gl_all_Cell(4, ni)

                    SS%gl_Gran_neighbour(1, node) = ni
                    SS%gl_Gran_neighbour(2, node) = -2

                    SS%gl_Cell_gran(1, ni) = node

                    node = node + 1
                end if

                ! Заполняем 4-ую грань для этой ячейки
                if (j /= 1) then
                    SS%gl_all_Gran(1, node) = SS%gl_all_Cell(4, ni)
                    SS%gl_all_Gran(2, node) = SS%gl_all_Cell(1, ni)

                    SS%gl_Gran_neighbour(1, node) = ni
                    SS%gl_Gran_neighbour(2, node) = SS%gl_Cell_B(i, j - 1)

                    SS%gl_Cell_gran(4, ni) = node
                    SS%gl_Cell_gran(3, SS%gl_Cell_B(i, j - 1)) = node

                    node = node + 1
                else !TODO добавил, не было в 3Д сетке
                    SS%gl_all_Gran(1, node) = SS%gl_all_Cell(4, ni)
                    SS%gl_all_Gran(2, node) = SS%gl_all_Cell(1, ni)

                    SS%gl_Gran_neighbour(1, node) = ni
                    SS%gl_Gran_neighbour(2, node) = -4

                    SS%gl_Cell_gran(4, ni) = node

                    node = node + 1
                end if

                if (j == N2) then  ! 3-яя грань для границы
                    SS%gl_all_Gran(1, node) = SS%gl_all_Cell(2, ni)
                    SS%gl_all_Gran(2, node) = SS%gl_all_Cell(3, ni)

                    SS%gl_Gran_neighbour(1, node) = ni

                    if(i < SS%par_n_TS) then
                        SS%gl_Gran_neighbour(2, node) = SS%gl_Cell_A(i, size(SS%gl_Cell_A(i,:)))
                        SS%gl_Cell_gran(3, SS%gl_Cell_A(i, size(SS%gl_Cell_A(i,:)))) = node
                    else
                        SS%gl_Gran_neighbour(2, node) = SS%gl_Cell_C(1, i - SS%par_n_TS + 1)
                        SS%gl_Cell_gran(2, SS%gl_Cell_C(1, i - SS%par_n_TS + 1)) = node
                    end if

                    SS%gl_Cell_gran(3, ni) = node

                    node = node + 1
                end if
            end do
        end do

        ! Ячейки группы C
        N2 = size(SS%gl_Cell_C(1, :))
        N1 = size(SS%gl_Cell_C(:, 1))

        do j = 1, N2
            do i = 1, N1

                ni = SS%gl_Cell_C(i, j)

                ! Заполняем 1-ую грань для этой ячейки
                if (i < N1) then
                    SS%gl_all_Gran(1, node) = SS%gl_all_Cell(2, ni)
                    SS%gl_all_Gran(2, node) = SS%gl_all_Cell(3, ni)

                    SS%gl_Gran_neighbour(1, node) = ni
                    SS%gl_Gran_neighbour(2, node) = SS%gl_Cell_C(i + 1, j)

                    SS%gl_Cell_gran(1, ni) = node
                    SS%gl_Cell_gran(2, SS%gl_Cell_C(i + 1, j)) = node

                    node = node + 1
                else
                    SS%gl_all_Gran(1, node) = SS%gl_all_Cell(2, ni)
                    SS%gl_all_Gran(2, node) = SS%gl_all_Cell(3, ni)

                    SS%gl_Gran_neighbour(1, node) = ni
                    SS%gl_Gran_neighbour(2, node) = -3    ! Верхний цилиндр

                    SS%gl_Cell_gran(1, ni) = node

                    node = node + 1
                end if

                ! Заполняем 4-ую грань для этой ячейки
                if (j /= 1) then
                    SS%gl_all_Gran(1, node) = SS%gl_all_Cell(1, ni)
                    SS%gl_all_Gran(2, node) = SS%gl_all_Cell(2, ni)

                    SS%gl_Gran_neighbour(1, node) = ni
                    SS%gl_Gran_neighbour(2, node) = SS%gl_Cell_C(i, j - 1)

                    SS%gl_Cell_gran(4, ni) = node
                    SS%gl_Cell_gran(3, SS%gl_Cell_C(i, j - 1)) = node

                    node = node + 1
                end if

                if (j == N2) then  ! Третья грань на границе
                    SS%gl_all_Gran(1, node) = SS%gl_all_Cell(3, ni)
                    SS%gl_all_Gran(2, node) = SS%gl_all_Cell(4, ni)

                    SS%gl_Gran_neighbour(1, node) = ni
                    SS%gl_Gran_neighbour(2, node) = -2

                    SS%gl_Cell_gran(3, ni) = node
                    node = node + 1
                end if


                if (j == 1) then  ! 4-яя грань для границы
                    SS%gl_all_Gran(1, node) = SS%gl_all_Cell(1, ni)
                    SS%gl_all_Gran(2, node) = SS%gl_all_Cell(2, ni)

                    SS%gl_Gran_neighbour(1, node) = ni

                    SS%gl_Gran_neighbour(2, node) = SS%gl_Cell_A(i + SS%par_n_TS - 1, size(SS%gl_Cell_A(i + SS%par_n_TS - 1, :)))
                    SS%gl_Cell_gran(3, SS%gl_Cell_A(i + SS%par_n_TS - 1, size(SS%gl_Cell_A(i + SS%par_n_TS - 1, :)))) = node

                    SS%gl_Cell_gran(4, ni) = node

                    node = node + 1
                end if

            end do
        end do

        ! Теперь для каждой ячейки добавим ей информацию какого она типа и какой у неё номер (тройка) в этом типе
        ! Начнём с группы А
        N2 = size(SS%gl_Cell_A(1, :))
        N1 = size(SS%gl_Cell_A(:, 1))

        do j = 1, N2
            do i = 1, N1
                num1 = SS%gl_Cell_A(i, j)
                SS%gl_Cell_type(num1) = "A"
                SS%gl_Cell_number(:, num1) = (/i, j/)
            end do
        end do
        
        ! Начнём с группы B
        N2 = size(SS%gl_Cell_B(1, :))
        N1 = size(SS%gl_Cell_B(:, 1))

        do j = 1, N2
            do i = 1, N1
                num1 = SS%gl_Cell_B(i, j)
                SS%gl_Cell_type(num1) = "B"
                SS%gl_Cell_number(:, num1) = (/i, j/)
            end do
        end do
        
        ! Начнём с группы C
        N2 = size(SS%gl_Cell_C(1, :))
        N1 = size(SS%gl_Cell_C(:, 1))

        do j = 1, N2
            do i = 1, N1
                num1 = SS%gl_Cell_C(i, j)
                SS%gl_Cell_type(num1) = "C"
                SS%gl_Cell_number(:, num1) = (/i, j/)
            end do
        end do

        call Geo_Find_Surface(SS)  ! Находим поверхности, которые выделяем
        call Geo_Culc_normal(SS, 1)
        call Geo_Culc_normal(SS, 2)
        call Culc_Cell_Centr(SS, 1)
        call Culc_Cell_Centr(SS, 2)

        call Belong_Init(SS)

        call Geo_Culc_length_area(SS, 1)
        call Geo_Culc_length_area(SS, 2)

        call Geo_culc_TVD_sosed(SS)
        call Geo_Culc_zone(SS)


        
    end subroutine Build_Setka_start

    subroutine Belong_Init(SS)
        ! Для каждой ячейки найдём коэффиценты каждой грани для определения принадлежности точки к ячейке
        TYPE (Setka), intent(in out) :: SS
        integer(4) :: N1, i, j, gran, yz1, yz2
        real(8) :: c, p(2), n(2), centr(2)

        N1 = size(SS%gl_Cell_gran(1, :))

        do i = 1, N1  ! По всем ячейкам
            do j = 1, 4  ! По всем граням
                gran = SS%gl_Cell_gran(j, i)

                if (gran == 0) then
                    SS%gl_Cell_belong(:, j, i) = 0.0
                    CYCLE
                end if

                yz1 = SS%gl_all_Gran(1, gran)
                yz2 = SS%gl_all_Gran(2, gran)
                if(yz1 == yz2) then
                    SS%gl_Cell_belong(:, j, i) = 0.0
                    CYCLE
                end if
                p = SS%gl_yzel(:, yz1, 1)
                n = SS%gl_Gran_normal(:, gran, 1)
                c = -DOT_PRODUCT(p, n)

                SS%gl_Cell_belong(1, j, i) = n(1)
                SS%gl_Cell_belong(2, j, i) = n(2)
                SS%gl_Cell_belong(3, j, i) = c
                centr = SS%gl_Cell_Centr(:, i, 1)

                if(centr(1) * n(1) + centr(2) * n(2) + c > 0) then
                    SS%gl_Cell_belong(:, j, i) = -SS%gl_Cell_belong(:, j, i)
                end if
            end do
        end do

    end subroutine Belong_Init

    subroutine Set_Ray_A(SS, i, j, R_TS, R_HP, R_BS, step)
        ! step - 1 или 2  показывает какой массив координат мы меняем
        TYPE (Setka), intent(in out) :: SS
        integer(4), INTENT(IN) :: i, j, step
        real(8), INTENT(IN) :: R_TS, R_HP, R_BS

        integer(4) :: node, N2
        real(8) :: the, r, rr, dd

        node = SS%gl_RAY_A(i, j)
        N2 = size(SS%gl_RAY_A(1, :))

        if (i == 1) then
            return
        end if

        ! Вычисляем координаты текущего луча в пространстве
        the = (j - 1) * par_pi/2.0/(N2 - 1)
        ! Вычисляем координаты точки на луче

        if (i <= SS%par_n_IB) then  ! NEW
            r =  2.0 * SS%par_R0 + (SS%par_R_inner - 2.0 * SS%par_R0) * (DBLE(i - 2)/(SS%par_n_IB - 2))**SS%par_kk1
		else if (i <= SS%par_n_TS) then  
			r =  SS%par_R_inner + (R_TS - SS%par_R_inner) * sgushenie_3( (DBLE(i - SS%par_n_IB)/(SS%par_n_TS - SS%par_n_IB)) , SS%par_kk12)
		else if (i <= SS%par_n_HP) then  
			r = R_TS + (R_HP - R_TS) * sgushenie_2(DBLE(i - SS%par_n_TS)/(SS%par_n_HP - SS%par_n_TS), SS%par_kk14)
		else if (i <= SS%par_n_BS) then 
			r = R_HP + (R_BS - R_HP) * (DBLE(i - SS%par_n_HP)/(SS%par_n_BS - SS%par_n_HP))**1.6
		else if (i <= SS%par_n_BS + 4) then 
            dd = (R_BS - R_HP) * ( 1.0 - (DBLE(SS%par_n_BS - 1 - SS%par_n_HP)/(SS%par_n_BS - SS%par_n_HP))**1.6)
            r = R_BS + (i - SS%par_n_BS) * dd
		else
            dd = (R_BS - R_HP) * ( 1.0 - (DBLE(SS%par_n_BS - 1 - SS%par_n_HP)/(SS%par_n_BS - SS%par_n_HP))**1.6)
            rr = R_BS + (4) * dd
			r = rr + (SS%par_R_END - rr) * (DBLE(i - SS%par_n_BS - 4)/(SS%par_n_END - SS%par_n_BS - 4))**(SS%par_kk2 * (0.55 + 0.45 * cos(the)) )
			!r = R_BS + (SS%par_R_END - R_BS) * (DBLE(i- SS%par_n_BS)/(SS%par_n_END - SS%par_n_BS))**(SS%par_kk2 * (0.55 + 0.45 * cos(the)) )
		end if

        SS%gl_yzel(1, node, step) = r * cos(the)
        SS%gl_yzel(2, node, step) = r * sin(the)
    end subroutine Set_Ray_A

    subroutine Set_Ray_B(SS, i, j, R_TS, R_HP, step)
        ! step - 1 или 2  показывает какой массив координат мы меняем
        ! R_TS - это расстояние от TS до центра
        !! R_HP - это расстояние по второй координате от TS до HP
        TYPE (Setka), intent(in out) :: SS
        integer(4), INTENT(IN) :: i, j, step
        real(8), INTENT(IN) :: R_TS, R_HP

        integer(4) :: node, N2
        real(8) :: the, r, the2

        node = SS%gl_RAY_B(i, j)
        N2 = size(SS%gl_RAY_B(1, :))

        if (i == 1) then
            return
        end if

        ! Вычисляем координаты текущего луча в пространстве
        the = par_pi/2.0 + (j) * SS%par_triple_point/(N2)
		the2 = par_pi/2.0 + (j) * SS%par_triple_point_2/(N2)
        ! Вычисляем координаты точки на луче

        ! до TS
        if (i <= SS%par_n_IB) then  
            r =  2.0 * SS%par_R0 + (SS%par_R_inner - 2.0 * SS%par_R0) * (DBLE(i - 2)/(SS%par_n_IB - 2))**SS%par_kk1
        else if (i <= SS%par_n_TS) then  
            r =  SS%par_R_inner + (R_TS - SS%par_R_inner) * sgushenie_3( (DBLE(i - SS%par_n_IB)/(SS%par_n_TS - SS%par_n_IB)) , SS%par_kk12)
        else if (i <= SS%par_n_HP) then  
            r =  R_HP * sgushenie_2(DBLE(i - SS%par_n_TS)/(SS%par_n_HP - SS%par_n_TS), SS%par_kk14)

            SS%gl_yzel(1, node, step) = R_TS * cos(the) + r * cos(the2)
            SS%gl_yzel(2, node, step) = R_TS * sin(the) + r * sin(the2)

            return
        end if

        SS%gl_yzel(1, node, step) = r * cos(the)
        SS%gl_yzel(2, node, step) = r * sin(the)

        return
    end subroutine Set_Ray_B

    subroutine Set_Ray_C(SS, i, j, step)
        ! step - 1 или 2  показывает какой массив координат мы меняем
        ! Этим лучам ничего не надо, они сами подстраиваются
        TYPE (Setka), intent(in out) :: SS
        integer(4), INTENT(IN) :: i, j, step

        integer(4) :: node, N2, N1
        real(8) :: the, r, rr, xx, R_BS, dd, rrr

        node = SS%gl_RAY_C(i, j)
        N2 = size(SS%gl_RAY_C(1, :))
        N1 = size(SS%gl_RAY_C(:, 1))

        if (i == 1) then
            return
        end if

        ! Вычисляем координаты точки на луче
        xx = SS%gl_yzel(1, SS%gl_RAY_B(SS%par_n_HP, j), step)
        rr = SS%gl_yzel(2, SS%gl_RAY_B(SS%par_n_HP, j), step)
        R_BS = SS%gl_yzel(2, SS%gl_RAY_A(SS%par_n_BS, SS%par_m_A), step)


        if (i <= SS%par_n_BS - SS%par_n_HP + 1) then
			r = rr + (R_BS - rr) * (DBLE(i)/(SS%par_n_BS - SS%par_n_HP + 1))**1.6
        else if (i <= SS%par_n_BS - SS%par_n_HP + 5) then
            dd = 1.07 * (R_BS - rr) * ( 1.0 - (DBLE(SS%par_n_BS - SS%par_n_HP)/(SS%par_n_BS - SS%par_n_HP + 1))**1.6 )
			r = R_BS + dd * (i - (SS%par_n_BS - SS%par_n_HP + 1))
		else
            dd = 1.07 * (R_BS - rr) * ( 1.0 - (DBLE(SS%par_n_BS - SS%par_n_HP)/(SS%par_n_BS - SS%par_n_HP + 1))**1.6 )
            rrr = R_BS + dd * 4
			r = rrr + (DBLE(i - (SS%par_n_BS - SS%par_n_HP + 5))/(N1 - (SS%par_n_BS - SS%par_n_HP + 5) ))**(0.55 * SS%par_kk2) * (SS%par_R_END - rrr)
		end if

        SS%gl_yzel(1, node, step) = xx
        SS%gl_yzel(2, node, step) = r

        return
    end subroutine Set_Ray_C

    subroutine Set_Ray_O(SS, i, j, R_HP, step)
        ! step - 1 или 2  показывает какой массив координат мы меняем
        ! R_HP - высота (y - координата)
        ! Распределение x - координат считается автоматически от x - координаты крайней точки B на гелиопаузе
        TYPE (Setka), intent(in out) :: SS
        integer(4), INTENT(IN) :: i, j, step
        real(8), INTENT(IN) :: R_HP

        integer(4) :: node, N2, N1
        real(8) :: the, r, rr, xx, R_BS, x, rrr, dd

        node = SS%gl_RAY_O(i, j)
        N2 = size(SS%gl_RAY_O(1, :))
        N1 = size(SS%gl_RAY_O(:, 1))

        R_BS = SS%gl_yzel(2, SS%gl_RAY_A(SS%par_n_BS, SS%par_m_A), step)
        xx = SS%gl_yzel(1, SS%gl_RAY_B(SS%par_n_HP, SS%par_m_BC), step)              ! Отталкиваемся от x - координаты крайней точки B на гелиопаузе в этой плоскости (k)
        

        x = xx - (DBLE(j)/N2)**SS%par_kk31 * (xx - SS%par_R_LEFT)

        if (i <= SS%par_n_BS - SS%par_n_HP + 1) then
			r = R_HP + (R_BS - R_HP) * (DBLE(i - 1)/(SS%par_n_BS - SS%par_n_HP))**1.6
		else if (i <= SS%par_n_BS - SS%par_n_HP + 5) then
            dd = 1.07 * (R_BS - R_HP) * (1.0 - (DBLE(SS%par_n_BS - SS%par_n_HP - 1)/(SS%par_n_BS - SS%par_n_HP))**1.6)
			!r = R_HP + (R_BS - R_HP) * (DBLE(i - 1)/(SS%par_n_BS - SS%par_n_HP))**1.6
            r = R_BS + dd * (i - (SS%par_n_BS - SS%par_n_HP + 1))
		else
            dd = 1.07 * (R_BS - R_HP) * (1.0 - (DBLE(SS%par_n_BS - SS%par_n_HP - 1)/(SS%par_n_BS - SS%par_n_HP))**1.6)
			rrr = R_BS + 4 * dd
            r = rrr + (DBLE(i - (SS%par_n_BS - SS%par_n_HP + 5))/(N1 - (SS%par_n_BS - SS%par_n_HP + 5) ))**(0.55 * SS%par_kk2) * (SS%par_R_END - rrr)
            !r = R_BS + (DBLE(i - (SS%par_n_BS - SS%par_n_HP + 1))/(N1 - (SS%par_n_BS - SS%par_n_HP + 1) ))**(0.55 * SS%par_kk2) * (SS%par_R_END - R_BS)
		end if

        SS%gl_yzel(1, node, step) = x
        SS%gl_yzel(2, node, step) = r

        return
    end subroutine Set_Ray_O

    subroutine Set_Ray_K(SS, i, j, R_TS, step)
        ! step - 1 или 2  показывает какой массив координат мы меняем
        TYPE (Setka), intent(in out) :: SS
        integer(4), INTENT(IN) :: i, j, step
        real(8), INTENT(IN) :: R_TS

        integer(4) :: node, N2, N1
        real(8) :: the, r, rr, xx, R_BS, x

        if (i == 1) return

        node = SS%gl_RAY_K(i, j)
        N2 = size(SS%gl_RAY_K(1, :))
        N1 = size(SS%gl_RAY_K(:, 1))

        the = par_pi/2.0 + SS%par_triple_point + (N2 - j + 1) * (par_pi/2.0 - SS%par_triple_point)/(N2)

        if (i <= SS%par_n_IB) then  ! NEW
            r =  2.0 * SS%par_R0 + (SS%par_R_inner - 2.0 * SS%par_R0) * (DBLE(i - 2)/(SS%par_n_IB - 2))**SS%par_kk1
        else 
            r =  SS%par_R_inner + (R_TS - SS%par_R_inner) * sgushenie_3(DBLE(i - SS%par_n_IB)/(SS%par_n_TS - SS%par_n_IB), SS%par_kk12)
		end if

        SS%gl_yzel(1, node, step) = r * cos(the)
        SS%gl_yzel(2, node, step) = r * sin(the)

        return
    end subroutine Set_Ray_K

    subroutine Set_Ray_D(SS, i, j, step)
        ! step - 1 или 2  показывает какой массив координат мы меняем
        TYPE (Setka), intent(in out) :: SS
        integer(4), INTENT(IN) :: i, j, step

        integer(4) :: node, N2, N1
        real(8) :: the, r, rr, xx, R_BS, x, y

        if (i == 1) return

        node = SS%gl_RAY_D(i, j)
        N2 = size(SS%gl_RAY_D(1, :))
        N1 = size(SS%gl_RAY_D(:, 1))

        if (j < N2) then
            xx = SS%gl_yzel(1, SS%gl_RAY_K(SS%par_n_TS, j), step)
            y = SS%gl_yzel(2, SS%gl_RAY_K(SS%par_n_TS, j), step)
        else
            xx = SS%gl_yzel(1, SS%gl_RAY_B(SS%par_n_TS, SS%par_m_BC), step)
            y = SS%gl_yzel(2, SS%gl_RAY_B(SS%par_n_TS, SS%par_m_BC), step)
        end if


        x = xx + (DBLE(i - 1)/(N1 - 1))**SS%par_kk3 * (SS%par_R_LEFT - xx)

        SS%gl_yzel(1, node, step) = x
        SS%gl_yzel(2, node, step) = y

        return
    end subroutine Set_Ray_D

    subroutine Set_Ray_E(SS, i, j, step)
        ! step - 1 или 2  показывает какой массив координат мы меняем
        TYPE (Setka), intent(in out) :: SS
        integer(4), INTENT(IN) :: i, j, step

        integer(4) :: node, N2, N1
        real(8) :: the, r, rr, xx, R_BS, x, y, yy

        if (i == 1) return

        node = SS%gl_RAY_E(i, j)
        N2 = size(SS%gl_RAY_E(1, :))
        N1 = size(SS%gl_RAY_E(:, 1))

        if (i == N1) return

        x = SS%gl_yzel(1, SS%gl_RAY_E(1, j), step)
        y = SS%gl_yzel(2, SS%gl_RAY_E(1, j), step)

        xx = SS%gl_yzel(1, SS%gl_RAY_O(1, j), step)
        yy = SS%gl_yzel(2, SS%gl_RAY_O(1, j), step)

        SS%gl_yzel(1, node, step) = x + (xx - x) * sgushenie_2(DBLE(i - 1)/(N1 - 1), SS%par_kk14)
        SS%gl_yzel(2, node, step) = y + (yy - y) * sgushenie_2(DBLE(i - 1)/(N1 - 1), SS%par_kk14)

        return
    end subroutine Set_Ray_E

    subroutine Geo_Find_Cell(SS, x, y, num, inzone)
        ! Поиск номера ячейки по её координатам
        ! num = предположительный изначальный номер (если не знаем пусть будет равен 1)
        TYPE (Setka), intent(in) :: SS
        real(8), intent(in) :: x, y
        integer(4), intent(in out) :: num
        integer(4) :: j, gran, sosed, max_num
        logical, intent(out), OPTIONAL :: inzone

        max_num = 0

        if(num < 1) then
            print*, "ERROR 987y843gobiufhwiuhluhfygoiwhf"
            pause
        end if


        loop1:do while(.TRUE.)
            max_num = max_num + 1
			if(present(inzone)) inzone = .True.
			
            if(max_num > 1000000) then
                STOP "ERROR Geo_Find_Cell 1091 784yrhfji348hr2ygytgfbibsvghvdcvgdscd"
            end if
            !print*, num
            !pause
            loop2:do j = 1, 4
                gran = SS%gl_Cell_gran(j, num)
                if(gran == 0) CYCLE loop2

                !print*, "Gran = ", j, " ", SS%gl_Cell_belong(:, j, num)
                !print*, "normal = ", SS%gl_Gran_normal(:, gran, 1)

                if(SS%gl_Cell_belong(1, j, num) * x + SS%gl_Cell_belong(2, j, num) * y + SS%gl_Cell_belong(3, j, num) > 0) then
                    if(present(inzone)) inzone = .False.
                    sosed = SS%gl_Gran_neighbour(1, gran)
                    if(sosed == num) sosed = SS%gl_Gran_neighbour(2, gran)

                    if(sosed > 0) then
                        num = sosed
                        cycle loop1
                    end if
                end if
            end do loop2

            EXIT loop1
        end do loop1

    end subroutine Geo_Find_Cell

    subroutine Geo_Belong_Cell(SS, x, y, num, inzone)
        ! Принадлежит ли точка ячейке
        ! num = предположительный изначальный номер (если не знаем пусть будет равен 1)
        TYPE (Setka), intent(in) :: SS
        real(8), intent(in) :: x, y
        integer(4), intent(in) :: num
        integer(4) :: j, gran, sosed, max_num
        logical, intent(out) :: inzone

        inzone = .True.

        do j = 1, 4
            gran = SS%gl_Cell_gran(j, num)
            if(gran == 0) CYCLE

            if(SS%gl_Cell_belong(1, j, num) * x + SS%gl_Cell_belong(2, j, num) * y + SS%gl_Cell_belong(3, j, num) > 0) then
                inzone = .False.
            end if
        end do


    end subroutine Geo_Belong_Cell

    subroutine Culc_Cell_Centr(SS, step)
        ! Считаем центры всех ячеек в массив
        ! step - 1 или 2  показывает какой массив координат мы меняем
        TYPE (Setka), intent(in out) :: SS
        integer(4), INTENT(IN) :: step
        integer(4) :: n1, i, a1, a2, a3, a4
        real(8) :: p1(2), p2(2), p3(2), p4(2)

        n1 = size(SS%gl_all_Cell(1, :))

        do i = 1, n1
            a1 = SS%gl_all_Cell(1, i)
            a2 = SS%gl_all_Cell(2, i)
            a3 = SS%gl_all_Cell(3, i)
            a4 = SS%gl_all_Cell(4, i)

            p1 = SS%gl_yzel(:, a1, step)
            p2 = SS%gl_yzel(:, a2, step)
            p3 = SS%gl_yzel(:, a3, step)
            p4 = SS%gl_yzel(:, a4, step)

            SS%gl_Cell_Centr(:, i, step) = (p1 + p2 + p3 + p4)/4.0
        end do
    end subroutine Culc_Cell_Centr

    subroutine Geo_Set_sxem(SS)
        ! Устанавливаем схему для грани
        TYPE (Setka), intent(in out) :: SS
        integer(4) :: N, i, s1, s2, t1, t2, yz1, yz2, j, HP
        real(8) :: c(2)

        SS%gl_Gran_shem = 2

        N = size(SS%gl_Gran_shem)
        do i = 1, N
            s1 = SS%gl_Gran_neighbour(1, i)
            s2 = SS%gl_Gran_neighbour(2, i)

            if(s1 < 1 .or. s2 < 1) then
                SS%gl_Gran_shem(i) = 3
                CYCLE
            end if

            t1 = SS%gl_all_Cell_zone(s1)
            t2 = SS%gl_all_Cell_zone(s2)

            if(t1 /= t2) then
                SS%gl_Gran_shem(i) = 3
                CYCLE
            end if

            if(t1 == 1 .or. t2 == 1) then
                SS%gl_Gran_shem(i) = 3
                CYCLE
            end if

            c = SS%gl_Gran_Center(:, i, 1)

            if(c(1) < -60) then
                SS%gl_Gran_shem(i) = 3
                CYCLE
            end if

            if(t1 == 2 .or. t2 == 2) then
                SS%gl_Gran_shem(i) = 1
                CYCLE
            end if

            if(t1 == 3 .or. t2 == 3) then
                SS%gl_Gran_shem(i) = 1
                CYCLE
            end if

            

        end do

        do j = 1, 1 !size(SS%gl_Cell_A(1, :))
            do i = SS%par_n_TS,  SS%par_n_HP - 1
                s1 = SS%gl_Cell_A(i, j)
                s2 = SS%gl_Cell_gran(3, s1)
                if(s2 > 0) SS%gl_Gran_shem(s2) = 1
                s2 = SS%gl_Cell_gran(4, s1)
                if(s2 > 0) SS%gl_Gran_shem(s2) = 1
            end do
        end do

        do j = 1, 1 !size(SS%gl_Cell_A(1, :))
            do i = SS%par_n_HP,  SS%par_n_BS - 1
                s1 = SS%gl_Cell_A(i, j)
                s2 = SS%gl_Cell_gran(3, s1)
                if(s2 > 0) SS%gl_Gran_shem(s2) = 1
                s2 = SS%gl_Cell_gran(4, s1)
                if(s2 > 0) SS%gl_Gran_shem(s2) = 1
            end do
        end do

        do j = 1, 1 !size(SS%gl_Cell_A(1, :))
            do i = SS%par_n_BS,  size(SS%gl_Cell_A(:, 1))
                s1 = SS%gl_Cell_A(i, j)
                s2 = SS%gl_Cell_gran(3, s1)
                if(s2 > 0) SS%gl_Gran_shem(s2) = 2
                s2 = SS%gl_Cell_gran(4, s1)
                if(s2 > 0) SS%gl_Gran_shem(s2) = 2
            end do
        end do

        ! do j = 1, 5 !size(SS%gl_Cell_A(1, :))
        !     do i = SS%par_n_HP - 8,  SS%par_n_HP + 7
        !         s1 = SS%gl_Cell_A(i, j)
        !         s2 = SS%gl_Cell_gran(1, s1)
        !         if(s2 > 0) SS%gl_Gran_shem(s2) = 0
        !         s2 = SS%gl_Cell_gran(2, s1)
        !         if(s2 > 0) SS%gl_Gran_shem(s2) = 0
        !     end do
        ! end do

        N = size(SS%gl_Gran_shem)
        do i = 1, N
            s1 = SS%gl_Gran_neighbour(1, i)
            s2 = SS%gl_Gran_neighbour(2, i)

            if(s1 < 1 .or. s2 < 1) then
                CYCLE
            end if

            t1 = SS%gl_all_Cell_zone(s1)
            t2 = SS%gl_all_Cell_zone(s2)

            if(t1 /= t2) then
                SS%gl_Gran_shem(i) = 3
                CYCLE
            end if
        end do

    end subroutine Geo_Set_sxem

    subroutine Geo_Time_fly(SS, XX, VV, time, cell, next)  ! Находим время time до вылета из ячейки
        TYPE (Setka), intent(in) :: SS
        real(8), intent(in out) :: XX(3)
        real(8), intent(in) :: VV(3)
        real(8), intent(out) :: time
        integer(4), intent(in) :: cell
        integer(4), intent(out) :: next

        real(8) :: A, B, C, a1, a2, a3
        real(8) :: b1, b2, b3, b4, b5, b6, b7, t1, t2, min_t, b33, r(2), rr(2)
        integer(4) :: i, gr, min_i, iter
        LOGICAL :: inzone

        iter = 0
        11   continue 
        iter = iter + 1

        min_t = 100000000000.0
        min_i = -1

        !print*, XX(1), norm2(XX(2:3))

        do i = 1, 4
            gr = SS%gl_Cell_gran(i, cell)

            if(SS%gl_Cell_neighbour(i, cell) == -4) CYCLE
            if(gr == 0) CYCLE  ! Этой грани нет в ячейке

            A = SS%gl_Cell_belong(1, i, cell)
            B = SS%gl_Cell_belong(2, i, cell)
            C = SS%gl_Cell_belong(3, i, cell)

            a3 = XX(2)**2 + XX(3)**2
            a2 = 2.0  * (XX(2) * VV(2) + XX(3) * VV(3))
            a1 = VV(2)**2 + VV(3)**2

            !print*, "abs = ", A, B, C, a1, a2, a3
            !print*, "X = ", XX
            !print*, "V - ", VV

            b1 = (A**2 * a3 * VV(1)**2 + a1 * (-a3 * B**2 + (C + A * XX(1))**2))
            b2 = a2**2 * B**2 - 4.0  * A * a2  * VV(1) * (C + A * XX(1))
            if(b2 + 4.0 * b1 < 0) CYCLE
            b33 = sqrt(b2 + 4.0 * b1)
            b3 = -a2 * B + b33
            b4 = a2 * B + b33
            b5 = 2 * A * VV(1) * (C + A * XX(1)) + B * b3
            b6 = 2 * A * VV(1) * (C + A * XX(1)) - B * b4
            b7 = 2 * a1 * B**2 - 2 * A**2 * VV(1)**2
            t1 = b5/b7
            t2 = b6/b7

            !print*, "t = ", t1, t2

            if(t1 > 0.0 .and. t1 < min_t .and. dabs(A * (XX(1) + t1 * VV(1)) + B * sqrt(a1 * t1**2 + a2 * t1 + a3) + C) < 0.00001  ) then
                min_t = t1
                min_i = i
            end if

            if(t2 > 0.0 .and. t2 < min_t .and. dabs(A * (XX(1) + t2 * VV(1)) + B * sqrt(a1 * t2**2 + a2 * t2 + a3) + C) < 0.00001  ) then
                min_t = t2
                min_i = i
            end if

        end do

        
        if(iter > 5) then
            print*, "ERROR iter > 5 345uimmnbvcw3rfbyunnrbv"
            print*, "_________________________"
            print*, XX
            print*, "_________________________ 1"
            print*, XX(1), norm2(XX(2:3))
            print*, "_________________________ 2"
            print*, SS%gl_Cell_Centr(:, cell, 1)
            print*, "_________________________"
            print*, VV
            print*, "_________________________"

            call Geo_Belong_Cell(SS, XX(1), sqrt(XX(2)**2 + XX(3)**2), cell, inzone)
            if(inzone == .False.) STOP "ERROR fhjgc2QFPcHioEkAowzgpxgP72Dfkf"

            ! Находим это время вручную
            b2 = 100.0/norm2(VV)
            b1 = 0.0

            call Geo_Belong_Cell(SS, XX(1) + b2 * VV(1), sqrt((XX(2) + b2 * VV(2))**2 + (XX(3) + b2 * VV(3))**2), cell, inzone)
            if(inzone == .True.) STOP "ERROR jmMxQH4OdkLjnuQXxVaetoqt2iX51f"

            do while( dabs(b1 - b2) > 0.00001)
                b3 = 0.5 * (b1 + b2)
                call Geo_Belong_Cell(SS, XX(1) + b3 * VV(1), sqrt((XX(2) + b3 * VV(2))**2 + (XX(3) + b3 * VV(3))**2), cell, inzone)
                if(inzone == .False.) then
                    b2 = b3
                else
                    b1 = b3
                end if
            end do

            time = b3
            next = cell
            print*, "YSPEX", b3
            return

        end if

        12 continue 

        if(min_i == -1) then
            
            !print*, "______________________________________"
            !print*, XX
            r = SS%gl_Cell_Centr(:, cell, 1)
            rr(1) = XX(1)
            rr(2) = norm2(XX(2:3))

            r = 0.995 * rr + 0.005 * r

            XX(1) = r(1)
            XX(2:3) = XX(2:3) * (r(2)/rr(2))
            !print*, XX
            !print*, "______________________________________"

            if(iter > 3) XX(2) = XX(2) + 0.0001
            GO TO 11
        end if

        t1 = min_t * 0.999
        t2 = min_t * 1.001

        call Geo_Belong_Cell(SS, XX(1) + t1 * VV(1), sqrt((XX(2) + t1 * VV(2))**2 + (XX(3) + t1 * VV(3))**2), cell, inzone)
        if(inzone /= .True.) then
            min_i = -1
            GO TO 12
        end if

        call Geo_Belong_Cell(SS, XX(1) + t2 * VV(1), sqrt((XX(2) + t2 * VV(2))**2 + (XX(3) + t2 * VV(3))**2), cell, inzone)
        if(inzone /= .False.) then
            min_i = -1
            GO TO 12
        end if

        call Geo_Belong_Cell(SS, XX(1) + 0.5 * min_t * VV(1), sqrt((XX(2) + 0.5 * min_t * VV(2))**2 + (XX(3) + 0.5 * min_t * VV(3))**2), cell, inzone)
        if(inzone /= .True.) then
            print*, "ERROR IwKfCrJlphvbfa5nXjye7CcDuiR7kT"
            min_i = -1
            GO TO 12
        end if



        
        next = SS%gl_Cell_neighbour(min_i, cell)
        time = t2

        if(next < 1) then  !! Значит не было соседа в этом направлении
            next = cell
        end if

        !print*, XX(1) + time * VV(1), norm2(XX(2:3) + time * VV(2:3))
        !pause
        return

    end subroutine Geo_Time_fly

    subroutine Geo_Culc_length_area(SS, step)
        ! Считаем длины граней и площади ячеек, используя координаты узлов на слое step
        TYPE (Setka), intent(in out) :: SS
        integer(4), INTENT(IN) :: step
        integer(4) :: i, n, y1, y2, y3, y4, j, gran
        real(8) :: p1(2), p2(2), S

        ! Считаем длины граней
        n = size(SS%gl_Gran_length(:, 1))

        do i = 1, n
            y1 = SS%gl_all_Gran(1, i)
            y2 = SS%gl_all_Gran(2, i)
            p1 = SS%gl_yzel(:, y1, step)
            p2 = SS%gl_yzel(:, y2, step)
            SS%gl_Gran_length(i, step) = norm2(p2 - p1)
            SS%gl_Gran_Center(:, i, step) = (p1 + p2)/2.0
        end do

        ! Считаем площади ячеек
        n = size(SS%gl_Cell_square(:, 1))

        do i = 1, n
            y1 = SS%gl_all_Cell(1, i)
            y2 = SS%gl_all_Cell(2, i)
            y3 = SS%gl_all_Cell(3, i)
            y4 = SS%gl_all_Cell(4, i)
            p1 = SS%gl_yzel(:, y1, step) - SS%gl_yzel(:, y3, step)
            p2 = SS%gl_yzel(:, y2, step) - SS%gl_yzel(:, y4, step)
            S = 0.5 * (p1(1) * p2(2) - p1(2) * p2(1))
            SS%gl_Cell_square(i, step) = S
            if(S <= 0.0) then
                print*, "S = ", S
                print*, "ERROR Geo_Culc_length_area 1178  ioyefdvenmlgp9yrtgb"
                !STOP
            end if

            p2 = SS%gl_Cell_Centr(:, i, step)
            do j = 1, 4
                gran = SS%gl_Cell_gran(j, i)
                if(gran == 0) then
                    SS%gl_Cell_gran_dist(j, i, step) = 0.0
                    CYCLE
                end if
                y1 = SS%gl_all_Gran(1, gran)
                y2 = SS%gl_all_Gran(2, gran)
                p1 = (SS%gl_yzel(:, y1, step) + SS%gl_yzel(:, y2, step))/2.0
                SS%gl_Cell_gran_dist(j, i, step) = norm2(p1 - p2)
            end do
        end do

    end subroutine Geo_Culc_length_area

    subroutine Geo_Culc_normal(SS, step)
        ! Считаем нормали всех ячеек
        !! Использует центры ячеек, поэтому они должны быть посчитаны
        TYPE (Setka), intent(in out) :: SS
        integer(4), INTENT(IN) :: step
        integer(4) :: N, i, a1, a2, s1, s2
        real(8) :: r1(2), r2(2), normal(2), c, norm, centr(2), vec(2), centr2(2)

        N = size(SS%gl_all_Gran(1, :))

        do i = 1, N
            a1 = SS%gl_all_Gran(1, i)
            a2 = SS%gl_all_Gran(2, i)
            r1 = SS%gl_yzel(:, a1, step)
            r2 = SS%gl_yzel(:, a2, step)

            normal = r2 - r1
            c = normal(2)
            normal(2) = -normal(1)
            normal(1) = c
            norm = norm2(normal)
            normal = normal/norm

            s1 = SS%gl_Gran_neighbour(1, i)
            s2 = SS%gl_Gran_neighbour(2, i)

            centr2 = SS%gl_Cell_Centr(:, s1, 1)
			if(s2 > 0) then
                centr = SS%gl_Cell_Centr(:, s2, 1)
			else
				centr = (r1 + r2)/2.0
			end if
            vec = centr - centr2

            if(DOT_PRODUCT(vec, normal) < 0.0) then
                normal = -normal
            end if

            SS%gl_Gran_normal(:, i, step) = normal

        end do

    end subroutine Geo_Culc_normal

    real(8) function Geo_Get_Volume_Rotate(SS, cell, angle)
        ! Угол подаётся в градусах
        TYPE (Setka), intent(in) :: SS
        integer(4), INTENT(IN) :: cell
        real(8), INTENT(IN) :: angle
        real(8) :: A, G, a1(2), a2(2), a3(2), a4(2), c1(2), c2(2), p1(2), p2(2), S1, S2, S
        integer(4) :: k, i, gr, yz1, yz2, yz3, yz4

        A = 0.0
        G = 0.0
        k = 0

        S = SS%gl_Cell_square(cell, 1)

        yz1 = SS%gl_all_Cell(1, cell)
        yz2 = SS%gl_all_Cell(2, cell)
        yz3 = SS%gl_all_Cell(3, cell)
        yz4 = SS%gl_all_Cell(4, cell)

        if(yz1 == yz4) then
            a1 = SS%gl_yzel(:, yz1, 1)
            a2 = SS%gl_yzel(:, yz2, 1)
            a3 = SS%gl_yzel(:, yz3, 1)

            c1 = (a1 + a2 + a3)/3.0
        else
            a1 = SS%gl_yzel(:, yz1, 1)
            a2 = SS%gl_yzel(:, yz2, 1)
            a3 = SS%gl_yzel(:, yz3, 1)
            a4 = SS%gl_yzel(:, yz4, 1)

            p1 = a3 - a1
            p2 = a4 - a1
            S1 = 0.5 * dabs(p1(1) * p2(2) - p1(2) * p2(1))
            p2 = a2 - a1
            S2 = 0.5 * dabs(p1(1) * p2(2) - p1(2) * p2(1))

            c1 = (a1 + a3 + a4)/3.0
            c2 = (a1 + a2 + a3)/3.0

            if(dabs(S - S1 - S2) > 0.0001) STOP "ERROR kfjoriuh9vn854ybevrw5tbyev5"

            c1 = (S1 * c1 + S2 * c2)/S
        end if

	    Geo_Get_Volume_Rotate = S * c1(2) * par_pi * angle / 180.0;   ! В градусах можно подавать угол
    end function Geo_Get_Volume_Rotate

    subroutine Geo_Find_Surface(SS)   ! Заполняет массивы поверхностей, которые выделяются
		use STORAGE
		implicit none
        TYPE (Setka), intent(in out) :: SS
		integer :: j, k, node, num, i

		! HP
		node = 1

		do j = 1, size( SS%gl_Cell_A(SS%par_n_HP - 1, :) )
			SS%gl_HP(node) = SS%gl_Cell_gran(1, SS%gl_Cell_A(SS%par_n_HP - 1, j))
			SS%gl_Gran_type(SS%gl_HP(node)) = 2
			node = node + 1
				
				
		end do

		do j = 1, size( SS%gl_Cell_C(SS%par_n_HP - SS%par_n_TS, :) )
			SS%gl_HP(node) = SS%gl_Cell_gran(1, SS%gl_Cell_C(SS%par_n_HP - SS%par_n_TS, j))
			SS%gl_Gran_type(SS%gl_HP(node)) = 2
			node = node + 1
				
		end do

		node = 1
		! TS
		do j = 1, size( SS%gl_Cell_A(SS%par_n_TS - 1, :) )
			SS%gl_TS(node) = SS%gl_Cell_gran(1, SS%gl_Cell_A(SS%par_n_TS - 1, j))
			SS%gl_Gran_type(SS%gl_TS(node)) = 1
			node = node + 1
		end do

		do j = 1, size( SS%gl_Cell_B(SS%par_n_TS - 1, :) )
			SS%gl_TS(node) = SS%gl_Cell_gran(1, SS%gl_Cell_B(SS%par_n_TS - 1, j))
			SS%gl_Gran_type(SS%gl_TS(node)) = 1
			node = node + 1
		end do


		node = 1
		! BS
		num = SS%par_m_A - 1
		do j = 1, num
			SS%gl_BS(node) = SS%gl_Cell_gran(1, SS%gl_Cell_A(SS%par_n_BS - 1, j))
			node = node + 1
		end do

    end subroutine Geo_Find_Surface

    subroutine Print_Point_from_Rays(SS)
        ! Печатаем узлы (но пробегаемся по лучам и печатаем узлы на лучах)
        ! Так проверяется связь лучей с точками
        TYPE (Setka), intent(in) :: SS
        integer :: N1, N2, i, j, n


        open(1, file = SS%name // '_Print_Point_from_Rays.txt')
        write(1,*) "TITLE = 'HP'  VARIABLES = 'X', 'Y', 'Type' "

        N1 = size(SS%gl_RAY_A(:, 1))
        N2 = size(SS%gl_RAY_A(1, :))
        do j = 1, N2
            do i = 1, N1
                n = SS%gl_RAY_A(i, j)
                if(n /= 0) then
                    write(1,*) SS%gl_yzel(:, n, 1), 1
                end if
            end do
        end do

        N1 = size(SS%gl_RAY_B(:, 1))
        N2 = size(SS%gl_RAY_B(1, :))
        do j = 1, N2
            do i = 1, N1
                n = SS%gl_RAY_B(i, j)
                if(n /= 0) then
                    write(1,*) SS%gl_yzel(:, n, 1), 2
                end if
            end do
        end do

        N1 = size(SS%gl_RAY_C(:, 1))
        N2 = size(SS%gl_RAY_C(1, :))
        do j = 1, N2
            do i = 1, N1
                n = SS%gl_RAY_C(i, j)
                if(n /= 0) then
                    write(1,*) SS%gl_yzel(:, n, 1), 3
                end if
            end do
        end do

        N1 = size(SS%gl_RAY_O(:, 1))
        N2 = size(SS%gl_RAY_O(1, :))
        do j = 1, N2
            do i = 1, N1
                n = SS%gl_RAY_O(i, j)
                if(n /= 0) then
                    write(1,*) SS%gl_yzel(:, n, 1), 4
                end if
            end do
        end do

        N1 = size(SS%gl_RAY_K(:, 1))
        N2 = size(SS%gl_RAY_K(1, :))
        do j = 1, N2
            do i = 1, N1
                n = SS%gl_RAY_K(i, j)
                if(n /= 0) then
                    write(1,*) SS%gl_yzel(:, n, 1), 5
                end if
            end do
        end do

        N1 = size(SS%gl_RAY_D(:, 1))
        N2 = size(SS%gl_RAY_D(1, :))
        do j = 1, N2
            do i = 1, N1
                n = SS%gl_RAY_D(i, j)
                if(n /= 0) then
                    write(1,*) SS%gl_yzel(:, n, 1), 6
                end if
            end do
        end do

        N1 = size(SS%gl_RAY_E(:, 1))
        N2 = size(SS%gl_RAY_E(1, :))
        do j = 1, N2
            do i = 1, N1
                n = SS%gl_RAY_E(i, j)
                if(n /= 0) then
                    write(1,*) SS%gl_yzel(:, n, 1), 7
                end if
            end do
        end do

        close(1)

    end subroutine Print_Point_from_Rays

    subroutine Print_Cell_Centr(SS)
        ! Печатаем центры всех ячеек
        TYPE (Setka), intent(in) :: SS
        integer :: i, j, node


        open(1, file = SS%name // '_Print_Cell_Centr.txt')
        write(1,*) "TITLE = 'HP'  VARIABLES = 'X', 'Y', 'Type', 'zone' "

        do j = 1, size(SS%gl_Cell_A(1, :))
            do i = 1, size(SS%gl_Cell_A(:, 1))
                node = SS%gl_Cell_A(i, j)
                write(1,*) SS%gl_Cell_Centr(:, node, 1), 1, SS%gl_all_Cell_zone(node)
            end do
        end do

        do j = 1, size(SS%gl_Cell_B(1, :))
            do i = 1, size(SS%gl_Cell_B(:, 1))
                node = SS%gl_Cell_B(i, j)
                write(1,*) SS%gl_Cell_Centr(:, node, 1), 2, SS%gl_all_Cell_zone(node)
            end do
        end do

        do j = 1, size(SS%gl_Cell_C(1, :))
            do i = 1, size(SS%gl_Cell_C(:, 1))
                node = SS%gl_Cell_C(i, j)
                write(1,*) SS%gl_Cell_Centr(:, node, 1), 3, SS%gl_all_Cell_zone(node)
            end do
        end do

        close(1)

    end subroutine Print_Cell_Centr

    subroutine Geo_request(SS)
        ! Печатаем центры всех ячеек
        ! Для того, чтобы считать в другой программе и вернуть значения газодинамических перменных
        TYPE (Setka), intent(in) :: SS
        integer :: i, j, node


        open(1, file = "request.bin", FORM = 'BINARY')

        write(1) size(SS%gl_all_Cell(1, :))

        do j = 1, size(SS%gl_all_Cell(1, :))
            write(1) SS%gl_Cell_Centr(:, j, 1)
        end do

        close(1)

    end subroutine Geo_request

    subroutine Geo_get_request(SS)
        ! Печатаем центры всех ячеек
        ! Для того, чтобы считать в другой программе и вернуть значения газодинамических перменных
        TYPE (Setka), intent(in out) :: SS
        integer :: i, j, n
        logical :: exists
        real(8) :: ro, u, v, p

        inquire(file="answer_request.bin", exist=exists)

        open(1, file = "answer_request.bin", FORM = 'BINARY', ACTION = "READ")

        read(1) n

        do j = 1, size(SS%gl_all_Cell(1, :))
            read(1) ro, u, v, p
            if(ro <= 0.0) then
                ro = 0.000001
            end if
            if(p <= 0.0) then
                p = 0.000001
            end if
            SS%gd(1, j, 1) = ro
            SS%gd(2, j, 1) = p
            SS%gd(3, j, 1) = u
            SS%gd(4, j, 1) = v

            read(1) ro, u, v, p
            if(ro <= 0.0) then
                ro = 0.000001
            end if
            if(p <= 0.0) then
                p = 0.000001
            end if
            SS%hydrogen(1, 1, j, 1) = ro
            SS%hydrogen(2, 1, j, 1) = p
            SS%hydrogen(3, 1, j, 1) = u
            SS%hydrogen(4, 1, j, 1) = v

            read(1) ro, u, v, p
            if(ro <= 0.0) then
                ro = 0.000001
            end if
            if(p <= 0.0) then
                p = 0.000001
            end if
            SS%hydrogen(1, 2, j, 1) = ro
            SS%hydrogen(2, 2, j, 1) = p
            SS%hydrogen(3, 2, j, 1) = u
            SS%hydrogen(4, 2, j, 1) = v

            read(1) ro, u, v, p
            if(ro <= 0.0) then
                ro = 0.000001
            end if
            if(p <= 0.0) then
                p = 0.000001
            end if
            SS%hydrogen(1, 3, j, 1) = ro
            SS%hydrogen(2, 3, j, 1) = p
            SS%hydrogen(3, 3, j, 1) = u
            SS%hydrogen(4, 3, j, 1) = v

            read(1) ro, u, v, p
            if(ro <= 0.0) then
                ro = 0.000001
            end if
            if(p <= 0.0) then
                p = 0.000001
            end if
            SS%hydrogen(1, 4, j, 1) = ro
            SS%hydrogen(2, 4, j, 1) = p
            SS%hydrogen(3, 4, j, 1) = u
            SS%hydrogen(4, 4, j, 1) = v

            SS%hydrogen(:, :, j, 2) = SS%hydrogen(:, :, j, 1)
            SS%gd(:, j, 2) = SS%gd(:, j, 1)
        end do

        close(1)

    end subroutine Geo_get_request

    subroutine Print_GD(SS)
        ! Печатаем центры всех ячеек
        TYPE (Setka), intent(in) :: SS
        integer :: i, j, node
        real(8) :: Mach

        open(1, file = SS%name // '_Print_GD.txt')
        write(1,*) "TITLE = 'HP'  VARIABLES = X, Y, Rho, p, u, v, Q, Mach"

        do j = 1, size(SS%gl_all_Cell(1, :))
            Mach = norm2(SS%gd(3:4, j, 1))/sqrt(SS%par_ggg * SS%gd(2, j, 1)/SS%gd(1, j, 1))
            write(1,*) SS%gl_Cell_Centr(:, j, 1), SS%gd(1:4, j, 1), SS%gd(5, j, 1)/SS%gd(1, j, 1), Mach
        end do

        close(1)

    end subroutine Print_GD

    subroutine Print_hydrogen(SS)
        ! Печатаем центры всех ячеек
        TYPE (Setka), intent(in) :: SS
        integer :: i, j, node
        real(8) :: Mach

        open(1, file = SS%name // '_Print_hydrogen.txt')
        write(1,*) "TITLE = 'HP'  VARIABLES = X, Y" 
        write(1,*)", Rho1, p1, u1, v1, T1"
        write(1,*)", Rho2, p2, u2, v2, T2"
        write(1,*)", Rho3, p3, u3, v3, T3"
        write(1,*)", Rho4, p4, u4, v4, T4"
        write(1,*)", k_u, k_v, k_T, In"

        do j = 1, size(SS%gl_all_Cell(1, :))
            write(1,*) SS%gl_Cell_Centr(:, j, 1), SS%hydrogen(1:5, 1, j, 1), &
            SS%hydrogen(1:5, 2, j, 1), SS%hydrogen(1:5, 3, j, 1), SS%hydrogen(1:5, 4, j, 1), SS%atom_source(1:4, j)
        end do

        close(1)

    end subroutine Print_hydrogen

    subroutine Print_GD_1D(SS)
        ! Печатаем центры всех ячеек
        TYPE (Setka), intent(in) :: SS
        integer :: i, j, node
        real(8) :: Mach, c(2)

        open(1, file = SS%name // '_Print_GD_1D.txt')
        write(1,*) "TITLE = 'HP'  VARIABLES = X, Rho, p, u, v, Q, Mach"

        do i = size(SS%gl_Cell_B(:, 1)), 1, -1
            j = SS%gl_Cell_B(i, 1)
            Mach = norm2(SS%gd(3:4, j, 1))/sqrt(SS%par_ggg * SS%gd(2, j, 1)/SS%gd(1, j, 1))
			c = SS%gl_Cell_Centr(1, j, 1)
            write(1,*) SS%gl_Cell_Centr(1, j, 1), SS%gd(1:4, j, 1), SS%gd(5, j, 1)/SS%gd(1, j, 1), Mach
        end do

        do i = 1, size(SS%gl_Cell_A(:, 1))
            j = SS%gl_Cell_A(i, 1)
            Mach = norm2(SS%gd(3:4, j, 1))/sqrt(SS%par_ggg * SS%gd(2, j, 1)/SS%gd(1, j, 1))
            write(1,*) SS%gl_Cell_Centr(1, j, 1), SS%gd(1:4, j, 1), SS%gd(5, j, 1)/SS%gd(1, j, 1), Mach
        end do

        close(1)

    end subroutine Print_GD_1D

    subroutine Print_hydrogen_1D(SS)
        ! Печатаем центры всех ячеек
        TYPE (Setka), intent(in) :: SS
        integer :: i, j, node

        if(SS%n_Hidrogen /= 4) then
            print*, "Error n_Hidrogen /= 4   Print_hydrogen_1D  1906  y87wtrguwvby8bej7yt7vwc8yynu6rbv"
        end if

        open(1, file = SS%name // '_Print_hydrogen_1D.txt')
        write(1,*) "TITLE = 'HP'  VARIABLES = X" 
        write(1,*)",Rho1, p1, u1, v1, T1"
        write(1,*)",Rho2, p2, u2, v2, T2"
        write(1,*)",Rho3, p3, u3, v3, T3"
        write(1,*)",Rho4, p4, u4, v4, T4"
        write(1,*)", k1, k2, k3"

        do i = size(SS%gl_Cell_B(:, 1)), 1, -1
            j = SS%gl_Cell_B(i, 1)
            
            write(1,*) SS%gl_Cell_Centr(1, j, 1), SS%hydrogen(1:5, 1, j, 1), &
                SS%hydrogen(1:5, 2, j, 1), SS%hydrogen(1:5, 3, j, 1), SS%hydrogen(1:5, 4, j, 1), SS%atom_source(1:3, j)
        end do

        do i = 1, size(SS%gl_Cell_A(:, 1))
            j = SS%gl_Cell_A(i, 1)

            write(1,*) SS%gl_Cell_Centr(1, j, 1), SS%hydrogen(1:5, 1, j, 1), &
            SS%hydrogen(1:5, 2, j, 1), SS%hydrogen(1:5, 3, j, 1), SS%hydrogen(1:5, 4, j, 1), SS%atom_source(1:3, j)
        end do

        close(1)

    end subroutine Print_hydrogen_1D

    subroutine Print_Cell(SS)
        ! Печатаем все ячейки (есть дублирование), каждая ячейка печатается отдельно
        TYPE (Setka), intent(in) :: SS
        integer :: N1, N2, i, j, N, node, yz

        N = size(SS%gl_Cell_A) + size(SS%gl_Cell_B) + size(SS%gl_Cell_C)
        open(1, file = SS%name // '_Print_Cell.txt')
        write(1,*) "TITLE = 'HP'  VARIABLES = 'X', 'Y'  ZONE T= 'HP', N= ", 4 * N, ", E =  ", 4 * N , ", F=FEPOINT, ET=LINESEG"

        do j = 1, size(SS%gl_Cell_A(1, :))
            do i = 1, size(SS%gl_Cell_A(:, 1))
                node = SS%gl_Cell_A(i, j)
                yz = SS%gl_all_Cell(1, node)
                write(1,*) SS%gl_yzel(:, yz, 1)
                yz = SS%gl_all_Cell(2, node)
                write(1,*) SS%gl_yzel(:, yz, 1)
                yz = SS%gl_all_Cell(3, node)
                write(1,*) SS%gl_yzel(:, yz, 1)
                yz = SS%gl_all_Cell(4, node)
                write(1,*) SS%gl_yzel(:, yz, 1)
            end do
        end do

        do j = 1, size(SS%gl_Cell_B(1, :))
            do i = 1, size(SS%gl_Cell_B(:, 1))
                node = SS%gl_Cell_B(i, j)
                yz = SS%gl_all_Cell(1, node)
                write(1,*) SS%gl_yzel(:, yz, 1)
                yz = SS%gl_all_Cell(2, node)
                write(1,*) SS%gl_yzel(:, yz, 1)
                yz = SS%gl_all_Cell(3, node)
                write(1,*) SS%gl_yzel(:, yz, 1)
                yz = SS%gl_all_Cell(4, node)
                write(1,*) SS%gl_yzel(:, yz, 1)
            end do
        end do

        do j = 1, size(SS%gl_Cell_C(1, :))
            do i = 1, size(SS%gl_Cell_C(:, 1))
                node = SS%gl_Cell_C(i, j)
                yz = SS%gl_all_Cell(1, node)
                write(1,*) SS%gl_yzel(:, yz, 1)
                yz = SS%gl_all_Cell(2, node)
                write(1,*) SS%gl_yzel(:, yz, 1)
                yz = SS%gl_all_Cell(3, node)
                write(1,*) SS%gl_yzel(:, yz, 1)
                yz = SS%gl_all_Cell(4, node)
                write(1,*) SS%gl_yzel(:, yz, 1)
            end do
        end do
        
        do j = 0, N
            write(1,*) 4 * j + 1, 4 * j + 2
            write(1,*) 4 * j + 2, 4 * j + 3
            write(1,*) 4 * j + 3, 4 * j + 4
            write(1,*) 4 * j + 4, 4 * j + 1
        end do

        close(1)

    end subroutine Print_Cell

    subroutine Print_Connect(SS)
        ! Печатаем все связи у ячеек
        TYPE (Setka), intent(in) :: SS
        integer(4) :: N, i, j, sosed
        real(8) :: norm
        
        N = size(SS%gl_Cell_A) + size(SS%gl_Cell_B) + size(SS%gl_Cell_C)

        open(1, file = SS%name // '_Print_Connect_Cell.txt')
        write(1,*) "TITLE = 'HP'  VARIABLES = 'X', 'Y'  ZONE T= 'HP', N= ", 5 * N, ", E =  ", 4 * N , ", F=FEPOINT, ET=LINESEG"

        do i = 1, N
            write(1,*) SS%gl_Cell_Centr(:, i, 1)
            do j = 1, 4
                sosed = SS%gl_Cell_neighbour(j, i)
                if(sosed >= 1) then
                    write(1,*) (SS%gl_Cell_Centr(:, sosed, 1) + SS%gl_Cell_Centr(:, i, 1))/2.0
                else if (sosed == -1) then
                    norm = norm2(SS%gl_Cell_Centr(:, i, 1))
                    write(1,*) SS%gl_Cell_Centr(:, i, 1) * SS%par_R_END/norm
                else if (sosed == -2) then
                    write(1,*) SS%par_R_LEFT, SS%gl_Cell_Centr(2, i, 1)
                else if (sosed == -3) then
                    write(1,*) SS%gl_Cell_Centr(1, i, 1), SS%par_R_END
                else if (sosed == -4) then
                    write(1,*) SS%gl_Cell_Centr(1, i, 1), 0.0
                else
                    write(1,*) SS%gl_Cell_Centr(:, i, 1)
                end if
            end do
        end do

        do j = 0, N
            write(1,*) 5 * j + 1, 5 * j + 2
            write(1,*) 5 * j + 1, 5 * j + 3
            write(1,*) 5 * j + 1, 5 * j + 4
            write(1,*) 5 * j + 1, 5 * j + 5
        end do

        close(1)

    end subroutine Print_Connect

    subroutine Print_Grans(SS, number)
        ! Печатаем все ячейки (есть дублирование), каждая ячейка печатается отдельно
        TYPE (Setka), intent(in) :: SS
        integer(4), intent(in), OPTIONAL :: number
        integer(4) :: N, i, a1, a2, j, num
        character(len=5) :: name

        num = 0
        if(PRESENT(number)) num = number
        write(unit=name,fmt='(i5.5)') num

        N = size(SS%gl_all_Gran(1, :))

        open(1, file = 'Print_Grans_' // name // '.txt')
        write(1,*) "TITLE = 'HP'  VARIABLES = 'X', 'Y', 'Sxem'  ZONE T= 'HP', N= ", 2 * N, ", E =  ", N , ", F=FEPOINT, ET=LINESEG"

        do i = 1, N
            a1 = SS%gl_all_Gran(1, i)
            a2 = SS%gl_all_Gran(2, i)
            write(1,*) SS%gl_yzel(:, a1, 1), SS%gl_Gran_shem(i)
            write(1,*) SS%gl_yzel(:, a2, 1), SS%gl_Gran_shem(i)
        end do

        do j = 0, N
            write(1,*) 2 * j + 1, 2 * j + 2
        end do

        close(1)
    end subroutine Print_Grans

    subroutine Print_TVD_Sosed(SS)
        ! Печатаем все ячейки (есть дублирование), каждая ячейка печатается отдельно
        TYPE (Setka), intent(in) :: SS
        integer(4) :: N, i, s1, s2, j
        real(8) :: c1(2), c2(2), c3(2)


        N = size(SS%gl_all_Gran(1, :))
        open(1, file = 'Print_TVD_sosed.txt')
        write(1,*) "TITLE = 'HP'  VARIABLES = 'X', 'Y'  ZONE T= 'HP', N= ", 3 * N, ", E =  ", 2 * N , ", F=FEPOINT, ET=LINESEG"

        do i = 1, N
            c1 = SS%gl_Gran_Center(:, i, 1)
            s1 = SS%gl_Gran_neighbour_TVD(1, i)
            s2 = SS%gl_Gran_neighbour_TVD(2, i)

            if(s1 > 0) then
                c2 = SS%gl_Cell_Centr(:, s1, 1)
            else
                c2 = c1
            end if

            if(s2 > 0) then
                c3 = SS%gl_Cell_Centr(:, s2, 1)
            else
                c3 = c1
            end if

            write(1, *) c1, c2, c3
        end do

        do j = 0, N
            write(1,*) 3 * j + 1, 3 * j + 2
            write(1,*) 3 * j + 1, 3 * j + 3
        end do

        close(1)
    end subroutine Print_TVD_Sosed

    subroutine Geo_Print_Surface(SS, number)
        ! Печатаем поверхности, которые выделяем
        TYPE (Setka), intent(in) :: SS
        integer(4), intent(in), OPTIONAL :: number
        integer(4) :: i, N, a1, a2, gr, j, num
        character(len=5) :: name

        num = 0
        if(PRESENT(number)) num = number
        write(unit=name,fmt='(i5.5)') num


        N = size(SS%gl_HP) + size(SS%gl_TS) + size(SS%gl_BS)

        open(1, file = SS%name // "_Print_Surface_" // name // ".txt")
        write(1,*) "TITLE = 'HP'  VARIABLES = 'X', 'Y'  ZONE T= 'HP', N= ", 2 * N, ", E =  ", N , ", F=FEPOINT, ET=LINESEG"

        do i = 1, size(SS%gl_HP)
            gr = SS%gl_HP(i)
            a1 = SS%gl_all_Gran(1, gr)
            a2 = SS%gl_all_Gran(2, gr)
            write(1, *) SS%gl_yzel(:, a1, 1)
            write(1, *) SS%gl_yzel(:, a2, 1)
        end do

        do i = 1, size(SS%gl_TS)
            gr = SS%gl_TS(i)
            a1 = SS%gl_all_Gran(1, gr)
            a2 = SS%gl_all_Gran(2, gr)
            write(1, *) SS%gl_yzel(:, a1, 1)
            write(1, *) SS%gl_yzel(:, a2, 1)
        end do

        do i = 1, size(SS%gl_BS)
            gr = SS%gl_BS(i)
            a1 = SS%gl_all_Gran(1, gr)
            a2 = SS%gl_all_Gran(2, gr)
            write(1, *) SS%gl_yzel(:, a1, 1)
            write(1, *) SS%gl_yzel(:, a2, 1)
        end do

        do j = 0, N
            write(1,*) 2 * j + 1, 2 * j + 2
        end do

        close(1)
	end subroutine Geo_Print_Surface

    subroutine Proverka_grans_sosed(SS)
        ! Автоматическая проверка граней (их соседей) и соседей (в ячейках)
        TYPE (Setka), intent(in) :: SS
        integer(4) :: i, j, sosed, gran, N, s1, s2
        real(8) :: normal(2), gran_center(2), gran_center2(2), z1(2), z2(2), the1, the2

        print*, "Proverka_grans_sosed start"

        N = size(SS%gl_Cell_Centr(1, :, 1))

        do i = 1, N  ! Пробегаемся по всем ячейкам
            do j = 1, 4
                sosed = SS%gl_Cell_neighbour(j, i)
                gran = SS%gl_Cell_gran(j, i)

                if(sosed == 0 .and. gran /= 0) STOP "ERROR Proverka_grans_sosed 1312 78iuygbnkuwf"  
                if(sosed /= 0 .and. gran == 0) STOP "ERROR Proverka_grans_sosed 1313 90olkjhgfdewazxcft"  

                if(sosed /= 0) then
                    if(sosed /= SS%gl_Gran_neighbour(1, gran) .and. sosed /= SS%gl_Gran_neighbour(2, gran)) STOP "ERROR Proverka_grans_sosed 1316 456yuhgrewsdhjkiu"  
                end if
            end do
        end do

        ! Проверка на то, чтобы нормали у поверхностей разрыва были по умолчанию ориентированы наружу

        N = size(SS%gl_TS)
        do i = 1, N
            gran = SS%gl_TS(i)
            normal = SS%gl_Gran_normal(:, gran, 1)
            gran_center = SS%gl_Gran_Center(:, gran, 1)

            if(DOT_PRODUCT(normal, gran_center) < 0.0) then
                print*, "Error TS normal Proverka_grans_sosed oihiuergfuevguebgeg"
                pause
                STOP
            end if
        end do

        N = size(SS%gl_HP)
        do i = 1, N
            gran = SS%gl_HP(i)
            normal = SS%gl_Gran_normal(:, gran, 1)
            gran_center = SS%gl_Gran_Center(:, gran, 1)
            gran_center2(1) = -120.0
            gran_center2(2) = -50.0

            if(DOT_PRODUCT(normal, gran_center - gran_center2) < 0.0) then
                print*, "Error HP normal Proverka_grans_sosed 08xcvgfdertyukopuy"
                pause
                STOP
            end if
        end do

        N = size(SS%gl_BS)
        do i = 1, N
            gran = SS%gl_BS(i)
            normal = SS%gl_Gran_normal(:, gran, 1)
            gran_center = SS%gl_Gran_Center(:, gran, 1)
            gran_center2(1) = -100.0
            gran_center2(2) = 0.0

            if(DOT_PRODUCT(normal, gran_center - gran_center2) < 0.0) then
                print*, "Error BS normal Proverka_grans_sosed ijtrogeuyvghcvhythjyt"
                pause
                STOP
            end if
        end do


        ! Проверка, чтобы узлы у граней-поверхностей шли в правильном порядке
        N = size(SS%gl_TS)
        do i = 1, N
            gran = SS%gl_TS(i)
            s1 = SS%gl_all_Gran(1, gran)
            s2 = SS%gl_all_Gran(2, gran)

            z1 = SS%gl_yzel(:, s1, 1)
            z2 = SS%gl_yzel(:, s2, 1)
            
            the1 = polar_angle(z1(1), z1(2))
            the2 = polar_angle(z2(1), z2(2))

            if(the2 < the1) then
                print*, "Error 1869 Proverka_grans_sosed  8w48tvw78wajv08waa"
                STOP
            end if
        end do

        ! Проверка, чтобы узлы у граней-поверхностей шли в правильном порядке
        N = size(SS%gl_HP)
        do i = 1, N
            gran = SS%gl_HP(i)
            s1 = SS%gl_all_Gran(1, gran)
            s2 = SS%gl_all_Gran(2, gran)

            z1 = SS%gl_yzel(:, s1, 1)
            z2 = SS%gl_yzel(:, s2, 1)
            
            the1 = polar_angle(z1(1), z1(2))
            the2 = polar_angle(z2(1), z2(2))

            if(the2 < the1) then
                print*, "Error 1888 Proverka_grans_sosed  lkfyrguyrdvsaefvrs"
                STOP
            end if
        end do

        ! Проверка, чтобы узлы у граней-поверхностей шли в правильном порядке
        N = size(SS%gl_BS)
        do i = 1, N
            gran = SS%gl_BS(i)
            s1 = SS%gl_all_Gran(1, gran)
            s2 = SS%gl_all_Gran(2, gran)

            z1 = SS%gl_yzel(:, s1, 1)
            z2 = SS%gl_yzel(:, s2, 1)
            
            the1 = polar_angle(z1(1), z1(2))
            the2 = polar_angle(z2(1), z2(2))

            if(the2 < the1) then
                print*, "Error 1907 Proverka_grans_sosed  56ub5345vwcqtvbyer"
                STOP
            end if
        end do


        print*, "Proverka_grans_sosed end"

    end subroutine Proverka_grans_sosed

    real(8) pure function sgushenie_2(x, all)
        ! x от 0 до 1 и возвращает функция от 0 до 1
        ! Сгущение точек к обоим концам отрезка (сильнее, чем предыдущая функция)
        implicit none
        real(8), intent(in) :: x, all

        sgushenie_2 = all * x - 10 * (-1 + all) * x**3 + 15 * (-1 + all) * x**4 - 6 * (-1 + all) * x**5
        
        return
    end function sgushenie_2

    real(8) pure function sgushenie_3(x, all)
        ! x от 0 до 1 и возвращает функция от 0 до 1
        ! Сгущение точек к 1
        implicit none
        real(8), intent(in) :: x, all
        sgushenie_3 = -(-x + 1)**all + 1.0
        return
    end function sgushenie_3
        
    real(8) pure function sgushenie_4(x, x0, y0)
        ! x от 0 до 1 и возвращает функция от 0 до 1
        ! Сгущение точек за BS
        implicit none
        real(8), intent(in) :: x, x0, y0
        
        if(x < x0) then
            sgushenie_4 = (x - x0) * y0/x0 + y0
        else
            sgushenie_4 = (1.0 - y0) * (x - x0)/(1.0 - x0) + y0
        end if
        
        return
    end function sgushenie_4

    subroutine Geo_culc_TVD_sosed(SS)
        TYPE (Setka), intent(in out) :: SS
        integer(4) :: Num, gr, cell, j, cell2, best
        real(8) :: normal(2), c1(2), c2(2), vec(2), xc, dotmax

        SS%gl_Gran_neighbour_TVD = 0

        Num = size(SS%gl_all_Gran(1, :))

        ! Пробегаемся по первым соседям каждой грани
        do gr = 1, Num
            normal = SS%gl_Gran_normal(:, gr, 1)
            cell = 0
            cell2 = 0
            dotmax = 0.0
            best = 0

            cell = SS%gl_Gran_neighbour(1, gr) 
            if(cell < 1) CYCLE

            c1 = SS%gl_Cell_Centr(:, cell, 1)

            do j = 1, 4
                cell2 = SS%gl_Cell_neighbour(j, cell)
                if(cell2 < 1) CYCLE
                c2 = SS%gl_Cell_Centr(:, cell2, 1)
                vec = c2 - c1
				vec = vec/norm2(vec)
                xc = DOT_PRODUCT(vec, -normal)
				if (dotmax < xc) then
					dotmax = xc
					best = j
				end if
            end do

			if (best /= 0 .and. dotmax > 0.3) then
				SS%gl_Gran_neighbour_TVD(1, gr) = SS%gl_Cell_neighbour(best, cell)
			end if

        end do

        ! Пробегаемся по вторым соседям каждой грани
        do gr = 1, Num
            normal = SS%gl_Gran_normal(:, gr, 1)
            cell = 0
            cell2 = 0
            dotmax = 0
            best = 0

            cell = SS%gl_Gran_neighbour(2, gr) 
            if(cell < 1) CYCLE

            c1 = SS%gl_Cell_Centr(:, cell, 1)

            do j = 1, 4
                cell2 = SS%gl_Cell_neighbour(j, cell)
                if(cell2 < 1) CYCLE
                c2 = SS%gl_Cell_Centr(:, cell2, 1)
                vec = c2 - c1
				vec = vec/norm2(vec)
                xc = DOT_PRODUCT(vec, normal)
				if (dotmax < xc) then
					dotmax = xc
					best = j
				end if
            end do

			if (best /= 0 .and. dotmax > 0.3) then
				SS%gl_Gran_neighbour_TVD(2, gr) = SS%gl_Cell_neighbour(best, cell)
			end if

        end do


    end subroutine Geo_culc_TVD_sosed

    subroutine Geo_Culc_zone(SS)
        TYPE (Setka), intent(in out) :: SS
        integer(4) :: N1, N2, i, j, cell

        ! Нужно распределить зоны ячеек
        SS%gl_all_Cell_zone = 4
        N2 = size(SS%gl_Cell_A(1, :))
        N1 = size(SS%gl_Cell_A(:, 1))
        do j = 1, N2
            do i = 1, N1
                cell = SS%gl_Cell_A(i, j)
                if(i < SS%par_n_TS) then
                    SS%gl_all_Cell_zone(cell) = 1 
                else if (i < SS%par_n_HP) then
                    SS%gl_all_Cell_zone(cell) = 2
                else if (i < SS%par_n_BS) then
                    SS%gl_all_Cell_zone(cell) = 3
                end if
            end do
        end do

        N2 = size(SS%gl_Cell_B(1, :))
        N1 = size(SS%gl_Cell_B(:, 1))
        do j = 1, N2
            do i = 1, N1
                cell = SS%gl_Cell_B(i, j)
                if(i < SS%par_n_TS) then
                    SS%gl_all_Cell_zone(cell) = 1 
                else 
                    SS%gl_all_Cell_zone(cell) = 2
                end if
            end do
        end do

        N2 = size(SS%gl_Cell_C(1, :))
        N1 = size(SS%gl_Cell_C(:, 1))
        do j = 1, N2
            do i = 1, N1
                cell = SS%gl_Cell_C(i, j)
                if(i < SS%par_n_HP - SS%par_n_TS + 1) then
                    SS%gl_all_Cell_zone(cell) = 2 
                else if(i < SS%par_n_BS - SS%par_n_TS + 1) then
                    SS%gl_all_Cell_zone(cell) = 3
                end if
            end do
        end do
    end subroutine Geo_Culc_zone

    subroutine Save_setka_bin(SS, name)  ! Сохранение сетки в бинарном файле
        TYPE (Setka), intent(in) :: SS
        CHARACTER(len = 5), intent(in) :: name

        open(1, file = "Save_all_" // name // ".bin", FORM = 'BINARY')

        write(1) SS%name

        write(1) SS%par_m_A
        write(1) SS%par_m_BC
        write(1) SS%par_m_O
        write(1) SS%par_m_K
        write(1) SS%par_triple_point
        write(1) SS%par_triple_point_2
        
        write(1) SS%par_n_TS 
        write(1) SS%par_n_HP 
        write(1) SS%par_n_BS 
        write(1) SS%par_n_END
        write(1) SS%par_n_IA 
        write(1) SS%par_n_IB 

        write(1) SS%par_R_character
        write(1) SS%par_R0
        write(1) SS%par_R_END
        write(1) SS%par_R_LEFT
        write(1) SS%par_R_inner

        !! Физические параметры ---------------------------------------
        write(1) SS%par_a_2 
        write(1) SS%par_ggg 
        write(1) SS%par_Velosity_inf
        write(1) SS%par_n_H_LISM 
        write(1) SS%par_Kn
        !! -------------------------------------------------------------

        !Набор параметров сгущения
        write(1) SS%par_kk1
        write(1) SS%par_kk2
        write(1) SS%par_kk3
        write(1) SS%par_kk31 
        write(1) SS%par_kk13
        write(1) SS%par_kk131
        write(1) SS%par_kk132
        write(1) SS%par_kk14  
        write(1) SS%par_kk12

        write(1) SS%par_n_points                         ! Всего точек в сетке (считается при инициализации сетки)

        write(1) SS%gl_yzel

        write(1) SS%gl_RAY_A
        write(1) SS%gl_RAY_B
        write(1) SS%gl_RAY_C
        write(1) SS%gl_RAY_O
        write(1) SS%gl_RAY_K
        write(1) SS%gl_RAY_D
        write(1) SS%gl_RAY_E

        ! Ячейки
        write(1) SS%gl_Cell_A
        write(1) SS%gl_Cell_B
        write(1) SS%gl_Cell_C

        write(1) SS%gl_all_Cell

        write(1) SS%gl_Cell_neighbour

        write(1) SS%gl_Cell_gran
        write(1) SS%gl_Cell_gran_dist

        write(1) SS%gl_Cell_Centr
        
        write(1) SS%gl_all_Gran
        write(1) SS%gl_Gran_neighbour
        write(1) SS%gl_Gran_normal
        write(1) SS%gl_Gran_length                     
        write(1) SS%gl_Gran_Center                      
        
        write(1) SS%gl_Cell_belong
        write(1) SS%gl_Cell_square

        write(1) SS%gl_Cell_type
        write(1) SS%gl_Cell_number

        ! Поверхности выделения
        write(1) SS%gl_HP
        write(1) SS%gl_TS
        write(1) SS%gl_BS

        write(1) SS%gl_Gran_type     

        !! ФИЗИКА
        write(1) SS%n_Hidrogen
        write(1) SS%n_par

        write(1) SS%gd
        write(1) SS%hydrogen


        write(1) 1;
        write(1) size(SS%atom_all_source(:, 1, 1)), size(SS%atom_all_source(1, :, 1)), size(SS%atom_all_source(1, 1, :))
        write(1) SS%atom_all_source
        write(1) size(SS%atom_source(:, 1)), size(SS%atom_source(1, :))
        write(1) SS%atom_source
        
        
        write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0
        write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0
        write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0
        write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0
        write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0
        write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0
        write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0
        write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0
        write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0
        write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0
        write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0
        write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0
        write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0
        write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0
        write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0
        write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0
        write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0
        write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0
        write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0; write(1) 0

        close(1)

    end subroutine Save_setka_bin

    subroutine Read_setka_bin(SS, name)  ! Сохранение сетки в бинарном файле
        TYPE (Setka), intent(in out) :: SS
        CHARACTER(len = 5), intent(in) :: name
        logical :: exists
        integer(4) :: i, j, k, n, n1, n2, n3

        inquire(file= "Save_all_" // name // ".bin", exist=exists)
    
        if (exists == .False.) then
            print*, "net faila 1898 tfgdhfwy4rfetrgfw4rwter!!!"
            STOP "net faila!!!"
        end if

        open(1, file = "Save_all_" // name // ".bin", FORM = 'BINARY', ACTION = "READ")

        read(1) SS%name

        read(1) SS%par_m_A
        read(1) SS%par_m_BC
        read(1) SS%par_m_O
        read(1) SS%par_m_K
        read(1) SS%par_triple_point
        read(1) SS%par_triple_point_2
        
        read(1) SS%par_n_TS 
        read(1) SS%par_n_HP 
        read(1) SS%par_n_BS 
        read(1) SS%par_n_END
        read(1) SS%par_n_IA 
        read(1) SS%par_n_IB 

        read(1) SS%par_R_character
        read(1) SS%par_R0
        read(1) SS%par_R_END
        read(1) SS%par_R_LEFT
        read(1) SS%par_R_inner

        !! Физические параметры ---------------------------------------
        read(1) SS%par_a_2 
        read(1) SS%par_ggg 
        read(1) SS%par_Velosity_inf
        read(1) SS%par_n_H_LISM 
        read(1) SS%par_Kn
        !! -------------------------------------------------------------

        !Набор параметров сгущения
        read(1) SS%par_kk1
        read(1) SS%par_kk2
        read(1) SS%par_kk3
        read(1) SS%par_kk31 
        read(1) SS%par_kk13
        read(1) SS%par_kk131
        read(1) SS%par_kk132
        read(1) SS%par_kk14  
        read(1) SS%par_kk12

        call Init_Setka(SS)

        read(1) SS%par_n_points                         ! Всего точек в сетке (считается при инициализации сетки)

        read(1) SS%gl_yzel

        read(1) SS%gl_RAY_A
        read(1) SS%gl_RAY_B
        read(1) SS%gl_RAY_C
        read(1) SS%gl_RAY_O
        read(1) SS%gl_RAY_K
        read(1) SS%gl_RAY_D
        read(1) SS%gl_RAY_E

        ! Ячейки
        read(1) SS%gl_Cell_A
        read(1) SS%gl_Cell_B
        read(1) SS%gl_Cell_C

        read(1) SS%gl_all_Cell

        read(1) SS%gl_Cell_neighbour

        read(1) SS%gl_Cell_gran
        read(1) SS%gl_Cell_gran_dist

        read(1) SS%gl_Cell_Centr
        

        read(1) SS%gl_all_Gran
        read(1) SS%gl_Gran_neighbour
        read(1) SS%gl_Gran_normal
        read(1) SS%gl_Gran_length                     
        read(1) SS%gl_Gran_Center                      
        

        
        read(1) SS%gl_Cell_belong
        read(1) SS%gl_Cell_square

        read(1) SS%gl_Cell_type
        read(1) SS%gl_Cell_number

        ! Поверхности выделения
        read(1) SS%gl_HP
        read(1) SS%gl_TS
        read(1) SS%gl_BS

        read(1) SS%gl_Gran_type   

        !! ФИЗИКА
        read(1) SS%n_Hidrogen
        read(1) SS%n_par

        read(1) SS%gd
        read(1) SS%hydrogen

        !if(name == "B0034") then
        read(1) n
        if(n == 1) then
            read(1) n1
            if(n1 /= size(SS%atom_all_source(:, 1, 1))) STOP "ERROR 2751 Geometry 87843y7t8uhihcw4uquiojyre5hgqwc4 "
            read(1) n2
            if(n2 /= size(SS%atom_all_source(1, :, 1))) STOP "ERROR 2752 Geometry 97087867jklpkjertvefervfgbytryhgf5 "
            read(1) n3
            if(n3 /= size(SS%atom_all_source(1, 1, :))) STOP "ERROR 2753 Geometry jytyjtyudreretsadeswcfewrcfrfeergege "
            read(1) SS%atom_all_source

            read(1) n1
            if(n1 /= size(SS%atom_source(:, 1))) STOP "ERROR 2754 Geometry ,imntbrtvecwxqervtyn7i6ue6y "
            read(1) n2
            if(n2 /= size(SS%atom_source(1, :))) STOP "ERROR 2755 Geometry cew3tyukb6vcx3e3eftegr6437iol9ek6i7u "
            read(1) SS%atom_source
        end if
        !end if

        call Geo_Culc_normal(SS, 1) 
        call Geo_Culc_length_area(SS, 1)
        call Culc_Cell_Centr(SS, 2)
        call Geo_Culc_normal(SS, 2) 
        call Geo_Culc_length_area(SS, 2)
        call Geo_culc_TVD_sosed(SS)
        call Geo_Culc_zone(SS)

        call Proverka_grans_sosed(SS)

        close(1)

        do i = 1, size(SS%hydrogen(1, 1, :, 1))
            do j = 1, size(SS%hydrogen(1, :, 1, 1))
                do k = 1, 2
                    if(SS%hydrogen(1, j, i, k) <= 0.0) SS%hydrogen(1, j, i, k) = 0.00000001
                    if(SS%hydrogen(2, j, i, k) <= 0.0) SS%hydrogen(2, j, i, k) = 0.00000001
                end do
            end do
        end do


        do i = 1, size(SS%atom_source(1, :))
            if(ieee_is_nan(SS%atom_source(1, i))) then
                SS%atom_source(1, i) = 0.0
            end if

            if(ieee_is_nan(SS%atom_source(2, i))) then
                SS%atom_source(2, i) = 0.0
            end if

            if(ieee_is_nan(SS%atom_source(3, i))) then
                SS%atom_source(3, i) = 0.0
            end if
        end do

    end subroutine Read_setka_bin

end module GEOMETRY