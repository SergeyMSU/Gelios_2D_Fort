module GEOMETRY
    use STORAGE 
    implicit none 

    contains

    subroutine Init_Setka(SS)        ! Выделение памяти под все массивы сетки
        ! Предполагается, что все параметры сетки определены
        TYPE (Setka), intent(in out) :: SS
        integer(4) :: N1, n2, n

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

        allocate(SS%gl_Cell_gran(4, N1))
        allocate(SS%gl_Cell_belong(3, 4, N1))
        allocate(SS%gl_Cell_square(N1, 2))

        ! Посчитаем число узлов в сетке
        SS%par_n_points = SS%par_n_END * (SS%par_m_A + SS%par_m_BC) + SS%par_m_K * (SS%par_n_TS + SS%par_m_O) + (SS%par_n_END - SS%par_n_TS + 1) * SS%par_m_O - &
                (SS%par_m_A + SS%par_m_BC + SS%par_m_K - 1)  ! Всего точек в сетке
        
        allocate(SS%gl_yzel(2, SS%par_n_points, 2))

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
        allocate(SS%gl_Gran_normal(2, n, 2))
        allocate(SS%gl_Gran_length(n, 2))

        allocate(SS%gl_Contact( (SS%par_m_O + SS%par_m_A + SS%par_m_BC - 1)  ))   ! Выделяем память под контакт
        allocate(SS%gl_TS( (SS%par_m_A + SS%par_m_BC + SS%par_m_K - 1) ))   ! Выделяем память под TS
        allocate(SS%gl_BS( (SS%par_m_A - 1) ))   ! Выделяем память под BS

        SS%gl_Contact = 0
        SS%gl_TS = 0
        SS%gl_BS = 0
        SS%gl_Gran_type = 0

        SS%gl_all_Gran = 0
        SS%gl_Gran_neighbour = 0
        SS%gl_Cell_gran = 0
        SS%gl_Cell_belong = 0.0
        SS%gl_Gran_normal = 0.0
        SS%gl_Gran_length = 0.0
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


    end subroutine Init_Setka

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
        real(8) :: the, r

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
		else
			r = R_BS + (SS%par_R_END - R_BS) * (DBLE(i- SS%par_n_BS)/(SS%par_n_END - SS%par_n_BS))**(SS%par_kk2 * (0.55 + 0.45 * cos(the)) )
		end if

        SS%gl_yzel(1, node, step) = r * cos(the)
        SS%gl_yzel(2, node, step) = r * sin(the)
    end subroutine Set_Ray_A

    subroutine Set_Ray_B(SS, i, j, R_TS, R_HP, step)
        ! step - 1 или 2  показывает какой массив координат мы меняем
        ! R_TS - это расстояние от TS до центра
        !! R_HP - это расстояние по второму лучу от TS до HP
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
        real(8) :: the, r, rr, xx, R_BS

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
		else
			r = R_BS + (DBLE(i - (SS%par_n_BS - SS%par_n_HP + 1))/(N1 - (SS%par_n_BS - SS%par_n_HP + 1) ))**(0.55 * SS%par_kk2) * (SS%par_R_END - R_BS)
		end if

        SS%gl_yzel(1, node, step) = xx
        SS%gl_yzel(2, node, step) = r

        return
    end subroutine Set_Ray_C

    subroutine Set_Ray_O(SS, i, j, R_HP, step)
        ! step - 1 или 2  показывает какой массив координат мы меняем
        ! R_HP - высота
        ! Распределение x - координат считается автоматически от x - координаты крайней точки B на гелиопаузе
        TYPE (Setka), intent(in out) :: SS
        integer(4), INTENT(IN) :: i, j, step
        real(8), INTENT(IN) :: R_HP

        integer(4) :: node, N2, N1
        real(8) :: the, r, rr, xx, R_BS, x

        node = SS%gl_RAY_O(i, j)
        N2 = size(SS%gl_RAY_O(1, :))
        N1 = size(SS%gl_RAY_O(:, 1))

        R_BS = SS%gl_yzel(2, SS%gl_RAY_A(SS%par_n_BS, SS%par_m_A), step)
        xx = SS%gl_yzel(1, SS%gl_RAY_B(SS%par_n_HP, SS%par_m_BC), step)              ! Отталкиваемся от x - координаты крайней точки B на гелиопаузе в этой плоскости (k)
        

        x = xx - (DBLE(j)/N2)**SS%par_kk31 * (xx - SS%par_R_LEFT)

        if (i <= SS%par_n_BS - SS%par_n_HP + 1) then
			r = R_HP + (R_BS - R_HP) * (DBLE(i - 1)/(SS%par_n_BS - SS%par_n_HP))**1.6
		else
			r = R_BS + (DBLE(i - (SS%par_n_BS - SS%par_n_HP + 1))/(N1 - (SS%par_n_BS - SS%par_n_HP + 1) ))**(0.55 * SS%par_kk2) * (SS%par_R_END - R_BS)
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

    subroutine Geo_Find_Cell(SS, x, y, num)
        ! Поиск номера ячейки по её координатам
        ! num = предположительный изначальный номер (если не знаем пусть будет равен 1)
        TYPE (Setka), intent(in) :: SS
        real(8), intent(in) :: x, y
        integer(4), intent(in out) :: num
        integer(4) :: j, gran, sosed, max_num

        max_num = 0

        loop1:do while(.TRUE.)
            max_num = max_num + 1

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

    subroutine Geo_Culc_length_area(SS, step)
        ! Считаем длины граней и площади ячеек, используя координаты узлов на слое step
        TYPE (Setka), intent(in out) :: SS
        integer(4), INTENT(IN) :: step
        integer(4) :: i, n, y1, y2, y3, y4
        real(8) :: p1(2), p2(2), S

        ! Считаем длины граней
        n = size(SS%gl_Gran_length(:, 1))

        do i = 1, n
            y1 = SS%gl_all_Gran(1, i)
            y2 = SS%gl_all_Gran(2, i)
            p1 = SS%gl_yzel(:, y1, step)
            p2 = SS%gl_yzel(:, y2, step)
            SS%gl_Gran_length(i, step) = norm2(p2 - p1)
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
                STOP "ERROR Geo_Culc_length_area 1178  ioyefdvenmlgp9yrtgb"
            end if
        end do

    end subroutine Geo_Culc_length_area

    subroutine Geo_Culc_normal(SS, step)
        ! Считаем нормали всех ячеек
        TYPE (Setka), intent(in out) :: SS
        integer(4), INTENT(IN) :: step
        integer(4) :: N, i, a1, a2
        real(8) :: r1(2), r2(2), normal(2), c, norm, centr(2), vec(2)

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

            centr = (r1 + r2)/2.0
            vec = centr - r1

            if(DOT_PRODUCT(vec, normal) < 0.0) then
                normal = -normal
            end if

            SS%gl_Gran_normal(:, i, step) = normal

        end do

    end subroutine Geo_Culc_normal

    subroutine Geo_Find_Surface(SS)   ! Заполняет массивы поверхностей, которые выделяются
		use STORAGE
		implicit none
        TYPE (Setka), intent(in out) :: SS
		integer :: j, k, node, num, i

		! HP
		node = 1

		do j = 1, size( SS%gl_Cell_A(SS%par_n_HP - 1, :) )
			SS%gl_Contact(node) = SS%gl_Cell_gran(1, SS%gl_Cell_A(SS%par_n_HP - 1, j))
			SS%gl_Gran_type(SS%gl_Contact(node)) = 2
			node = node + 1
				
				
		end do

		do j = 1, size( SS%gl_Cell_C(SS%par_n_HP - SS%par_n_TS, :) )
			SS%gl_Contact(node) = SS%gl_Cell_gran(1, SS%gl_Cell_C(SS%par_n_HP - SS%par_n_TS, j))
			SS%gl_Gran_type(SS%gl_Contact(node)) = 2
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
        write(1,*) "TITLE = 'HP'  VARIABLES = 'X', 'Y', 'Type' "

        do j = 1, size(SS%gl_Cell_A(1, :))
            do i = 1, size(SS%gl_Cell_A(:, 1))
                node = SS%gl_Cell_A(i, j)
                write(1,*) SS%gl_Cell_Centr(:, node, 1), 1
            end do
        end do

        do j = 1, size(SS%gl_Cell_B(1, :))
            do i = 1, size(SS%gl_Cell_B(:, 1))
                node = SS%gl_Cell_B(i, j)
                write(1,*) SS%gl_Cell_Centr(:, node, 1), 2
            end do
        end do

        do j = 1, size(SS%gl_Cell_C(1, :))
            do i = 1, size(SS%gl_Cell_C(:, 1))
                node = SS%gl_Cell_C(i, j)
                write(1,*) SS%gl_Cell_Centr(:, node, 1), 3
            end do
        end do

        close(1)

    end subroutine Print_Cell_Centr

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

    subroutine Print_Grans(SS)
        ! Печатаем все ячейки (есть дублирование), каждая ячейка печатается отдельно
        TYPE (Setka), intent(in) :: SS
        integer(4) :: N, i, a1, a2, j

        N = size(SS%gl_all_Gran(1, :))

        open(1, file = SS%name // '_Print_Grans.txt')
        write(1,*) "TITLE = 'HP'  VARIABLES = 'X', 'Y'  ZONE T= 'HP', N= ", 2 * N, ", E =  ", N , ", F=FEPOINT, ET=LINESEG"

        do i = 1, N
            a1 = SS%gl_all_Gran(1, i)
            a2 = SS%gl_all_Gran(2, i)
            write(1,*) SS%gl_yzel(:, a1, 1)
            write(1,*) SS%gl_yzel(:, a2, 1)
        end do

        do j = 0, N
            write(1,*) 2 * j + 1, 2 * j + 2
        end do

        close(1)
    end subroutine Print_Grans

    subroutine Geo_Print_Surface(SS)
        ! Печатаем поверхности, которые выделяем
        TYPE (Setka), intent(in) :: SS
        integer(4) :: i, N, a1, a2, gr, j

        N = size(SS%gl_Contact) + size(SS%gl_TS) + size(SS%gl_BS)

        open(1, file = SS%name // '_Print_Surface.txt')
        write(1,*) "TITLE = 'HP'  VARIABLES = 'X', 'Y'  ZONE T= 'HP', N= ", 2 * N, ", E =  ", N , ", F=FEPOINT, ET=LINESEG"

        do i = 1, size(SS%gl_Contact)
            gr = SS%gl_Contact(i)
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
        integer(4) :: i, j, sosed, gran, N

        print*, "Proverka_grans_sosed start"

        N = size(SS%gl_Cell_A) + size(SS%gl_Cell_B) + size(SS%gl_Cell_C)

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

end module GEOMETRY