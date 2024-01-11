
module STORAGE
    implicit none 

    !! Набор общепринятых констант (которые никогда не поменяются)
    real(8), parameter :: par_pi = acos(-1.0_8) 

    !! Модуль хранит всю сетку со всеми параметрами
    TYPE Setka 

        character(len=5) :: name = "00000"
        logical :: init_geo = .False.   ! Инициализирована ли данная сетка (выделена ли память под массивы геометрии)


        ! Набор переменных, определяющих структуру сетки
        integer(4) :: par_m_A = 7! 30      ! Количество лучей A в плоскости
        integer(4) :: par_m_BC = 6! 18      ! Количество лучей B/C в плоскости
        integer(4) :: par_m_O = 7! 17      ! Количество лучей O в плоскости
        integer(4) :: par_m_K = 4! 7      ! Количество лучей K в плоскости
        real(8) :: par_triple_point = 13.0 * par_pi/40.0     ! До какого угла начиная от pi/2 (с положительного x) тройная точка
        real(8) :: par_triple_point_2 = 7.0 * par_pi/40.0     ! Под каким углом выходит луч после тройной точки начиная от pi/2 (с положительного x) 
        
        ! Количество точек по лучам A
        integer(4) :: par_n_TS =  15! 26                    ! Количество точек до TS (TS включается)
        integer(4) :: par_n_HP =  25! 40                 ! Количество точек до HP (HP включается)  всё от 0 считается
        integer(4) :: par_n_BS =  35! 60! 5                 ! Количество точек BS (BS включается)
        integer(4) :: par_n_END = 45! 72! 6                ! Количество точек до конца сетки (конец включается)
        integer(4) :: par_n_IA =  8! 12                   ! Количество точек, которые входят во внутреннюю область
        integer(4) :: par_n_IB =  10! 14                   ! Количество точек, которые входят во внутреннюю область (с зазором)

        ! Набор параметров, задающих размеры сетки
        real(8) :: par_R_character = 35.0         ! Характерный размер в задаче (расстояние до TS на начальном этапе построения сетки)
        real(8) :: par_R0 = 0.197035         ! Характерный размер 1 а.е. (внутренней сферы) Там находится вторая точка на лучах от цетра (первая находится в нуле)
        real(8) :: par_R_END = 300.0         !  
        real(8) :: par_R_LEFT = -240.0 ! -390.0         !  Левая граница
        real(8) :: par_R_inner = 9.0! 5.0_8     ! До какого расстояния внутренняя сфера

        !Набор параметров сгущения
        real(8) :: par_kk1 = 1.7_8     ! Степень сгущения сетки к нулю в области до TS: 1 - линейное, 2 - квадратичное и т.д.
        real(8) :: par_kk2 = 2.0_8     ! Степень сгущения в головной области на бесконечности
        real(8) :: par_kk3 = 1.8_8     ! Степень сгущения в хвосте
        real(8) :: par_kk31 = 1.0_8     ! Степень сгущения в хвосте для точек на контакте (первая точка в О - луче)
        real(8) :: par_kk13 = 1.8_8     ! Степень сгущения точек в головной области во внешнем ударном слое  от 0 до 1
        real(8) :: par_kk131 = 0.1_8
        real(8) :: par_kk132 = 1.5_8
        real(8) :: par_kk14 = 1.0_8     ! Степень сгущения точек в головной области во внутреннем ударном слое  от 0 до 1
        ! (сгущение сразу к TS и HP)  
        real(8) :: par_kk12 = 1.4_8     ! Степень сгущения точек до TS к ударной волне  >= 1
        ! Должно делиться на 4 для удобного вывода результатов в плоскостях

        integer(4) :: par_n_points                         ! Всего точек в сетке (считается при инициализации сетки)

        real(8), allocatable :: gl_yzel(:, :)   ! (2, :) набор координат узлов сетки
        real(8), allocatable :: gl_yzel_2(:, :)   ! (2, :) набор координат узлов сетки на следующем временном шаге

        ! Лучи, на которых распологаются точки сетки
        integer(4), allocatable :: gl_RAY_A(:,:)   ! Набор А-лучей размерности 3 (на этом луче, в этой плоскости)
        integer(4), allocatable :: gl_RAY_B(:,:)   ! Набор B-лучей размерности 3 (на этом луче, в этой плоскости)
        integer(4), allocatable :: gl_RAY_C(:,:)   ! Набор C-лучей размерности 3 (на этом луче, в этой плоскости)
        integer(4), allocatable :: gl_RAY_O(:,:)   ! Набор O-лучей размерности 3 (на этом луче, в этой плоскости)
        integer(4), allocatable :: gl_RAY_K(:,:)   ! Набор K-лучей размерности 3 (на этом луче, в этой плоскости)
        integer(4), allocatable :: gl_RAY_D(:,:)   ! Набор D-лучей размерности 3 (на этом луче, в этой плоскости)
        integer(4), allocatable :: gl_RAY_E(:,:)   ! Набор E-лучей размерности 3 (на этом луче, в этой плоскости)

        ! Ячейки
        integer(4), allocatable :: gl_Cell_A(:,:)   ! Набор A-ечеек размерности 3 (на этом луче, в этой плоскости)
        integer(4), allocatable :: gl_Cell_B(:,:)   ! Набор B-ечеек размерности 3 (на этом луче, в этой плоскости)
        integer(4), allocatable :: gl_Cell_C(:,:)   ! Набор C-ечеек размерности 3 (на этом луче, в этой плоскости)

        integer(4), allocatable :: gl_all_Cell(:,:)   ! Весь набор ячеек (4, :) - первая координата массива - это набор узлов ячейки

        integer(4), allocatable :: gl_Cell_neighbour(:,:)   ! (4, :) Набор из 4 соседей для каждой ячейки 
        ! -1   ! Граница (набегающий поток)
        ! -3   ! Граница верхний цилиндр
        ! -4   ! Ось симметрии
	END TYPE Setka
	
	!! Набор глобальных переменных 
    TYPE (Setka):: gl_S1


    contains 

    subroutine Init_Setka(SS)        ! Выделение памяти под все массивы в сетки
        ! Предполагается, что все параметры сетки определены
        TYPE (Setka), intent(in out) :: SS

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
        
        allocate(SS%gl_all_Cell(4, size(SS%gl_Cell_A(:,:)) + size(SS%gl_Cell_B(:,:)) + size(SS%gl_Cell_C(:,:)) ) )
        allocate(SS%gl_Cell_neighbour(4, size(SS%gl_Cell_A(:,:)) + size(SS%gl_Cell_B(:,:)) + size(SS%gl_Cell_C(:,:))))

        ! Посчитаем число узлов в сетке
        SS%par_n_points = SS%par_n_END * (SS%par_m_A + SS%par_m_BC) + SS%par_m_K * (SS%par_n_TS + SS%par_m_O) + (SS%par_n_END - SS%par_n_TS + 1) * SS%par_m_O - &
                (SS%par_m_A + SS%par_m_BC + SS%par_m_K - 1)  ! Всего точек в сетке
        
        allocate(SS%gl_yzel(2, SS%par_n_points))
        allocate(SS%gl_yzel_2(2, SS%par_n_points))

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
        SS%gl_yzel_2 = 0.0_8


    end subroutine Init_Setka

    subroutine Build_Setka_start(SS)        ! Начальное построение сетки
        TYPE (Setka), intent(in out) :: SS

        integer(4) :: i, j, N1, N2, i1, node, kk2, ni, num1
        real(8) :: r, phi, the, xx, x, y, z, rr, x2, y2, z2, R_HP, the2

        node = 1
        ! Задаём первый узел
        SS%gl_yzel(:, node) = (/ 0.0_8, 0.0_8 /)
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
        
    end subroutine Build_Setka_start

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

        if(step == 1) then
            SS%gl_yzel(1, node) = r * cos(the)
            SS%gl_yzel(2, node) = r * sin(the)
        else
            SS%gl_yzel_2(1, node) = r * cos(the)
            SS%gl_yzel_2(2, node) = r * sin(the)
        end if
    end subroutine Set_Ray_A

    subroutine Set_Ray_B(SS, i, j, R_TS, R_HP, step)
        ! step - 1 или 2  показывает какой массив координат мы меняем
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

            if(step == 1) then
                SS%gl_yzel(1, node) = R_TS * cos(the) + r * cos(the2)
                SS%gl_yzel(2, node) = R_TS * sin(the) + r * sin(the2)
            else
                SS%gl_yzel_2(1, node) = R_TS * cos(the) + r * cos(the2)
                SS%gl_yzel_2(2, node) = R_TS * sin(the) + r * sin(the2)
            end if
            return
        end if

        if(step == 1) then
            SS%gl_yzel(1, node) = r * cos(the)
            SS%gl_yzel(2, node) = r * sin(the)
        else
            SS%gl_yzel_2(1, node) = r * cos(the)
            SS%gl_yzel_2(2, node) = r * sin(the)
        end if

        return
    end subroutine Set_Ray_B

    subroutine Set_Ray_C(SS, i, j, step)
        ! step - 1 или 2  показывает какой массив координат мы меняем
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
        if(step == 1) then
            xx = SS%gl_yzel(1, SS%gl_RAY_B(SS%par_n_HP, j))
            rr = SS%gl_yzel(2, SS%gl_RAY_B(SS%par_n_HP, j))
            R_BS = SS%gl_yzel(2, SS%gl_RAY_A(SS%par_n_BS, SS%par_m_A))
        else
            xx = SS%gl_yzel_2(1, SS%gl_RAY_B(SS%par_n_HP, j))
            rr = SS%gl_yzel_2(2, SS%gl_RAY_B(SS%par_n_HP, j))
            R_BS = SS%gl_yzel_2(2, SS%gl_RAY_A(SS%par_n_BS, SS%par_m_A))
        end if


        if (i <= SS%par_n_BS - SS%par_n_HP + 1) then
			r = rr + (R_BS - rr) * (DBLE(i)/(SS%par_n_BS - SS%par_n_HP + 1))**1.6
		else
			r = R_BS + (DBLE(i - (SS%par_n_BS - SS%par_n_HP + 1))/(N1 - (SS%par_n_BS - SS%par_n_HP + 1) ))**(0.55 * SS%par_kk2) * (SS%par_R_END - R_BS)
		end if

        if(step == 1) then
            SS%gl_yzel(1, node) = xx
            SS%gl_yzel(2, node) = r
        else
            SS%gl_yzel_2(1, node) = xx
            SS%gl_yzel_2(2, node) = r
        end if

        return
    end subroutine Set_Ray_C

    subroutine Set_Ray_O(SS, i, j, R_HP, step)
        ! step - 1 или 2  показывает какой массив координат мы меняем
        TYPE (Setka), intent(in out) :: SS
        integer(4), INTENT(IN) :: i, j, step
        real(8), INTENT(IN) :: R_HP

        integer(4) :: node, N2, N1
        real(8) :: the, r, rr, xx, R_BS, x

        node = SS%gl_RAY_O(i, j)
        N2 = size(SS%gl_RAY_O(1, :))
        N1 = size(SS%gl_RAY_O(:, 1))

        if(step == 1) then
            R_BS = SS%gl_yzel(2, SS%gl_RAY_A(SS%par_n_BS, SS%par_m_A))
            xx = SS%gl_yzel(1, SS%gl_RAY_B(SS%par_n_HP, SS%par_m_BC))              ! Отталкиваемся от x - координаты крайней точки B на гелиопаузе в этой плоскости (k)
        else
            R_BS = SS%gl_yzel_2(2, SS%gl_RAY_A(SS%par_n_BS, SS%par_m_A))
            xx = SS%gl_yzel_2(1, SS%gl_RAY_B(SS%par_n_HP, SS%par_m_BC))
        end if

        x = xx - (DBLE(j)/N2)**SS%par_kk31 * (xx - SS%par_R_LEFT)

        if (i <= SS%par_n_BS - SS%par_n_HP + 1) then
			r = R_HP + (R_BS - R_HP) * (DBLE(i - 1)/(SS%par_n_BS - SS%par_n_HP))**1.6
		else
			r = R_BS + (DBLE(i - (SS%par_n_BS - SS%par_n_HP + 1))/(N1 - (SS%par_n_BS - SS%par_n_HP + 1) ))**(0.55 * SS%par_kk2) * (SS%par_R_END - R_BS)
		end if

        if(step == 1) then
            SS%gl_yzel(1, node) = x
            SS%gl_yzel(2, node) = r
        else
            SS%gl_yzel_2(1, node) = x
            SS%gl_yzel_2(2, node) = r
        end if

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

        if(step == 1) then
            SS%gl_yzel(1, node) = r * cos(the)
            SS%gl_yzel(2, node) = r * sin(the)
        else
            SS%gl_yzel_2(1, node) = r * cos(the)
            SS%gl_yzel_2(2, node) = r * sin(the)
        end if

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

        if(step == 1) then
            if (j < N2) then
                xx = SS%gl_yzel(1, SS%gl_RAY_K(SS%par_n_TS, j))
                y = SS%gl_yzel(2, SS%gl_RAY_K(SS%par_n_TS, j))
            else
                xx = SS%gl_yzel(1, SS%gl_RAY_B(SS%par_n_TS, SS%par_m_BC))
                y = SS%gl_yzel(2, SS%gl_RAY_B(SS%par_n_TS, SS%par_m_BC))
            end if
        else
            if (j < N2) then
                xx = SS%gl_yzel_2(1, SS%gl_RAY_K(SS%par_n_TS, j))
                y = SS%gl_yzel_2(2, SS%gl_RAY_K(SS%par_n_TS, j))
            else
                xx = SS%gl_yzel_2(1, SS%gl_RAY_B(SS%par_n_TS, SS%par_m_BC))
                y = SS%gl_yzel_2(2, SS%gl_RAY_B(SS%par_n_TS, SS%par_m_BC))
            end if
        end if


        x = xx + (DBLE(i - 1)/(N1 - 1))**SS%par_kk3 * (SS%par_R_LEFT - xx)

        if(step == 1) then
            SS%gl_yzel(1, node) = x
            SS%gl_yzel(2, node) = y
        else
            SS%gl_yzel_2(1, node) = x
            SS%gl_yzel_2(2, node) = y
        end if

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

        if(step == 1) then
            x = SS%gl_yzel(1, SS%gl_RAY_E(1, j))
            y = SS%gl_yzel(2, SS%gl_RAY_E(1, j))

            xx = SS%gl_yzel(1, SS%gl_RAY_O(1, j))
            yy = SS%gl_yzel(2, SS%gl_RAY_O(1, j))
        else
            x = SS%gl_yzel_2(1, SS%gl_RAY_E(1, j))
            y = SS%gl_yzel_2(2, SS%gl_RAY_E(1, j))

            xx = SS%gl_yzel_2(1, SS%gl_RAY_O(1, j))
            yy = SS%gl_yzel_2(2, SS%gl_RAY_O(1, j))
        end if

        if(step == 1) then
            SS%gl_yzel(1, node) = x + (xx - x) * sgushenie_2(DBLE(i - 1)/(N1 - 1), SS%par_kk14)
            SS%gl_yzel(2, node) = y + (yy - y) * sgushenie_2(DBLE(i - 1)/(N1 - 1), SS%par_kk14)
        else
            SS%gl_yzel_2(1, node) = x + (xx - x) * sgushenie_2(DBLE(i - 1)/(N1 - 1), SS%par_kk14)
            SS%gl_yzel_2(2, node) = y + (yy - y) * sgushenie_2(DBLE(i - 1)/(N1 - 1), SS%par_kk14)
        end if

        return
    end subroutine Set_Ray_E

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
                    write(1,*) SS%gl_yzel(:, n), 1
                end if
            end do
        end do

        N1 = size(SS%gl_RAY_B(:, 1))
        N2 = size(SS%gl_RAY_B(1, :))
        do j = 1, N2
            do i = 1, N1
                n = SS%gl_RAY_B(i, j)
                if(n /= 0) then
                    write(1,*) SS%gl_yzel(:, n), 2
                end if
            end do
        end do

        N1 = size(SS%gl_RAY_C(:, 1))
        N2 = size(SS%gl_RAY_C(1, :))
        do j = 1, N2
            do i = 1, N1
                n = SS%gl_RAY_C(i, j)
                if(n /= 0) then
                    write(1,*) SS%gl_yzel(:, n), 3
                end if
            end do
        end do

        N1 = size(SS%gl_RAY_O(:, 1))
        N2 = size(SS%gl_RAY_O(1, :))
        do j = 1, N2
            do i = 1, N1
                n = SS%gl_RAY_O(i, j)
                if(n /= 0) then
                    write(1,*) SS%gl_yzel(:, n), 4
                end if
            end do
        end do

        N1 = size(SS%gl_RAY_K(:, 1))
        N2 = size(SS%gl_RAY_K(1, :))
        do j = 1, N2
            do i = 1, N1
                n = SS%gl_RAY_K(i, j)
                if(n /= 0) then
                    write(1,*) SS%gl_yzel(:, n), 5
                end if
            end do
        end do

        N1 = size(SS%gl_RAY_D(:, 1))
        N2 = size(SS%gl_RAY_D(1, :))
        do j = 1, N2
            do i = 1, N1
                n = SS%gl_RAY_D(i, j)
                if(n /= 0) then
                    write(1,*) SS%gl_yzel(:, n), 6
                end if
            end do
        end do

        N1 = size(SS%gl_RAY_E(:, 1))
        N2 = size(SS%gl_RAY_E(1, :))
        do j = 1, N2
            do i = 1, N1
                n = SS%gl_RAY_E(i, j)
                if(n /= 0) then
                    write(1,*) SS%gl_yzel(:, n), 7
                end if
            end do
        end do

        close(1)

    end subroutine Print_Point_from_Rays

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

    end subroutine Print_Cell

    subroutine Print_Connect(SS)
        ! Печатаем все связи у ячеек
        TYPE (Setka), intent(in) :: SS

    end subroutine Print_Connect

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

end module STORAGE