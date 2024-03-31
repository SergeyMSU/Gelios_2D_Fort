
module STORAGE
    USE OMP_LIB
    implicit none 

    !! Набор общепринятых констант (которые никогда не поменяются)
    real(8), parameter :: par_pi = acos(-1.0_8) 
    real(8), parameter :: par_sqrtpi = sqrt(par_pi)

    integer(4), parameter :: par_n_zone = 6! 7  !  Количество радиусов (но есть ещё внешняя зона)
	integer(4), parameter :: par_m_zone = 7! 6  !  Количество лучей по углу (от 0 до 180)
    integer(4), parameter :: par_n_potok = 32! 32! 24! 32  ! Число потоков (у каждого потока свой стек)
    integer(4), parameter :: par_n_claster = 1  ! Число компьютеров (для MPI)
    integer(4), parameter :: par_n_parallel = 20! 20  ! Для распараллеливания цикла (т.е. каждый поток будет в среднем обрабатывать такое число итераций
    integer(4), parameter :: par_stek = 1000  ! Глубина стека (заранее выделяется память под него)
    integer(4), parameter :: par_n_sort = 4!4!6  !  4 Количество сортов атомов

    logical, parameter :: MK_is_NaN = .False.    ! Нужны ли проверки на nan
	logical, parameter :: MK_Mu_stat = .True.    ! Нужно ли накапливать веса для статистики и весовых каэффициентов
	logical, parameter :: MK_photoionization = .True.    ! Нужна ли фотоионизация
	logical, parameter :: MK_el_impact = .False.    ! Нужен ли электронный удар
	logical, parameter :: MK_statistik_file = .True.    ! Какой файл со статистикой использовать?  true - это первый файл

	logical, parameter :: par_Hydro = .True.    ! Нужен ли водород? Или чисто газовую-динамику считаем?
	logical, parameter :: par_TVD_linear_HP = .True.! .False. !.True.    ! Нужно ли сносить на контакт линейно только с одной стороны
	logical, parameter :: par_move_setka = .True.    ! Нужно ли сносить на контакт линейно только с одной стороны
    real(8), parameter :: par_Rmax = 220.0 !300.0! 220.0  !  Радиус сферы, с которой запускаем частицы


    ! Число частиц у каждого потока!
	! Число должно быть кратно par_n_parallel
	integer(4), parameter :: MK_k_multiply = 12 * 6! 12 * 2!11 * 8!12 * 10! 12 * 7!14 * 8!6 * 3! * 6 * 9! * 6 * 8!6 * 6 * 2  !   ! 6 = 20 минут счёта (с пикапами 30 минут)
    ! 9 сейчас с пикапами
    ! 12 (14) - это 1 час с пикапами
    ! 18 - это 1 час без пикапов
	integer(4), parameter :: MK_k_mul1 = 6 * MK_k_multiply! 6
	integer(4), parameter :: MK_k_mul2 = 1 * MK_k_multiply! 
	integer(4), parameter :: MK_N1 = MK_k_mul1 * 60/par_n_parallel   ! 60 Число исходных частиц первого типа (с полусферы)
	integer(4), parameter :: MK_N2 = MK_k_mul1 * 20/par_n_parallel   ! 20
	integer(4), parameter :: MK_N3 = MK_k_mul2 * 20/par_n_parallel   ! 20
	integer(4), parameter :: MK_N4 = MK_k_mul1 * 20/par_n_parallel   ! 20

    !! Модуль хранит всю сетку со всеми параметрами
    TYPE Setka 

        character(len=5) :: name = "00000"
        logical :: init_geo = .False.   ! Инициализирована ли данная сетка (выделена ли память под массивы геометрии)


        ! Набор переменных, определяющих структуру сетки
        integer(4) :: par_m_A = 20! 30      ! Количество лучей A в плоскости
        integer(4) :: par_m_BC = 10! 18      ! Количество лучей B/C в плоскости
        integer(4) :: par_m_O = 10! 17      ! Количество лучей O в плоскости
        integer(4) :: par_m_K = 8! 7      ! Количество лучей K в плоскости
        real(8) :: par_triple_point = 13.0 * par_pi/40.0     ! До какого угла начиная от pi/2 (с положительного x) тройная точка
        real(8) :: par_triple_point_2 = 7.0 * par_pi/40.0     ! Под каким углом выходит луч после тройной точки начиная от pi/2 (с положительного x) 
        
        ! Количество точек по лучам A
        integer(4) :: par_n_TS =  33! 26                    ! Количество точек до TS (TS включается)
        integer(4) :: par_n_HP =  63! 40                 ! Количество точек до HP (HP включается)  всё от 0 считается
        integer(4) :: par_n_BS =  89! 60! 5                 ! Количество точек BS (BS включается)
        integer(4) :: par_n_END = 98! 72! 6                ! Количество точек до конца сетки (конец включается)
        integer(4) :: par_n_IA =  20! 12                   ! Количество точек, которые входят во внутреннюю область
        integer(4) :: par_n_IB =  22! 14                   ! Количество точек, которые входят во внутреннюю область (с зазором)

        ! Набор параметров, задающих размеры сетки
        real(8) :: par_R_character = 35.0         ! Характерный размер в задаче (расстояние до TS на начальном этапе построения сетки)
        real(8) :: par_R0 = 0.198956         ! Характерный размер 1 а.е. (внутренней сферы) Там находится вторая точка на лучах от цетра (первая находится в нуле)
        real(8) :: par_R_END = 300.0         !  
        real(8) :: par_R_LEFT = -240.0 ! -390.0         !  Левая граница
        real(8) :: par_R_inner = 9.0! 5.0_8     ! До какого расстояния внутренняя сфера

        !! Физические параметры ---------------------------------------
        real(8) :: par_a_2 = 0.130735_8        ! Параметр в сечении перезарядки  !! ЗАДАН ТАКЖЕ ГЛОБАЛЬНО, надо переделать
        real(8) :: par_ggg = 5.0/3.0                 ! До какого расстояния внутренняя сфера
        real(8) :: par_Velosity_inf = -2.54385_8 !-2.54278_8
        real(8) :: par_n_H_LISM = 3.0_8
        real(8) :: par_Kn = 50.3858   !0.4326569808         ! в перезарядке
        real(8) :: par_nu_ph = 12.0969 
        real(8) :: par_E_ph = 0.10878
        real(8) :: par_chi = 41.0391
        real(8) :: par_rho_E = 1.0                 !? Сейчас эта переменная нигде не используется (но сохраняется в файл)
        real(8) :: par_Max_e = 10.0!! 5.91662
        real(8) :: par_poglosh = 0.618589! 0.389274        !! На что обезразмериваем поглощение
        real(8) :: par_rho_LISM = 1.60063        !! Коэффициент увеличения плотности из-за гелия на бесконечности и т.д.
        real(8) :: par_p_LISM = 1.15        !! Коэффициент увеличения давления из-за гелия на бесконечности и т.д.
        real(8) :: par_rho_He_Lism = 0.6_8         !! Какая часть от общей плотности - это гелий
        real(8) :: par_rho_He_E = 0.155327_8       

        ! Для электронного удара
        real(8) :: nu_e_impact = 0.214207_8! 3.86624! 0.15465_8       
        real(8) :: lambda_e = 12.1458_8       

        !! -------------------------------------------------------------

        !! Набор параметров сгущения
        real(8) :: par_kk1 = 2.0_8     ! Степень сгущения сетки к нулю в области до TS: 1 - линейное, 2 - квадратичное и т.д.
        real(8) :: par_kk2 = 2.08_8 !1.7_8     ! Степень сгущения в головной области на бесконечности
        real(8) :: par_kk3 = 1.8_8     ! Степень сгущения в хвосте
        real(8) :: par_kk31 = 1.0_8     ! Степень сгущения в хвосте для точек на контакте (первая точка в О - луче)
        real(8) :: par_kk13 = 1.8_8     ! Степень сгущения точек в головной области во внешнем ударном слое  от 0 до 1
        real(8) :: par_kk131 = 0.1_8
        real(8) :: par_kk132 = 1.5_8
        real(8) :: par_kk14 = 1.0_8     ! Степень сгущения точек в головной области во внутреннем ударном слое  от 0 до 1
        ! (сгущение сразу к TS и HP)  
        real(8) :: par_kk12 = 1.0_8     ! Степень сгущения точек до TS к ударной волне  >= 1
        real(8) :: par_kk113 = 1.6_8     ! 1.6
        ! Должно делиться на 4 для удобного вывода результатов в плоскостях

        !! Поверхностное натяжение
        real(8) :: par_nat_TS = 0.02_8 ! 0.003_8   ! Поверхностное натяжение
        real(8) :: par_nat_HP = 0.005_8   ! Поверхностное натяжение
        real(8) :: par_nat_BS = 0.003_8   ! 0.002 Поверхностное натяжение

        !! Скорость движения поверхностей
        real(8) :: par_koeff_TS = 0.002_8 
        real(8) :: par_koeff_HP = 0.1_8   
        real(8) :: par_koeff_BS = 0.1_8   

        !! Параметры для Монте-Карло ----------------------------------------------------------

        integer(4), allocatable :: sensor(:, :, :)  !(3, 2, : par_n_potok - число потоков)  ! датчики случайных чисел 
	    ! Каждому потоку по два датчика

        integer(4), allocatable :: stek(:)   ! (: число потоков) Переменная чтения и записи в стек
	    ! Где стоит переменная, там что-то лежит, чтобы записать, нужно увеличить значение на 1

        real(8), allocatable :: M_K_particle(:, :, :)   ! Частицы (8, par_stek, число потоков)
        ! (три координаты, три скорости, вес, радиус перегелия)
        integer(4), allocatable :: M_K_particle_2(:, :, :)  ! Частицы (5, par_stek, число потоков)
        ! (в какой ячейке частица, сорт, зона назначения по r, зона назначения по углу, в какой ячейке в интерполяционной сетке)
        logical(4), allocatable :: M_K_particle_3(:, :, :, :)  ! Частицы (par_n_zone + 1, par_m_zone + 1, par_stek, число потоков)
        ! Массив для набора статистики весов по зонам 

        integer(4) :: par_n_moment = 13 !9  !  Сколько различных моментов считаем (длинна массива)
        real(8), allocatable :: M_K_Moment(:, :, :, :)  ! (19, par_n_sort, :, par_n_potok) То, что накапливаем в ячейках (по каждому сорту отдельно)
        !(rho, u, v, T, In, Iu, Iv, IT, Huu, Huv, Hvv, Huuu, Hvvv)
        !(1  , 2, 3, 4, 5,  6,  7,  8,  9,   10,  11,  12,   13)

        real(8) :: MK_R_zone(par_n_zone)   ! Радиусы зон
        real(8) :: MK_al_zone(par_m_zone)   ! Лучи зон
        real(8) :: MK_SINKR(par_m_zone + 1)   ! Критические синусы для каждой зоны по углу
        real(8), allocatable :: MK_Mu(:, :, :)   ! Веса зон (par_n_zone + 1, par_m_zone + 1, сортов par_n_sort)
        real(8), allocatable :: MK_Mu_statistic(:, :, :)   ! Веса зон (par_n_zone + 1, par_m_zone + 1, сортов par_n_sort)
        ! Для накапливания весов зон

        real(8) :: MK_gam_zone(par_n_zone)   ! Параметр гамма для зон
        real(8) :: MK_A0_, MK_A1_   ! Параметры для начального запуска

        real(8) :: sqv_1, sqv_2, sqv_3, sqv_4, sqv   ! Потоки частиц через 
        real(8) :: MK_mu1, MK_mu2, MK_mu3, MK_mu4
        integer(4) :: MK_N                  ! Сколько всего частиц запущено (сумма по всем потокам)

        real(8) :: MK_Mu_mult = 100.0_8  ! На что домножаем веса для избежания потери точности

        real(8) :: par_Rleft   ! Левая стенка для Монте-Карло (она правее сеточной)
        real(8) :: par_Rup     ! Верзняя стенка для Монте-Карло (она ниже чем у сетки)

        !! -------------------------------------------------------------------------------------
        

        integer(4) :: par_n_points                         ! Всего точек в сетке (считается при инициализации сетки)

        real(8), allocatable :: gl_yzel(:, :, :)   ! (2, :, 2) набор координат узлов сетки
        real(8), allocatable :: gl_yzel_Vel(:, :)   ! (2, :) Скорость движения узлов сетки  !TODO NO-SAVE
        ! Этот массив не надо сохранять при сохранении сетки в файл - это промежуточный рабочий массив
        integer(4), allocatable :: gl_Point_num(:)   ! Сколько раз скорость записана в узел  !TODO NO-SAVE

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
        !? Гарантируется ли расположение точек в ячейке по кругу?

        integer(4), allocatable :: gl_Cell_neighbour(:,:)   ! (4, :) Набор из 4 соседей для каждой ячейки  !! соседи согласованы с соседями для граней!
        ! -1   ! Граница (набегающий поток)
        ! -2   ! Выходная граница
        ! -3   ! Граница верхний цилиндр
        ! -4   ! Ось симметрии

        integer(4), allocatable :: gl_all_Cell_zone(:)   ! зона ячейки
        ! 1, 2, 3, 4

        integer(4), allocatable :: gl_Cell_gran(:,:)        ! (4, :) Набор из 4 граней для каждой ячейки (если номер = 0, то грани нет в этом направлении)
        ! 0 - нет грани (такие ячейки есть, например возде нуля)
        !! Отрицательных номеров у грани нет!
        real(8), allocatable :: gl_Cell_gran_dist(:, :, :)      ! (4, :, 2) Расстояния от ценра ячеек до центров граней

        real(8), allocatable :: gl_Cell_Centr(:, :, :)   ! (2, : число ячеек, 2) набор координат центров ячеек
        real(8), allocatable :: gl_Cell_alpha_center(:)   ! (число ячеек) полярный угол центра ячейки  !? Не сохранять (всегда можно посчитать)
        

        integer(4), allocatable :: gl_all_Gran(:,:)       ! Все грани (2,:) имеют по 2 узла
        integer(4), allocatable :: gl_Gran_neighbour(:,:) ! Соседи каждой грани (2,:) имеют по 2 соседа, нормаль ведёт от первого ко второму
        real(8), allocatable :: gl_Gran_normal(:, :, :)       ! (2, :, 2) Нормаль грани     !! Нормаль идёт от первой ячейки ко второй обязательно!                   
        real(8), allocatable :: gl_Gran_length(:, :)       ! (:, 2) Длина грани                       
        real(8), allocatable :: gl_Gran_Center(:, :, :)       ! (2, :, 2) Центр грани       

        real(8), allocatable :: gl_Gran_POTOK(:)       ! ПОТОК ГРАНИ УБРАТЬ - проверял что потоки через грани одинаковые с обеих сторон        
        
        integer(4), allocatable :: gl_Gran_neighbour_TVD(:,:) ! TVD-Соседи каждой грани (2,:) имеют по 2 соседа
        ! 0 - значит соседа нет
        ! При этом первый TVD-сосед - это сосед первого обычного соседа. Т.е. нормаль грани тоже ведёт от первого ко второму TVD-соседу

        
        real(8), allocatable :: gl_Cell_belong(:,:,:)      ! (3, 4, :)  ! Определяются для координат на первом слое по времени
        ! для каждой грани коэффициенты A, B, C   Ax + By + C = 0  (если больше 0, то точка вне ячейки! Т.е. нормаль внешняя)
        real(8), allocatable :: gl_Cell_square(:,:)      ! (:, 2)

        character, allocatable :: gl_Cell_type(:)           ! Тип каждой ячейки А, Б, С
        integer(4), allocatable :: gl_Cell_number(:, :)     ! (2, :) номер каждой ячейки внутри своего типа

        ! Поверхности выделения
        integer(4), allocatable :: gl_HP(:)            ! Контакт - номеры граней, которые образуют контактную поверхность
        ! Контакт состоит из номеров граней
        integer(4), allocatable :: gl_TS(:)
        integer(4), allocatable :: gl_BS(:)
        !? ПРОВЕРЕНО, чтобы нормали у граней-поверхностей были ориентированы "наружу"
        !? ПРОВЕРЕНО, что у граней первый и второй узел правильно расположены (с права на лево)
        !! Проверить что сами грани в правильном порядке

        integer(4), allocatable :: gl_Gran_type(:)      ! Показывает тип грани 
        ! (0 - обычная, 1 - TS, 2 - HP, 3 - BS)     

        integer(4), allocatable :: gl_Gran_shem(:)      ! Показывает схему для грани
        ! 0 - Lax
        ! 1 - HLL
        ! 2 - HLLC
        ! 3 - Godunov

        !! ФИЗИКА
        integer(4) :: n_Hidrogen = par_n_sort  ! Число сортов атомов водорода
        integer(4) :: n_atom_source = 7  ! Число источников атомов
        integer(4) :: n_par = 5! 6  !! Число физических параметров в задаче   5 - если нет гелия
        ! 5 газодинамических
        ! 4 * n_Hidrogen - Водород

        real(8), allocatable :: gd(:, :, :)  ! (n_par, :, 2 временной слой)
        ! (rho p u v Q He)
        !   1  2 3 4 5 6
        real(8), allocatable :: hydrogen(:, :, :, :)  ! (5, n_Hidrogen, :, 2 временной слой)
        ! rho p u v T

        real(8), allocatable :: atom_all_source(:, :, :)  ! (4, n_Hidrogen, : число ячеек)
        ! (In, Iu, Iv, IT)

        real(8), allocatable :: atom_source(:, :)  ! (n_atom_source, : число ячеек)
        ! (k_u, k_v, k_T, In, Iu, Iv, IT)
        ! (1     2    3    4   5  6    7)

        LOGICAL :: pogl_ = .True.  ! Считаем ли поглощение
        real(8) :: pogl_v_min = -15.0
        real(8) :: pogl_v_max = 15.0
        integer(4) :: pogl_iter = 300
        real(8) :: pogl_ddd
        real(8), allocatable :: pogloshenie(:, :, :)   !  (n_Hidrogen, рабиений по скорости, ячеек)


        !! PUI 
        ! PUI - тяжёловесный блок, поэтому память выделяется не в "Geometry", а в "PUI" только при необходимости работы с пикапами
        LOGICAL :: culc_pui = .True.   ! Считаем ли PUI ? 
        integer :: pui_nW = 60      ! 50   !TODO СОХРАНЯЕТСЯ
        real(8) :: pui_wR = 150.0    ! 150.0  !TODO СОХРАНЯЕТСЯ
        integer :: pui_size            ! Сколько ячеек содержат pui
        integer :: pui_n_par = 3            ! Сколько ячеек содержат pui    !TODO СОХРАНЯЕТСЯ
        real(8), allocatable :: f_pui(:, :)            ! (pui_nW, : pui_size) !TODO СОХРАНЯЕТСЯ
        integer, allocatable :: f_pui_num(:)           ! По номеру в массиве пуи, определяем номер ячейки в сетке  !TODO СОХРАНЯЕТСЯ
        integer, allocatable :: f_pui_num2(:)		   ! По номеру ячейки в сетке, определяем номер в массиве PUI (если он есть - если область внутренняя) !TODO СОХРАНЯЕТСЯ
        ! 0 - ячейка не содержит pui
        
        real(8), allocatable :: par_pui(:, :)               ! (3 pui_n_par, : pui_size)	   Параметры пикапов !TODO СОХРАНЯЕТСЯ
        ! (n_pui, T_pui, p_pui)
        real(8), allocatable :: pui_Sm(:, :)           ! (pui_nW, : pui_size)  !TODO СОХРАНЯЕТСЯ
	    real(8), allocatable :: pui_Sp(:, :)           ! (pui_nW, : pui_size)  !TODO СОХРАНЯЕТСЯ

        !? функция h0(U_H) - см. документацию PUI  !TODO ДАЛЕЕ НИЧЕГО НЕ СОХРАНЯЕТСЯ
        integer :: pui_h0_n = 1000      
        real(8) :: pui_h0_wc = 100.0    
        real(8), allocatable :: h0_pui(:)   ! Для функции отказов при розыгрыше - её максимум для нормировки (см. документацию)
        

        !? Розыгрышь пикапов (заранее считаем функцию розыгрыша)
        integer :: pui_F_n = 300      ! На сколько частей мы разбиваем первообразную для розыгрыша PUI 
        ! т.е. мы будем разыгрывать ksi от 0 до 1 и брать значения скорости при этой ksi. 
        ! интеграл посчитан для dksi = 1.0/pui_F_n а в значениях между придётся линейно интерполировать
        real(8), allocatable :: F_integr_pui(:, :)           ! (pui_F_n, :) первообразная для розыгрыша
        real(8), allocatable :: nu_integr_pui(:, :)           ! (pui_F_n, :)   частота перезарядки
        real(8), allocatable :: Mz_integr_pui(:, :)           ! (pui_F_n, :)   источник импульса
        real(8), allocatable :: E_integr_pui(:, :)           ! (pui_F_n, :)	   источник энергии
        integer (kind=omp_lock_kind), allocatable :: pui_lock(:)  ! Для openMP

    END TYPE Setka

    TYPE Inter_Setka  ! Сетка для интерполяции

        logical :: init = .False.   ! Инициализирована ли данная сетка (выделена ли память под массивы геометрии)
        real(8), allocatable :: gl_yzel(:, :)   ! (2, :) набор координат узлов сетки

        !? Наблюдательный факт
        ! Почти все ячейки состояит из 4 граней, но есть несколько треугольников
        ! В этом случае четвёрный узел равен первому (для всех кроме одного)
        ! В единсвенном варианте четвёртый равен третьему

        ! Ячейки
        !! Здесь А и B ячейки пересекаются (как раз на стыке A и B)
        integer(4), allocatable :: gl_Cell_A(:,:)   ! Набор A-ечеек размерности 3 (на этом луче, в этой плоскости)
        integer(4), allocatable :: gl_Cell_B(:,:)   ! Набор B-ечеек размерности 3 (на этом луче, в этой плоскости)
        integer(4), allocatable :: gl_Cell_C(:,:)   ! Набор C-ечеек размерности 3 (на этом луче, в этой плоскости)

        integer(4), allocatable :: gl_all_Cell(:,:)   ! Весь набор ячеек (4, :) - первая координата массива - это набор узлов ячейки
        real(8), allocatable :: gl_Cell_center(:,:)   ! (2, :) 

        logical, allocatable :: gl_all_triangle(:)
        ! .True. - если ячейка является треугольником и .False. в противном случае

        integer(4), allocatable :: gl_Cell_neighbour(:,:)   ! (4, :) Набор из 4 соседей для каждой ячейки 
        real(8), allocatable :: gl_Cell_Belong(:, :, :)       ! Все грани (3, 4, :) имеют по A B C, 4 грани в ячейке, 
        ! Первая грань - это 1 и 2 узел, вторая - 2 и 3 и т.д.
        ! Ax + By + C = 0  (если > 0 то за пределами ячейки)

        real(4), allocatable :: gl_Cell_interpol_matrix(:, :, :)   ! (4, 4, :) Набор из 4 соседей для каждой ячейки
        ! Ax + By + Cxy + D
        ! Ax + By + C
        ! 
        ! интерполяционная матрица для каждой ячейки
        ! Эта интерполяция не работает для кривых областей! Точнее работает не правильно (интерполируемое значения бывает больше, чем значения в узлах)

        !! ФИЗИКА
        integer(4) :: n_Hidrogen = 4  ! Число сортов атомов водорода
        integer(4) :: n_par = 5  !! Число физических параметров в задаче
        ! 5 газодинамических
        ! 4 * n_Hidrogen - Водород

        real(8), allocatable :: gd(:, :)  ! (n_par, :)
        ! (rho p u v Q He)
        real(8), allocatable :: hydrogen(:, :, :)  ! (5, n_Hidrogen, :)
        ! rho p u v T

        real(8), allocatable :: atom_all_source(:, :, :)  ! (4, n_Hidrogen, : число ячеек)
        ! (In, Iu, Iv, IT)

        real(8), allocatable :: atom_source(:, :)  ! (7, : число ячеек)
        ! (k_u, k_v, k_T, In, Iu, Iv, IT)

    END TYPE Inter_Setka

    TYPE Surfaces
        ! Модуль хранения поверхностей для движения сетки к этим поверхностям
        ! Такое приём нужен, например, для перестройки текущей сетки (так как добавлять ячейки нельзя, но можно построить другую сетку)
        ! и подвинуть её поверхности
        logical :: init = .False.

        real(8), allocatable :: TS(:, :)  ! (2, :) угол, радиус
        real(8), allocatable :: HP(:, :)  ! (4, :) угол, радиус, x, y
        real(8), allocatable :: BS(:, :)  ! (4, :) угол, радиус, x, y

    END TYPE Surfaces
	
	!! Набор глобальных переменных 
    TYPE (Setka):: gl_S1
    TYPE (Setka):: gl_S3
    TYPE (Setka):: gl_S4

    TYPE (Inter_Setka):: gl_S2, gl_I1

    TYPE (Surfaces):: gl_surf1


    contains 

    real(8) pure function MK_sigma(SS, x)
	    TYPE (Setka), intent(in) :: SS
        real(8), intent (in) :: x
        MK_sigma = (1.0 - SS%par_a_2 * log(x))**2
    end function MK_sigma

    real(8) pure function MK_sigma2(SS, x, y)
	    TYPE (Setka), intent(in) :: SS
        real(8), intent (in) :: x, y
        MK_sigma2 = (1.0 - SS%par_a_2 * log(x * y))**2
    end function MK_sigma2

    subroutine PUI_SET(SS)
        TYPE (Setka), intent(in out) :: SS
        integer(4) :: n, i, j

        if( SS%culc_pui == .False.) then
            print*, "ERROR 9uy87tyg83h8ou9t4358y30q9pug9"
            pause
            STOP
        end if

        ! Проверяем выделена ли память под массивы PUI 
        if(ALLOCATED(SS%f_pui) == .False.) then 
            ! Посчитаем, сколько ячеек будут содержать pui
            n = 0
            do i = 1, size(SS%gl_all_Cell_zone(:))
                j = SS%gl_all_Cell_zone(i)
                if(j <= 2) n = n + 1
            end do
            
            SS%pui_size = n
            ALLOCATE(SS%f_pui(SS%pui_nW, n))
            ALLOCATE(SS%f_pui_num(n))
            ALLOCATE(SS%f_pui_num2, mold = SS%gl_all_Cell_zone)
            ALLOCATE(SS%par_pui(SS%pui_n_par, n))
            ALLOCATE(SS%pui_Sm(SS%pui_nW, n))
            ALLOCATE(SS%pui_Sp(SS%pui_nW, n))
            ALLOCATE(SS%pui_lock(n))


            SS%f_pui = 0.0
            SS%f_pui_num = 0.0
            SS%f_pui_num2 = 0.0
            SS%par_pui = 0.0
            SS%pui_Sm = 0.0
            SS%pui_Sp = 0.0

            do i = 1, n
                call omp_init_lock(SS%pui_lock(i))
            end do

            ! Теперь надо установить связь между номерами ячееки пуи и номерами в сетке
            n = 1
            do i = 1, size(SS%gl_all_Cell_zone(:))
                j = SS%gl_all_Cell_zone(i)
                if(j <= 2) then
                    SS%f_pui_num(n) = i
                    SS%f_pui_num2(i) = n
                    n = n + 1
                else
                    SS%f_pui_num2(i) = 0
                end if
            end do
        end if

    end subroutine PUI_SET


end module STORAGE


