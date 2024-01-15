module Phys_parameter
    use STORAGE 
    use GEOMETRY
    implicit none 

    contains

    subroutine Phys_Init(par, SS, info, grid)
        ! Функция инициализации физических параметров
        TYPE (Setka), intent(in) :: SS
        TYPE (Phys_par), intent(in out) :: par
        character(len=4), intent(in) :: info
        character(len=5), intent(in) :: grid
        integer(4) :: n

        par%info = info
        par%grid = grid

        if(par%grid /= SS%name) then
            STOP "Phys_parameter  Phys_Init  error name grid  67892iuytr6y43"  
        end if

        if(SS%init_geo == .False.) then
            STOP "Phys_parameter  Phys_Init  error init_geo  90-plkjuytg"  
        end if


        if(par%info == "cent") then
            n = size(SS%gl_all_Cell(1, :))
        else if(par%info == "yzel") then
            n = size(SS%gl_yzel(1, :, 1))
        else
            STOP "Phys_parameter  Phys_Init  error info  87432qasdrtyhgfd"  
        end if

        ALLOCATE(par%gd(4, n ,2))
        par%init = .True.

    end subroutine Phys_Init

    subroutine Inner_Conditions(SS, x, y, rho, p, u, v)
        !! Задаём внутренние граничные условия
        TYPE (Setka), intent(in) :: SS
        real(8), intent(in) :: x, y
        real(8), intent(out) :: rho, p, u, v
        real(8) :: r, p_0

        r = sqrt(x**2 + y**2)

        rho = 6.0 * (SS%par_R0/r)**2
        u = x/r * 430.0
        v = y/r * 430.0
        p_0 = 430.0**2 * 6.0**2/(SS%par_ggg * 10.0**2)
        p = p_0 * (SS%par_R0/r)**(2.0 * SS%par_ggg)
    end subroutine Inner_Conditions

end module Phys_parameter