module SURFACE
    USE STORAGE
    USE GEOMETRY
    USE My_func
    implicit none 

    contains

    subroutine SUR_Download(SURF, name)
        TYPE (Surfaces), intent(in out) :: SURF
        character(len=5), intent(in) :: name
        logical :: exists
        integer(4) :: n, i
        real(8) :: alp, r, x, y, dist

        dist = 78.0  !! Коэффициент масштабирования

        print*, "Start SUR_Download"
        inquire(file= "surface_" // name // ".bin", exist=exists)
        if (exists == .False.) then
            STOP "ERROR SUR_Download net faila!  953gfiyfgefiguewehrywe4f"
        end if

        if(SURF%init == .True.) then
            STOP "ERROR SUR_Download SURF%init!  -irhgty6wtrfdvbfje"
        end if

        SURF%init = .True.

        open(1, file = "surface_" // name // ".bin", FORM = 'BINARY', ACTION = "READ")

        ! ----------------------------- TS
        read(1) n
        print*, "n1 = ", n

        ALLOCATE(SURF%TS(2, n))

        do i = 1, n
            read(1) alp
            read(1) r
            SURF%TS(1, i) = alp
            SURF%TS(2, i) = r * dist
        end do

        ! ----------------------------- HP
        read(1) n
        print*, "n2 = ", n

        ALLOCATE(SURF%HP(4, n))

        do i = 1, n
            read(1) x
            read(1) y
            read(1) alp
            read(1) r
            SURF%HP(1, i) = alp
            SURF%HP(2, i) = r * dist
            SURF%HP(3, i) = x * dist
            SURF%HP(4, i) = y * dist
        end do

        ! ----------------------------- HP
        read(1) n
        print*, "n3 = ", n

        ALLOCATE(SURF%BS(4, n))

        do i = 1, n
            read(1) x
            read(1) y
            read(1) alp
            read(1) r
            SURF%BS(1, i) = alp
            SURF%BS(2, i) = r * dist
            SURF%BS(3, i) = x * dist
            SURF%BS(4, i) = y * dist
        end do

        close(1)

        print*, "END SUR_Download"

	end subroutine SUR_Download

    subroutine SUR_Save(SURF, name)
        TYPE (Surfaces), intent(in) :: SURF
        character(len=5), intent(in) :: name
        logical :: exists
        integer(4) :: n, i
        real(8) :: alp, r, x, y


        print*, "Start SUR_Save"

        if(SURF%init == .False.) then
            STOP "ERROR SUR_Save SURF%init!  -irhgty6wtrfdvbfje"
        end if

        open(1, file = "surface_" // name // ".bin", FORM = 'BINARY')

        ! ----------------------------- TS
        n = size(SURF%TS(1, :))
        write(1) n

        do i = 1, n
            write(1) SURF%TS(1, i)
            write(1) SURF%TS(2, i)
        end do

        ! ----------------------------- HP
        n = size(SURF%HP(1, :))
        write(1) n

        do i = 1, n
            write(1) SURF%HP(3, i)
            write(1) SURF%HP(4, i)
            write(1) SURF%HP(1, i)
            write(1) SURF%HP(2, i)
        end do

        ! ----------------------------- HP
        n = size(SURF%BS(1, :))
        write(1) n

        do i = 1, n
            write(1) SURF%BS(3, i) 
            write(1) SURF%BS(4, i) 
            write(1) SURF%BS(1, i) 
            write(1) SURF%BS(2, i) 
        end do

        close(1)

        print*, "END SUR_Save"

	end subroutine SUR_Save

    subroutine SUR_init(SURF, SS)
        TYPE (Surfaces), intent(in out) :: SURF
        TYPE (Setka), intent(in) :: SS

        integer(4) :: i, j, k, n, gr
        real(8) :: c(2), r, alp, x, y

        n = size(SS%gl_TS)
        allocate(SURF%TS(2, n))

        do i = 1, n
            gr = SS%gl_TS(i)
            c = SS%gl_Gran_Center(:, gr, 1)
            SURF%TS(1, i) = polar_angle(c(1), c(2))
            SURF%TS(2, i) = norm2(c)
        end do


        n = size(SS%gl_HP)
        allocate(SURF%HP(4, n))

        do i = 1, n
            gr = SS%gl_HP(i)
            c = SS%gl_Gran_Center(:, gr, 1)
            SURF%HP(1, i) = polar_angle(c(1), c(2))
            SURF%HP(2, i) = norm2(c)
            SURF%HP(3:4, i) = c
        end do


        n = size(SS%gl_BS)
        allocate(SURF%BS(4, n))

        do i = 1, n
            gr = SS%gl_BS(i)
            c = SS%gl_Gran_Center(:, gr, 1)
            SURF%BS(1, i) = polar_angle(c(1), c(2))
            SURF%BS(2, i) = norm2(c)
            SURF%BS(3:4, i) = c
        end do

    end subroutine SUR_init

	! real(8) pure function SUR_GET_TS(SURF, alp)
    real(8) function SUR_GET_TS(SURF, alp)
        ! Получить r_TS по углу
        TYPE (Surfaces), intent(in) :: SURF
        real(8), intent(in) :: alp
        real(8) :: al1, al2, r1, r2
        integer(4) n, i

        n = size(SURF%TS(1, :))

        if (SURF%TS(1, 1) >= alp) then
            SUR_GET_TS = SURF%TS(2, 1)
            return
        else if (SURF%TS(1, n) <= alp) then
            SUR_GET_TS = SURF%TS(2, n)
            return
        else
            do i = 1, n
                if (alp < SURF%TS(1, i)) then
                    al1 = SURF%TS(1, i - 1)
                    al2 = SURF%TS(1, i)
                    r1 = SURF%TS(2, i - 1)
                    r2 = SURF%TS(2, i)
                    EXIT
                end if
            end do
            SUR_GET_TS = linear1(al1, r1, al2, r2, alp)
            return
        end if

        print*, "ERROR SUR_GET_TS 112 0987ytghjkowiuhfwefw"
    end function SUR_GET_TS

    real(8) function SUR_GET_HP_alp(SURF, alp)
        ! Получить r_TS по углу
        TYPE (Surfaces), intent(in) :: SURF
        real(8), intent(in) :: alp
        real(8) :: al1, al2, r1, r2
        integer(4) n, i

        n = size(SURF%HP(1, :))

        if (SURF%HP(1, 1) >= alp) then
            SUR_GET_HP_alp = SURF%HP(2, 1)
            return
        else if (SURF%HP(1, n) <= alp) then
            SUR_GET_HP_alp = SURF%HP(2, n)
            return
        else
            do i = 1, n
                if (alp < SURF%HP(1, i)) then
                    al1 = SURF%HP(1, i - 1)
                    al2 = SURF%HP(1, i)
                    r1 = SURF%HP(2, i - 1)
                    r2 = SURF%HP(2, i)
                    EXIT
                end if
            end do
            SUR_GET_HP_alp = linear1(al1, r1, al2, r2, alp)
            return
        end if

        print*, "ERROR SUR_GET_HP_alp 112 5egiouyfeewfwef"
    end function SUR_GET_HP_alp

    real(8) function SUR_GET_HP_x(SURF, x)
        ! Получить r_TS по углу
        TYPE (Surfaces), intent(in) :: SURF
        real(8), intent(in) :: x
        real(8) :: x1, x2, y1, y2
        integer(4) n, i

        n = size(SURF%HP(1, :))

        if (SURF%HP(3, 1) <= x) then
            SUR_GET_HP_x = SURF%HP(4, 1)
            return
        else if (SURF%HP(3, n) >= x) then
            SUR_GET_HP_x = SURF%HP(4, n)
            return
        else
            do i = 1, n
                if (x > SURF%HP(3, i)) then
                    x1 = SURF%HP(3, i - 1)
                    x2 = SURF%HP(3, i)
                    y1 = SURF%HP(4, i - 1)
                    y2 = SURF%HP(4, i)
                    EXIT
                end if
            end do
            SUR_GET_HP_x = linear1(x1, y1, x2, y2, x)
            return
        end if

        print*, "ERROR SUR_GET_HP_x 112 5egiouyfeewfwef"
    end function SUR_GET_HP_x

    real(8) function SUR_GET_BS_alp(SURF, alp)
        ! Получить r_TS по углу
        TYPE (Surfaces), intent(in) :: SURF
        real(8), intent(in) :: alp
        real(8) :: al1, al2, r1, r2
        integer(4) n, i

        n = size(SURF%BS(1, :))

        if (SURF%BS(1, 1) >= alp) then
            SUR_GET_BS_alp = SURF%BS(2, 1)
            return
        else if (SURF%BS(1, n) <= alp) then
            SUR_GET_BS_alp = SURF%BS(2, n)
            return
        else
            do i = 1, n
                if (alp < SURF%BS(1, i)) then
                    al1 = SURF%BS(1, i - 1)
                    al2 = SURF%BS(1, i)
                    r1 = SURF%BS(2, i - 1)
                    r2 = SURF%BS(2, i)
                    EXIT
                end if
            end do
            SUR_GET_BS_alp = linear1(al1, r1, al2, r2, alp)
            return
        end if

        print*, "ERROR SUR_GET_BS_alp 112 5egiouyfeewfwef"
    end function SUR_GET_BS_alp

end module SURFACE