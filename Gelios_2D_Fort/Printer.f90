
module Printer
    ! Модуль для печати всего и вся
    use STORAGE 
    use My_func
    USE ieee_arithmetic
    Use Phys_parameter

    contains 

    subroutine Print_hydrogen_1D(SS, number)
        ! Печатаем центры всех ячеек
        USE Phys_parameter
        TYPE (Setka), intent(in) :: SS
        integer :: i, j, node, num
        character(len=1) :: name
        character(len=5) :: name2
        integer(4), intent(in), OPTIONAL :: number
        real(8) :: source(4)

        num = 0
        if(PRESENT(number)) num = number
        write(unit=name2,fmt='(i5.5)') num

        open(1, file = 'Print_hydrogen_1D_' // name2 // '.txt')
        write(1,*) "TITLE = 'HP'  VARIABLES = X, rho, p, u, v" 

        do i = 1, SS%n_Hidrogen
            write(unit=name,fmt='(i1.1)') i
            write(1,*) ", Rho" // name // ", p" // name // ", u" // name // ", v" // name // ", T" // name
        end do

        write(1,*)", k1, k2, k3, In, Iu, Iv, IT, Mn, Mu, Mv, MT"

        do i = 1, SS%n_Hidrogen
            write(unit=name,fmt='(i1.1)') i
            write(1,*) ", IT" // name
        end do

        do i = size(SS%gl_Cell_B(:, 1)), 1, -1
            j = SS%gl_Cell_B(i, 1)
            call Calc_sourse_MF(SS, j, source, 1, use_koeff_ = .False.)
            write(1,*) SS%gl_Cell_Centr(1, j, 1), SS%gd(1:4, j, 1), SS%hydrogen(1:5, :, j, 1), SS%atom_source(:, j), source, SS%atom_all_source(4, :, j)
        end do

        do i = 1, size(SS%gl_Cell_A(:, 1))
            j = SS%gl_Cell_A(i, 1)
            ! if(SS%gl_Cell_Centr(1, j, 1) > 80) CYCLE
            ! if(SS%gl_Cell_Centr(1, j, 1) < 78) CYCLE

            call Calc_sourse_MF(SS, j, source, 1, use_koeff_ = .False.)
            ! print*, "source = ", source
            ! pause
            write(1,*) SS%gl_Cell_Centr(1, j, 1), SS%gd(1:4, j, 1), SS%hydrogen(1:5, :, j, 1),  SS%atom_source(:, j), source, SS%atom_all_source(4, :, j)
        end do

        close(1)

    end subroutine Print_hydrogen_1D

    subroutine Print_GD_1D(SS, number)
        ! Печатаем центры всех ячеек
        integer(4), intent(in), OPTIONAL :: number
        TYPE (Setka), intent(in) :: SS
        integer :: i, j, node, al, num, gr
        real(8) :: Mach, c(2), ro_He, par1(SS%n_par), par2(SS%n_par)
        character(len=5) :: name
        real(8) :: source(4)

        source = 0.0

        num = 0
        if(PRESENT(number)) num = number
        write(unit=name,fmt='(i5.5)') num

        open(1, file = name // '_Print_GD_1D.txt')
        write(1,*) "TITLE = 'HP'  VARIABLES = X, Rho, p, u, v, Q, Mach, Rho_He, In, Iu, Iv, IT"

        do i = size(SS%gl_Cell_B(:, 1)), 1, -1
            j = SS%gl_Cell_B(i, 1)
            Mach = norm2(SS%gd(3:4, j, 1))/sqrt(SS%par_ggg * SS%gd(2, j, 1)/SS%gd(1, j, 1))
			c = SS%gl_Cell_Centr(1, j, 1)
            call Calc_sourse_MF(SS, j, source, 1, use_koeff_ = .False.)

            ro_He = 0.0
            if(SS%n_par >= 6) ro_He = SS%gd(6, j, 1)
            al = 1
            if(SS%gl_all_Cell_zone(j) <= 2) al = 2
            call Sootnosheniya(SS%gd(1, j, 1), SS%gd(2, j, 1), ro_He, 0.0_8, 0.0_8, al)

            write(1,*) SS%gl_Cell_Centr(1, j, 1), SS%gd(1:4, j, 1), SS%gd(5, j, 1)/SS%gd(1, j, 1), Mach, ro_He, source
        end do

        do i = 1, SS%par_n_TS - 1 ! size(SS%gl_Cell_A(:, 1))
            j = SS%gl_Cell_A(i, 1)
            Mach = norm2(SS%gd(3:4, j, 1))/sqrt(SS%par_ggg * SS%gd(2, j, 1)/SS%gd(1, j, 1))
            call Calc_sourse_MF(SS, j, source, 1, use_koeff_ = .False.)
            ro_He = 0.0
            if(SS%n_par >= 6) ro_He = SS%gd(6, j, 1)
            al = 1
            if(SS%gl_all_Cell_zone(j) <= 2) al = 2
            call Sootnosheniya(SS%gd(1, j, 1), SS%gd(2, j, 1), ro_He, 0.0_8, 0.0_8, al)

            write(1,*) SS%gl_Cell_Centr(1, j, 1), SS%gd(1:4, j, 1), SS%gd(5, j, 1)/SS%gd(1, j, 1), Mach, ro_He, source
        end do

        ! Отдельно печатаем TS
        do i = SS%par_n_TS - 1, SS%par_n_TS - 1 ! size(SS%gl_Cell_A(:, 1))
            j = SS%gl_Cell_A(i, 1)
            gr = SS%gl_Cell_gran(1, j)
            call Calc_sourse_MF(SS, j, source, 1, use_koeff_ = .False.)
            call Get_gran_parameter(SS, gr, j, par1, par2, 1)

            Mach = norm2(par1(3:4))/sqrt(SS%par_ggg * par1(2)/par1(1))
            ro_He = 0.0
            if(SS%n_par >= 6) ro_He = par1(6)
            al = 1
            if(SS%gl_all_Cell_zone(j) <= 2) al = 2
            call Sootnosheniya(par1(1), par1(2), ro_He, 0.0_8, 0.0_8, al)
            write(1,*) SS%gl_Gran_Center(1, gr, 1), par1(1:4), par1(5)/par1(1), Mach, ro_He, source
            
            call Calc_sourse_MF(SS, SS%gl_Cell_A(i + 1, 1), source, 1, use_koeff_ = .False.)
            Mach = norm2(par2(3:4))/sqrt(SS%par_ggg * par2(2)/par2(1))
            ro_He = 0.0
            if(SS%n_par >= 6) ro_He = par2(6)
            al = 1
            if(SS%gl_all_Cell_zone(j) <= 2) al = 2
            call Sootnosheniya(par2(1), par2(2), ro_He, 0.0_8, 0.0_8, al)
            write(1,*) SS%gl_Gran_Center(1, gr, 1), par2(1:4), par2(5)/par2(1), Mach, ro_He, source
        end do

        do i = SS%par_n_TS, SS%par_n_HP - 1 ! size(SS%gl_Cell_A(:, 1))
            j = SS%gl_Cell_A(i, 1)
            Mach = norm2(SS%gd(3:4, j, 1))/sqrt(SS%par_ggg * SS%gd(2, j, 1)/SS%gd(1, j, 1))
            call Calc_sourse_MF(SS, j, source, 1, use_koeff_ = .False.)
            ro_He = 0.0
            if(SS%n_par >= 6) ro_He = SS%gd(6, j, 1)
            al = 1
            if(SS%gl_all_Cell_zone(j) <= 2) al = 2
            call Sootnosheniya(SS%gd(1, j, 1), SS%gd(2, j, 1), ro_He, 0.0_8, 0.0_8, al)

            write(1,*) SS%gl_Cell_Centr(1, j, 1), SS%gd(1:4, j, 1), SS%gd(5, j, 1)/SS%gd(1, j, 1), Mach, ro_He, source
        end do

        ! Отдельно печатаем HP
        do i = SS%par_n_HP - 1, SS%par_n_HP - 1 ! size(SS%gl_Cell_A(:, 1))
            j = SS%gl_Cell_A(i, 1)
            gr = SS%gl_Cell_gran(1, j)
            
            call Get_gran_parameter(SS, gr, j, par1, par2, 1)
            call Calc_sourse_MF(SS, j, source, 1, use_koeff_ = .False.)
            Mach = norm2(par1(3:4))/sqrt(SS%par_ggg * par1(2)/par1(1))
            ro_He = 0.0
            if(SS%n_par >= 6) ro_He = par1(6)
            al = 1
            if(SS%gl_all_Cell_zone(j) <= 2) al = 2
            call Sootnosheniya(par1(1), par1(2), ro_He, 0.0_8, 0.0_8, al)
            write(1,*) SS%gl_Gran_Center(1, gr, 1), par1(1:4), par1(5)/par1(1), Mach, ro_He, source

            Mach = norm2(par2(3:4))/sqrt(SS%par_ggg * par2(2)/par2(1))
            ro_He = 0.0
            if(SS%n_par >= 6) ro_He = par2(6)
            al = 1
            if(SS%gl_all_Cell_zone(j) <= 2) al = 2
            call Sootnosheniya(par2(1), par2(2), ro_He, 0.0_8, 0.0_8, al)
            write(1,*) SS%gl_Gran_Center(1, gr, 1), par2(1:4), par2(5)/par2(1), Mach, ro_He, source
        end do

        do i = SS%par_n_HP, SS%par_n_BS - 1 ! size(SS%gl_Cell_A(:, 1))
            j = SS%gl_Cell_A(i, 1)
            Mach = norm2(SS%gd(3:4, j, 1))/sqrt(SS%par_ggg * SS%gd(2, j, 1)/SS%gd(1, j, 1))
            call Calc_sourse_MF(SS, j, source, 1, use_koeff_ = .False.)
            ro_He = 0.0
            if(SS%n_par >= 6) ro_He = SS%gd(6, j, 1)
            al = 1
            if(SS%gl_all_Cell_zone(j) <= 2) al = 2
            call Sootnosheniya(SS%gd(1, j, 1), SS%gd(2, j, 1), ro_He, 0.0_8, 0.0_8, al)

            write(1,*) SS%gl_Cell_Centr(1, j, 1), SS%gd(1:4, j, 1), SS%gd(5, j, 1)/SS%gd(1, j, 1), Mach, ro_He, source
        end do

        ! Отдельно печатаем BS
        do i = SS%par_n_BS - 1, SS%par_n_BS - 1 ! size(SS%gl_Cell_A(:, 1))
            j = SS%gl_Cell_A(i, 1)
            gr = SS%gl_Cell_gran(1, j)
            
            call Get_gran_parameter(SS, gr, j, par1, par2, 1)
            call Calc_sourse_MF(SS, j, source, 1, use_koeff_ = .False.)
            Mach = norm2(par1(3:4))/sqrt(SS%par_ggg * par1(2)/par1(1))
            ro_He = 0.0
            if(SS%n_par >= 6) ro_He = par1(6)
            al = 1
            if(SS%gl_all_Cell_zone(j) <= 2) al = 2
            call Sootnosheniya(par1(1), par1(2), ro_He, 0.0_8, 0.0_8, al)
            write(1,*) SS%gl_Gran_Center(1, gr, 1), par1(1:4), par1(5)/par1(1), Mach, ro_He, source

            Mach = norm2(par2(3:4))/sqrt(SS%par_ggg * par2(2)/par2(1))
            ro_He = 0.0
            if(SS%n_par >= 6) ro_He = par2(6)
            al = 1
            if(SS%gl_all_Cell_zone(j) <= 2) al = 2
            call Sootnosheniya(par2(1), par2(2), ro_He, 0.0_8, 0.0_8, al)
            write(1,*) SS%gl_Gran_Center(1, gr, 1), par2(1:4), par2(5)/par2(1), Mach, ro_He, source
        end do

        do i = SS%par_n_BS, size(SS%gl_Cell_A(:, 1))
            j = SS%gl_Cell_A(i, 1)
            Mach = norm2(SS%gd(3:4, j, 1))/sqrt(SS%par_ggg * SS%gd(2, j, 1)/SS%gd(1, j, 1))
            call Calc_sourse_MF(SS, j, source, 1, use_koeff_ = .False.)
            ro_He = 0.0
            if(SS%n_par >= 6) ro_He = SS%gd(6, j, 1)
            al = 1
            if(SS%gl_all_Cell_zone(j) <= 2) al = 2
            call Sootnosheniya(SS%gd(1, j, 1), SS%gd(2, j, 1), ro_He, 0.0_8, 0.0_8, al)

            write(1,*) SS%gl_Cell_Centr(1, j, 1), SS%gd(1:4, j, 1), SS%gd(5, j, 1)/SS%gd(1, j, 1), Mach, ro_He, source
        end do



        close(1)

    end subroutine Print_GD_1D


end module Printer