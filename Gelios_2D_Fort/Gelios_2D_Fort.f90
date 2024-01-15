!  Gelios_2D_Fort.f90 
!
!  FUNCTIONS:
!  Gelios_2D_Fort - Entry point of console application.
!
	
include "Storage.f90"
include "Geometry.f90"
include "Phys_parameter.f90"
include "2D_algoritm.f90"
	
	
program Gelios_2D_Fort
    USE STORAGE
	USE GEOMETRY
    USE Phys_parameter
    USE Algoritm

    implicit none

    integer(4) :: num

    gl_S1%name = "00001"
    call Init_Setka(gl_S1)
    call Build_Setka_start(gl_S1)
    call Print_Point_from_Rays(gl_S1)
    call Print_Cell(gl_S1)
    ! call Culc_Cell_Centr(gl_S1, 1)
    ! call Culc_Cell_Centr(gl_S1, 2)
    call Print_Cell_Centr(gl_S1)
    call Print_Connect(gl_S1)
    call Print_Grans(gl_S1)
    call Proverka_grans_sosed(gl_S1)
    call Geo_Print_Surface(gl_S1)

    call Phys_Init(gl_par1, gl_S1, "cent", "00001")

    call Gas_dynamic_algoritm(gl_S1, gl_par1)

    num = 1
    call Geo_Find_Cell(gl_S1, 3.0_8, 3.0_8, num)
    print*, "Num = ", num

    ! print*, gl_S1%gl_Cell_gran(1, 1), gl_S1%gl_Cell_gran(2, 1), gl_S1%gl_Cell_gran(3, 1), gl_S1%gl_Cell_gran(4, 1)
    ! print*, gl_S1%gl_Cell_Centr(:, 1, 1)

	pause

end program Gelios_2D_Fort

