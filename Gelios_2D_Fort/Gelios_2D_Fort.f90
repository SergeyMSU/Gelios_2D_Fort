!  Gelios_2D_Fort.f90 
!
!  FUNCTIONS:
!  Gelios_2D_Fort - Entry point of console application.
!
	
include "Storage.f90"
include "Geometry.f90"
	
	
program Gelios_2D_Fort
    USE STORAGE
	USE GEOMETRY
    implicit none

    gl_S1%name = "00001"
    call Init_Setka(gl_S1)
    call Build_Setka_start(gl_S1)
    call Print_Point_from_Rays(gl_S1)
    call Print_Cell(gl_S1)
    call Culc_Cell_Centr(gl_S1, 1)
    call Culc_Cell_Centr(gl_S1, 2)
    call Print_Cell_Centr(gl_S1)
    call Print_Connect(gl_S1)
    call Print_Grans(gl_S1)
    call Proverka_grans_sosed(gl_S1)

    ! print*, gl_S1%gl_Cell_gran(1, 1), gl_S1%gl_Cell_gran(2, 1), gl_S1%gl_Cell_gran(3, 1), gl_S1%gl_Cell_gran(4, 1)
    ! print*, gl_S1%gl_Cell_Centr(:, 1, 1)

    print *, gl_S1%init_geo
	pause

end program Gelios_2D_Fort

