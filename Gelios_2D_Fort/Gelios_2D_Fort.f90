!  Gelios_2D_Fort.f90 
!
!  FUNCTIONS:
!  Gelios_2D_Fort - Entry point of console application.
!
	
include "Storage.f90"
	
	
program Gelios_2D_Fort
    USE STORAGE
    implicit none

    gl_S1%name = "00001"
    call Init_Setka(gl_S1)
    call Build_Setka_start(gl_S1)
    call Print_Point_from_Rays(gl_S1)
    call Print_Cell(gl_S1)

    print *, gl_S1%init_geo
	pause

end program Gelios_2D_Fort

