!  Gelios_2D_Fort.f90 
!
!  FUNCTIONS:
!  Gelios_2D_Fort - Entry point of console application.
!
	
include "cgod_3D.f90"
include "Storage.f90"
include "My_func.f90"
include "Solvers.f90"
include "Geometry.f90"
include "Phys_parameter.f90"
include "Surfaces.f90"
include "Interpol.f90"
include "M-K.f90"
include "2D_algoritm.f90"

	
	
program Gelios_2D_Fort
    USE STORAGE
	USE GEOMETRY
    USE SURFACE
    USE Phys_parameter
    USE Interpol
    USE Monte_Karlo
    USE Algoritm
	USE My_func

    implicit none

    ! integer(4) :: num, mi, ni
    ! real(8) :: M(4, 4)
    ! real(8) :: B(4, 4)

    ! M = 0.0
    ! M(1, 1) = 1.0
    ! M(2, 1) = 1.0
    ! M(2, 2) = 1.0
    ! M(3, 3) = 1.0
    ! M(4, 4) = 1.0
    ! M(2, 3) = 1.0

    ! B = matinv4(M)
	
	! do mi = 1, 4
    !     do ni = 1, 4
    !         write(*,"(F12.2,$)") M(mi,ni)
    !     end do
    !     write (*,*) ''
    ! end do

    ! print*, "__________________"

    ! do mi = 1, 4
    !     do ni = 1, 4
    !         write(*,"(F12.2,$)") B(mi,ni)
    !     end do
    !     write (*,*) ''
    ! end do
	
	
	! pause


    ! gl_S1%name = "00001"
    ! call Init_Setka(gl_S1)
    ! call Build_Setka_start(gl_S1)
    ! call Print_Point_from_Rays(gl_S1)
    !call Print_Cell(gl_S1)
    ! call Culc_Cell_Centr(gl_S1, 1)
    ! call Culc_Cell_Centr(gl_S1, 2)
    ! call Print_Cell_Centr(gl_S1)
    ! call Print_Connect(gl_S1)
    ! call Print_Grans(gl_S1)
    ! call Proverka_grans_sosed(gl_S1)
    ! call Geo_Print_Surface(gl_S1)


    ! call SUR_Download(gl_surf1, "00000")
	
	! call Algoritm_ReMove_Surface(gl_S1, gl_surf1)
    ! call Print_Grans(gl_S1)
	! call Geo_Print_Surface(gl_S1)

    ! call Int_Init(gl_S2, gl_S1)
    ! call Int_Print_all_point(gl_S2)
    ! call Int_Print_Cell(gl_S2)

    ! call Int_Print_center(gl_S2)
    ! call Int_Print_connect(gl_S2)

    ! call Proverka_grans_sosed(gl_S1)
    
    ! call Geo_request(gl_S1)

    !call Gas_dynamic_algoritm2(gl_S1)
    call MK_algoritm(gl_S1)
    !call Perestroika_algoritm(gl_S1)

    !call Print_Cell(gl_S1)

    ! num = 1
    ! call Geo_Find_Cell(gl_S1, -50.0_8, 170.0_8, num)
    ! print*, "Num 1 = ", num

    ! num = 1
    ! call Int_Find_Cell(gl_S2, 8.0_8, 8.0_8, num)
    ! print*, "Num 2 = ", num

    ! print*, gl_S1%gl_Cell_gran(1, 1), gl_S1%gl_Cell_gran(2, 1), gl_S1%gl_Cell_gran(3, 1), gl_S1%gl_Cell_gran(4, 1)
    ! print*, gl_S1%gl_Cell_Centr(:, 1, 1)


end program Gelios_2D_Fort

