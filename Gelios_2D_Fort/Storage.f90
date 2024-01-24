
module STORAGE
    implicit none 

    !! ����� ������������ �������� (������� ������� �� ����������)
    real(8), parameter :: par_pi = acos(-1.0_8) 

    !! ������ ������ ��� ����� �� ����� �����������
    TYPE Setka 

        character(len=5) :: name = "00000"
        logical :: init_geo = .False.   ! ���������������� �� ������ ����� (�������� �� ������ ��� ������� ���������)


        ! ����� ����������, ������������ ��������� �����
        integer(4) :: par_m_A = 20! 30      ! ���������� ����� A � ���������
        integer(4) :: par_m_BC = 10! 18      ! ���������� ����� B/C � ���������
        integer(4) :: par_m_O = 10! 17      ! ���������� ����� O � ���������
        integer(4) :: par_m_K = 8! 7      ! ���������� ����� K � ���������
        real(8) :: par_triple_point = 13.0 * par_pi/40.0     ! �� ������ ���� ������� �� pi/2 (� �������������� x) ������� �����
        real(8) :: par_triple_point_2 = 7.0 * par_pi/40.0     ! ��� ����� ����� ������� ��� ����� ������� ����� ������� �� pi/2 (� �������������� x) 
        
        ! ���������� ����� �� ����� A
        integer(4) :: par_n_TS =  33! 26                    ! ���������� ����� �� TS (TS ����������)
        integer(4) :: par_n_HP =  63! 40                 ! ���������� ����� �� HP (HP ����������)  �� �� 0 ���������
        integer(4) :: par_n_BS =  89! 60! 5                 ! ���������� ����� BS (BS ����������)
        integer(4) :: par_n_END = 98! 72! 6                ! ���������� ����� �� ����� ����� (����� ����������)
        integer(4) :: par_n_IA =  20! 12                   ! ���������� �����, ������� ������ �� ���������� �������
        integer(4) :: par_n_IB =  22! 14                   ! ���������� �����, ������� ������ �� ���������� ������� (� �������)

        ! ����� ����������, �������� ������� �����
        real(8) :: par_R_character = 35.0         ! ����������� ������ � ������ (���������� �� TS �� ��������� ����� ���������� �����)
        real(8) :: par_R0 = 0.197035         ! ����������� ������ 1 �.�. (���������� �����) ��� ��������� ������ ����� �� ����� �� ����� (������ ��������� � ����)
        real(8) :: par_R_END = 300.0         !  
        real(8) :: par_R_LEFT = -240.0 ! -390.0         !  ����� �������
        real(8) :: par_R_inner = 9.0! 5.0_8     ! �� ������ ���������� ���������� �����

        !! ���������� ��������� ---------------------------------------
        real(8) :: par_a_2 = 0.130738_8        ! �������� � ������� �����������
        real(8) :: par_ggg = 5.0/3.0                 ! �� ������ ���������� ���������� �����
        real(8) :: par_Velosity_inf = -2.54279_8
        real(8) :: par_n_H_LISM = 3.5_8
        real(8) :: par_Kn = 49.9018   !0.4326569808         ! � �����������
        !! -------------------------------------------------------------

        !����� ���������� ��������
        real(8) :: par_kk1 = 2.0_8     ! ������� �������� ����� � ���� � ������� �� TS: 1 - ��������, 2 - ������������ � �.�.
        real(8) :: par_kk2 = 1.7_8     ! ������� �������� � �������� ������� �� �������������
        real(8) :: par_kk3 = 1.8_8     ! ������� �������� � ������
        real(8) :: par_kk31 = 1.0_8     ! ������� �������� � ������ ��� ����� �� �������� (������ ����� � � - ����)
        real(8) :: par_kk13 = 1.8_8     ! ������� �������� ����� � �������� ������� �� ������� ������� ����  �� 0 �� 1
        real(8) :: par_kk131 = 0.1_8
        real(8) :: par_kk132 = 1.5_8
        real(8) :: par_kk14 = 1.0_8     ! ������� �������� ����� � �������� ������� �� ���������� ������� ����  �� 0 �� 1
        ! (�������� ����� � TS � HP)  
        real(8) :: par_kk12 = 1.0_8     ! ������� �������� ����� �� TS � ������� �����  >= 1
        ! ������ �������� �� 4 ��� �������� ������ ����������� � ����������

        real(8) :: par_nat_TS = 0.001_8   ! ������������� ���������
        real(8) :: par_nat_HP = 0.005_8   ! ������������� ���������
        real(8) :: par_nat_BS = 0.002_8   ! ������������� ���������

        integer(4) :: par_n_points                         ! ����� ����� � ����� (��������� ��� ������������� �����)

        real(8), allocatable :: gl_yzel(:, :, :)   ! (2, :, 2) ����� ��������� ����� �����
        real(8), allocatable :: gl_yzel_Vel(:, :)   ! (2, :) �������� �������� ����� �����  !TODO NO-SAVE
        ! ���� ������ �� ���� ��������� ��� ���������� ����� � ���� - ��� ������������� ������� ������
        integer(4), allocatable :: gl_Point_num(:)   ! ������� ��� �������� �������� � ����  !TODO NO-SAVE

        ! ����, �� ������� ������������� ����� �����
        integer(4), allocatable :: gl_RAY_A(:,:)   ! ����� �-����� ����������� 3 (�� ���� ����, � ���� ���������)
        integer(4), allocatable :: gl_RAY_B(:,:)   ! ����� B-����� ����������� 3 (�� ���� ����, � ���� ���������)
        integer(4), allocatable :: gl_RAY_C(:,:)   ! ����� C-����� ����������� 3 (�� ���� ����, � ���� ���������)
        integer(4), allocatable :: gl_RAY_O(:,:)   ! ����� O-����� ����������� 3 (�� ���� ����, � ���� ���������)
        integer(4), allocatable :: gl_RAY_K(:,:)   ! ����� K-����� ����������� 3 (�� ���� ����, � ���� ���������)
        integer(4), allocatable :: gl_RAY_D(:,:)   ! ����� D-����� ����������� 3 (�� ���� ����, � ���� ���������)
        integer(4), allocatable :: gl_RAY_E(:,:)   ! ����� E-����� ����������� 3 (�� ���� ����, � ���� ���������)

        ! ������
        integer(4), allocatable :: gl_Cell_A(:,:)   ! ����� A-����� ����������� 3 (�� ���� ����, � ���� ���������)
        integer(4), allocatable :: gl_Cell_B(:,:)   ! ����� B-����� ����������� 3 (�� ���� ����, � ���� ���������)
        integer(4), allocatable :: gl_Cell_C(:,:)   ! ����� C-����� ����������� 3 (�� ���� ����, � ���� ���������)

        integer(4), allocatable :: gl_all_Cell(:,:)   ! ���� ����� ����� (4, :) - ������ ���������� ������� - ��� ����� ����� ������

        integer(4), allocatable :: gl_Cell_neighbour(:,:)   ! (4, :) ����� �� 4 ������� ��� ������ ������  !! ������ ����������� � �������� ��� ������!
        ! -1   ! ������� (���������� �����)
        ! -2   ! �������� �������
        ! -3   ! ������� ������� �������
        ! -4   ! ��� ���������

        integer(4), allocatable :: gl_all_Cell_zone(:)   ! ���� ������
        ! 1, 2, 3, 4

        integer(4), allocatable :: gl_Cell_gran(:,:)        ! (4, :) ����� �� 4 ������ ��� ������ ������ (���� ����� = 0, �� ����� ��� � ���� �����������)
        ! 0 - ��� ����� (����� ������ ����, �������� ����� ����)
        !! ������������� ������� � ����� ���!
        real(8), allocatable :: gl_Cell_gran_dist(:, :, :)      ! (4, :, 2) ���������� �� ����� ����� �� ������� ������

        real(8), allocatable :: gl_Cell_Centr(:, :, :)   ! (2, : ����� �����, 2) ����� ��������� ������� �����
        

        integer(4), allocatable :: gl_all_Gran(:,:)       ! ��� ����� (2,:) ����� �� 2 ����
        integer(4), allocatable :: gl_Gran_neighbour(:,:) ! ������ ������ ����� (2,:) ����� �� 2 ������, ������� ���� �� ������� �� �������
        real(8), allocatable :: gl_Gran_normal(:, :, :)       ! (2, :, 2) ������� �����     !! ������� ��� �� ������ ������ �� ������ �����������!                   
        real(8), allocatable :: gl_Gran_length(:, :)       ! (:, 2) ����� �����                       
        real(8), allocatable :: gl_Gran_Center(:, :, :)       ! (2, :, 2) ����� �����                       
        
        integer(4), allocatable :: gl_Gran_neighbour_TVD(:,:) ! TVD-������ ������ ����� (2,:) ����� �� 2 ������
        ! 0 - ������ ������ ���
        ! ��� ���� ������ TVD-����� - ��� ����� ������� �������� ������. �.�. ������� ����� ���� ���� �� ������� �� ������� TVD-������

        
        real(8), allocatable :: gl_Cell_belong(:,:,:)      ! (3, 4, :)  ! ������������ ��� ��������� �� ������ ���� �� �������
        ! ��� ������ ����� ������������ A, B, C   Ax + By + C = 0  (���� ������ 0, �� ����� ��� ������! �.�. ������� �������)
        real(8), allocatable :: gl_Cell_square(:,:)      ! (:, 2)

        character, allocatable :: gl_Cell_type(:)           ! ��� ������ ������ �, �, �
        integer(4), allocatable :: gl_Cell_number(:, :)     ! (2, :) ����� ������ ������ ������ ������ ����

        ! ����������� ���������
        integer(4), allocatable :: gl_HP(:)            ! ������� - ������ ������, ������� �������� ���������� �����������
        ! ������� ������� �� ������� ������
        integer(4), allocatable :: gl_TS(:)
        integer(4), allocatable :: gl_BS(:)
        !? ���������, ����� ������� � ������-������������ ���� ������������� "������"
        !? ���������, ��� � ������ ������ � ������ ���� ��������� ����������� (� ����� �� ����)

        integer(4), allocatable :: gl_Gran_type(:)      ! ���������� ��� ����� 
        ! (0 - �������, 1 - TS, 2 - HP, 3 - BS)     

        !! ������
        integer(4) :: n_Hidrogen = 4  ! ����� ������ ������ ��������
        integer(4) :: n_par = 5  !! ����� ���������� ���������� � ������
        ! 5 ����������������
        ! 4 * n_Hidrogen - �������

        real(8), allocatable :: gd(:, :, :)  ! (n_par, :, 2 ��������� ����)
        ! (rho p u v Q)
        real(8), allocatable :: hydrogen(:, :, :, :)  ! (5, n_Hidrogen, :, 2 ��������� ����)
        ! rho p u v T

    END TYPE Setka

    TYPE Inter_Setka  ! ����� ��� ������������

        logical :: init = .False.   ! ���������������� �� ������ ����� (�������� �� ������ ��� ������� ���������)
        real(8), allocatable :: gl_yzel(:, :)   ! (2, :) ����� ��������� ����� �����

        !? �������������� ����
        ! ����� ��� ������ �������� �� 4 ������, �� ���� ��������� �������������
        ! � ���� ������ �������� ���� ����� ������� (��� ���� ����� ������)
        ! � ����������� �������� �������� ����� ��������

        ! ������
        !! ����� � � B ������ ������������ (��� ��� �� ����� A � B)
        integer(4), allocatable :: gl_Cell_A(:,:)   ! ����� A-����� ����������� 3 (�� ���� ����, � ���� ���������)
        integer(4), allocatable :: gl_Cell_B(:,:)   ! ����� B-����� ����������� 3 (�� ���� ����, � ���� ���������)
        integer(4), allocatable :: gl_Cell_C(:,:)   ! ����� C-����� ����������� 3 (�� ���� ����, � ���� ���������)

        integer(4), allocatable :: gl_all_Cell(:,:)   ! ���� ����� ����� (4, :) - ������ ���������� ������� - ��� ����� ����� ������
        real(8), allocatable :: gl_Cell_center(:,:)   ! ���� ����� ����� (2, :) - ������ ���������� ������� - ��� ����� ����� ������

        logical, allocatable :: gl_all_triangle(:)
        ! .True. - ���� ������ �������� ������������� � .False. � ��������� ������

        integer(4), allocatable :: gl_Cell_neighbour(:,:)   ! (4, :) ����� �� 4 ������� ��� ������ ������ 
        real(8), allocatable :: gl_Cell_Belong(:, :, :)       ! ��� ����� (3, 4, :) ����� �� A B C, 4 ����� � ������, 
        ! Ax + By + C = 0  (���� > 0 �� �� ��������� ������)

        real(4), allocatable :: gl_Cell_interpol_matrix(:, :, :)   ! (4, 4, :) ����� �� 4 ������� ��� ������ ������
        ! Ax + By + Cxy + D
        ! Ax + By + C
        ! 
        ! ���������������� ������� ��� ������ ������

        !! ������
        integer(4) :: n_Hidrogen = 4  ! ����� ������ ������ ��������
        integer(4) :: n_par = 5  !! ����� ���������� ���������� � ������
        ! 5 ����������������
        ! 4 * n_Hidrogen - �������

        real(8), allocatable :: gd(:, :)  ! (n_par, :)
        ! (rho p u v Q)
        real(8), allocatable :: hydrogen(:, :, :)  ! (5, n_Hidrogen, :)
        ! rho p u v T

    END TYPE Inter_Setka


    TYPE Surfaces
        ! ������ �������� ������������ ��� �������� ����� � ���� ������������
        ! ����� ���� �����, ��������, ��� ����������� ������� ����� (��� ��� ��������� ������ ������, �� ����� ��������� ������ �����)
        ! � ��������� � �����������
        logical :: init = .False.

        real(8), allocatable :: TS(:, :)  ! (2, :) ����, ������
        real(8), allocatable :: HP(:, :)  ! (4, :) ����, ������, x, y
        real(8), allocatable :: BS(:, :)  ! (4, :) ����, ������, x, y

    END TYPE Surfaces
	
	!! ����� ���������� ���������� 
    TYPE (Setka):: gl_S1
    TYPE (Setka):: gl_S3

    TYPE (Inter_Setka):: gl_S2

    TYPE (Surfaces):: gl_surf1


    contains 


end module STORAGE



! N2 = size(SS%gl_Cell_A(1, :))
! N1 = size(SS%gl_Cell_A(:, 1))

! 2.0 * N1 * N2 - N1

! N2 = size(SS%gl_Cell_B(1, :))
! N1 = size(SS%gl_Cell_B(:, 1))

! 2.0 * N1 * N2

! N2 = size(SS%gl_Cell_C(1, :))
! N1 = size(SS%gl_Cell_C(:, 1))

! 2.0 * N1 * N2 + N1