
module STORAGE
    USE OMP_LIB
    implicit none 

    !! ����� ������������ �������� (������� ������� �� ����������)
    real(8), parameter :: par_pi = acos(-1.0_8) 
    real(8), parameter :: par_sqrtpi = sqrt(par_pi)

    integer(4), parameter :: par_n_zone = 6! 7  !  ���������� �������� (�� ���� ��� ������� ����)
	integer(4), parameter :: par_m_zone = 7! 6  !  ���������� ����� �� ���� (�� 0 �� 180)
    integer(4), parameter :: par_n_potok = 32! 32! 24! 32  ! ����� ������� (� ������� ������ ���� ����)
    integer(4), parameter :: par_n_claster = 1  ! ����� ����������� (��� MPI)
    integer(4), parameter :: par_n_parallel = 20! 20  ! ��� ����������������� ����� (�.�. ������ ����� ����� � ������� ������������ ����� ����� ��������
    integer(4), parameter :: par_stek = 1000  ! ������� ����� (������� ���������� ������ ��� ����)
    integer(4), parameter :: par_n_sort = 4!4!6  !  4 ���������� ������ ������

    logical, parameter :: MK_is_NaN = .False.    ! ����� �� �������� �� nan
	logical, parameter :: MK_Mu_stat = .True.    ! ����� �� ����������� ���� ��� ���������� � ������� �������������
	logical, parameter :: MK_photoionization = .True.    ! ����� �� �������������
	logical, parameter :: MK_el_impact = .False.    ! ����� �� ����������� ����
	logical, parameter :: MK_statistik_file = .True.    ! ����� ���� �� ����������� ������������?  true - ��� ������ ����

	logical, parameter :: par_Hydro = .True.    ! ����� �� �������? ��� ����� �������-�������� �������?
	logical, parameter :: par_TVD_linear_HP = .True.! .False. !.True.    ! ����� �� ������� �� ������� ������� ������ � ����� �������
	logical, parameter :: par_move_setka = .True.    ! ����� �� ������� �� ������� ������� ������ � ����� �������
    real(8), parameter :: par_Rmax = 220.0 !300.0! 220.0  !  ������ �����, � ������� ��������� �������


    ! ����� ������ � ������� ������!
	! ����� ������ ���� ������ par_n_parallel
	integer(4), parameter :: MK_k_multiply = 12 * 6! 12 * 2!11 * 8!12 * 10! 12 * 7!14 * 8!6 * 3! * 6 * 9! * 6 * 8!6 * 6 * 2  !   ! 6 = 20 ����� ����� (� �������� 30 �����)
    ! 9 ������ � ��������
    ! 12 (14) - ��� 1 ��� � ��������
    ! 18 - ��� 1 ��� ��� �������
	integer(4), parameter :: MK_k_mul1 = 6 * MK_k_multiply! 6
	integer(4), parameter :: MK_k_mul2 = 1 * MK_k_multiply! 
	integer(4), parameter :: MK_N1 = MK_k_mul1 * 60/par_n_parallel   ! 60 ����� �������� ������ ������� ���� (� ���������)
	integer(4), parameter :: MK_N2 = MK_k_mul1 * 20/par_n_parallel   ! 20
	integer(4), parameter :: MK_N3 = MK_k_mul2 * 20/par_n_parallel   ! 20
	integer(4), parameter :: MK_N4 = MK_k_mul1 * 20/par_n_parallel   ! 20

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
        real(8) :: par_R0 = 0.198956         ! ����������� ������ 1 �.�. (���������� �����) ��� ��������� ������ ����� �� ����� �� ����� (������ ��������� � ����)
        real(8) :: par_R_END = 300.0         !  
        real(8) :: par_R_LEFT = -240.0 ! -390.0         !  ����� �������
        real(8) :: par_R_inner = 9.0! 5.0_8     ! �� ������ ���������� ���������� �����

        !! ���������� ��������� ---------------------------------------
        real(8) :: par_a_2 = 0.130735_8        ! �������� � ������� �����������  !! ����� ����� ���������, ���� ����������
        real(8) :: par_ggg = 5.0/3.0                 ! �� ������ ���������� ���������� �����
        real(8) :: par_Velosity_inf = -2.54385_8 !-2.54278_8
        real(8) :: par_n_H_LISM = 3.0_8
        real(8) :: par_Kn = 50.3858   !0.4326569808         ! � �����������
        real(8) :: par_nu_ph = 12.0969 
        real(8) :: par_E_ph = 0.10878
        real(8) :: par_chi = 41.0391
        real(8) :: par_rho_E = 1.0                 !? ������ ��� ���������� ����� �� ������������ (�� ����������� � ����)
        real(8) :: par_Max_e = 10.0!! 5.91662
        real(8) :: par_poglosh = 0.618589! 0.389274        !! �� ��� ��������������� ����������
        real(8) :: par_rho_LISM = 1.60063        !! ����������� ���������� ��������� ��-�� ����� �� ������������� � �.�.
        real(8) :: par_p_LISM = 1.15        !! ����������� ���������� �������� ��-�� ����� �� ������������� � �.�.
        real(8) :: par_rho_He_Lism = 0.6_8         !! ����� ����� �� ����� ��������� - ��� �����
        real(8) :: par_rho_He_E = 0.155327_8       

        ! ��� ������������ �����
        real(8) :: nu_e_impact = 0.214207_8! 3.86624! 0.15465_8       
        real(8) :: lambda_e = 12.1458_8       

        !! -------------------------------------------------------------

        !! ����� ���������� ��������
        real(8) :: par_kk1 = 2.0_8     ! ������� �������� ����� � ���� � ������� �� TS: 1 - ��������, 2 - ������������ � �.�.
        real(8) :: par_kk2 = 2.08_8 !1.7_8     ! ������� �������� � �������� ������� �� �������������
        real(8) :: par_kk3 = 1.8_8     ! ������� �������� � ������
        real(8) :: par_kk31 = 1.0_8     ! ������� �������� � ������ ��� ����� �� �������� (������ ����� � � - ����)
        real(8) :: par_kk13 = 1.8_8     ! ������� �������� ����� � �������� ������� �� ������� ������� ����  �� 0 �� 1
        real(8) :: par_kk131 = 0.1_8
        real(8) :: par_kk132 = 1.5_8
        real(8) :: par_kk14 = 1.0_8     ! ������� �������� ����� � �������� ������� �� ���������� ������� ����  �� 0 �� 1
        ! (�������� ����� � TS � HP)  
        real(8) :: par_kk12 = 1.0_8     ! ������� �������� ����� �� TS � ������� �����  >= 1
        real(8) :: par_kk113 = 1.6_8     ! 1.6
        ! ������ �������� �� 4 ��� �������� ������ ����������� � ����������

        !! ������������� ���������
        real(8) :: par_nat_TS = 0.02_8 ! 0.003_8   ! ������������� ���������
        real(8) :: par_nat_HP = 0.005_8   ! ������������� ���������
        real(8) :: par_nat_BS = 0.003_8   ! 0.002 ������������� ���������

        !! �������� �������� ������������
        real(8) :: par_koeff_TS = 0.002_8 
        real(8) :: par_koeff_HP = 0.1_8   
        real(8) :: par_koeff_BS = 0.1_8   

        !! ��������� ��� �����-����� ----------------------------------------------------------

        integer(4), allocatable :: sensor(:, :, :)  !(3, 2, : par_n_potok - ����� �������)  ! ������� ��������� ����� 
	    ! ������� ������ �� ��� �������

        integer(4), allocatable :: stek(:)   ! (: ����� �������) ���������� ������ � ������ � ����
	    ! ��� ����� ����������, ��� ���-�� �����, ����� ��������, ����� ��������� �������� �� 1

        real(8), allocatable :: M_K_particle(:, :, :)   ! ������� (8, par_stek, ����� �������)
        ! (��� ����������, ��� ��������, ���, ������ ���������)
        integer(4), allocatable :: M_K_particle_2(:, :, :)  ! ������� (5, par_stek, ����� �������)
        ! (� ����� ������ �������, ����, ���� ���������� �� r, ���� ���������� �� ����, � ����� ������ � ���������������� �����)
        logical(4), allocatable :: M_K_particle_3(:, :, :, :)  ! ������� (par_n_zone + 1, par_m_zone + 1, par_stek, ����� �������)
        ! ������ ��� ������ ���������� ����� �� ����� 

        integer(4) :: par_n_moment = 13 !9  !  ������� ��������� �������� ������� (������ �������)
        real(8), allocatable :: M_K_Moment(:, :, :, :)  ! (19, par_n_sort, :, par_n_potok) ��, ��� ����������� � ������� (�� ������� ����� ��������)
        !(rho, u, v, T, In, Iu, Iv, IT, Huu, Huv, Hvv, Huuu, Hvvv)
        !(1  , 2, 3, 4, 5,  6,  7,  8,  9,   10,  11,  12,   13)

        real(8) :: MK_R_zone(par_n_zone)   ! ������� ���
        real(8) :: MK_al_zone(par_m_zone)   ! ���� ���
        real(8) :: MK_SINKR(par_m_zone + 1)   ! ����������� ������ ��� ������ ���� �� ����
        real(8), allocatable :: MK_Mu(:, :, :)   ! ���� ��� (par_n_zone + 1, par_m_zone + 1, ������ par_n_sort)
        real(8), allocatable :: MK_Mu_statistic(:, :, :)   ! ���� ��� (par_n_zone + 1, par_m_zone + 1, ������ par_n_sort)
        ! ��� ������������ ����� ���

        real(8) :: MK_gam_zone(par_n_zone)   ! �������� ����� ��� ���
        real(8) :: MK_A0_, MK_A1_   ! ��������� ��� ���������� �������

        real(8) :: sqv_1, sqv_2, sqv_3, sqv_4, sqv   ! ������ ������ ����� 
        real(8) :: MK_mu1, MK_mu2, MK_mu3, MK_mu4
        integer(4) :: MK_N                  ! ������� ����� ������ �������� (����� �� ���� �������)

        real(8) :: MK_Mu_mult = 100.0_8  ! �� ��� ��������� ���� ��� ��������� ������ ��������

        real(8) :: par_Rleft   ! ����� ������ ��� �����-����� (��� ������ ��������)
        real(8) :: par_Rup     ! ������� ������ ��� �����-����� (��� ���� ��� � �����)

        !! -------------------------------------------------------------------------------------
        

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
        !? ������������� �� ������������ ����� � ������ �� �����?

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
        real(8), allocatable :: gl_Cell_alpha_center(:)   ! (����� �����) �������� ���� ������ ������  !? �� ��������� (������ ����� ���������)
        

        integer(4), allocatable :: gl_all_Gran(:,:)       ! ��� ����� (2,:) ����� �� 2 ����
        integer(4), allocatable :: gl_Gran_neighbour(:,:) ! ������ ������ ����� (2,:) ����� �� 2 ������, ������� ���� �� ������� �� �������
        real(8), allocatable :: gl_Gran_normal(:, :, :)       ! (2, :, 2) ������� �����     !! ������� ��� �� ������ ������ �� ������ �����������!                   
        real(8), allocatable :: gl_Gran_length(:, :)       ! (:, 2) ����� �����                       
        real(8), allocatable :: gl_Gran_Center(:, :, :)       ! (2, :, 2) ����� �����       

        real(8), allocatable :: gl_Gran_POTOK(:)       ! ����� ����� ������ - �������� ��� ������ ����� ����� ���������� � ����� ������        
        
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
        !! ��������� ��� ���� ����� � ���������� �������

        integer(4), allocatable :: gl_Gran_type(:)      ! ���������� ��� ����� 
        ! (0 - �������, 1 - TS, 2 - HP, 3 - BS)     

        integer(4), allocatable :: gl_Gran_shem(:)      ! ���������� ����� ��� �����
        ! 0 - Lax
        ! 1 - HLL
        ! 2 - HLLC
        ! 3 - Godunov

        !! ������
        integer(4) :: n_Hidrogen = par_n_sort  ! ����� ������ ������ ��������
        integer(4) :: n_atom_source = 7  ! ����� ���������� ������
        integer(4) :: n_par = 5! 6  !! ����� ���������� ���������� � ������   5 - ���� ��� �����
        ! 5 ����������������
        ! 4 * n_Hidrogen - �������

        real(8), allocatable :: gd(:, :, :)  ! (n_par, :, 2 ��������� ����)
        ! (rho p u v Q He)
        !   1  2 3 4 5 6
        real(8), allocatable :: hydrogen(:, :, :, :)  ! (5, n_Hidrogen, :, 2 ��������� ����)
        ! rho p u v T

        real(8), allocatable :: atom_all_source(:, :, :)  ! (4, n_Hidrogen, : ����� �����)
        ! (In, Iu, Iv, IT)

        real(8), allocatable :: atom_source(:, :)  ! (n_atom_source, : ����� �����)
        ! (k_u, k_v, k_T, In, Iu, Iv, IT)
        ! (1     2    3    4   5  6    7)

        LOGICAL :: pogl_ = .True.  ! ������� �� ����������
        real(8) :: pogl_v_min = -15.0
        real(8) :: pogl_v_max = 15.0
        integer(4) :: pogl_iter = 300
        real(8) :: pogl_ddd
        real(8), allocatable :: pogloshenie(:, :, :)   !  (n_Hidrogen, �������� �� ��������, �����)


        !! PUI 
        ! PUI - ����������� ����, ������� ������ ���������� �� � "Geometry", � � "PUI" ������ ��� ������������� ������ � ��������
        LOGICAL :: culc_pui = .True.   ! ������� �� PUI ? 
        integer :: pui_nW = 60      ! 50   !TODO �����������
        real(8) :: pui_wR = 150.0    ! 150.0  !TODO �����������
        integer :: pui_size            ! ������� ����� �������� pui
        integer :: pui_n_par = 3            ! ������� ����� �������� pui    !TODO �����������
        real(8), allocatable :: f_pui(:, :)            ! (pui_nW, : pui_size) !TODO �����������
        integer, allocatable :: f_pui_num(:)           ! �� ������ � ������� ���, ���������� ����� ������ � �����  !TODO �����������
        integer, allocatable :: f_pui_num2(:)		   ! �� ������ ������ � �����, ���������� ����� � ������� PUI (���� �� ���� - ���� ������� ����������) !TODO �����������
        ! 0 - ������ �� �������� pui
        
        real(8), allocatable :: par_pui(:, :)               ! (3 pui_n_par, : pui_size)	   ��������� ������� !TODO �����������
        ! (n_pui, T_pui, p_pui)
        real(8), allocatable :: pui_Sm(:, :)           ! (pui_nW, : pui_size)  !TODO �����������
	    real(8), allocatable :: pui_Sp(:, :)           ! (pui_nW, : pui_size)  !TODO �����������

        !? ������� h0(U_H) - ��. ������������ PUI  !TODO ����� ������ �� �����������
        integer :: pui_h0_n = 1000      
        real(8) :: pui_h0_wc = 100.0    
        real(8), allocatable :: h0_pui(:)   ! ��� ������� ������� ��� ��������� - � �������� ��� ���������� (��. ������������)
        

        !? ��������� ������� (������� ������� ������� ���������)
        integer :: pui_F_n = 300      ! �� ������� ������ �� ��������� ������������� ��� ��������� PUI 
        ! �.�. �� ����� ����������� ksi �� 0 �� 1 � ����� �������� �������� ��� ���� ksi. 
        ! �������� �������� ��� dksi = 1.0/pui_F_n � � ��������� ����� ������� ������� ���������������
        real(8), allocatable :: F_integr_pui(:, :)           ! (pui_F_n, :) ������������� ��� ���������
        real(8), allocatable :: nu_integr_pui(:, :)           ! (pui_F_n, :)   ������� �����������
        real(8), allocatable :: Mz_integr_pui(:, :)           ! (pui_F_n, :)   �������� ��������
        real(8), allocatable :: E_integr_pui(:, :)           ! (pui_F_n, :)	   �������� �������
        integer (kind=omp_lock_kind), allocatable :: pui_lock(:)  ! ��� openMP

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
        real(8), allocatable :: gl_Cell_center(:,:)   ! (2, :) 

        logical, allocatable :: gl_all_triangle(:)
        ! .True. - ���� ������ �������� ������������� � .False. � ��������� ������

        integer(4), allocatable :: gl_Cell_neighbour(:,:)   ! (4, :) ����� �� 4 ������� ��� ������ ������ 
        real(8), allocatable :: gl_Cell_Belong(:, :, :)       ! ��� ����� (3, 4, :) ����� �� A B C, 4 ����� � ������, 
        ! ������ ����� - ��� 1 � 2 ����, ������ - 2 � 3 � �.�.
        ! Ax + By + C = 0  (���� > 0 �� �� ��������� ������)

        real(4), allocatable :: gl_Cell_interpol_matrix(:, :, :)   ! (4, 4, :) ����� �� 4 ������� ��� ������ ������
        ! Ax + By + Cxy + D
        ! Ax + By + C
        ! 
        ! ���������������� ������� ��� ������ ������
        ! ��� ������������ �� �������� ��� ������ ��������! ������ �������� �� ��������� (��������������� �������� ������ ������, ��� �������� � �����)

        !! ������
        integer(4) :: n_Hidrogen = 4  ! ����� ������ ������ ��������
        integer(4) :: n_par = 5  !! ����� ���������� ���������� � ������
        ! 5 ����������������
        ! 4 * n_Hidrogen - �������

        real(8), allocatable :: gd(:, :)  ! (n_par, :)
        ! (rho p u v Q He)
        real(8), allocatable :: hydrogen(:, :, :)  ! (5, n_Hidrogen, :)
        ! rho p u v T

        real(8), allocatable :: atom_all_source(:, :, :)  ! (4, n_Hidrogen, : ����� �����)
        ! (In, Iu, Iv, IT)

        real(8), allocatable :: atom_source(:, :)  ! (7, : ����� �����)
        ! (k_u, k_v, k_T, In, Iu, Iv, IT)

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

        ! ��������� �������� �� ������ ��� ������� PUI 
        if(ALLOCATED(SS%f_pui) == .False.) then 
            ! ���������, ������� ����� ����� ��������� pui
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

            ! ������ ���� ���������� ����� ����� �������� ������ ��� � �������� � �����
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


