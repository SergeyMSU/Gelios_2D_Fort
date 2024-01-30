module Monte_Karlo
    use STORAGE 
    use GEOMETRY
    USE ieee_arithmetic
    USE My_func
	USE OMP_LIB
    implicit none 


    contains

    subroutine M_K_start(SS, SS_int)
        ! Variables
        TYPE (Setka), intent(in out) :: SS
        TYPE (Inter_Setka), intent(in) :: SS_int
        real(8) :: start_time, end_time
        integer(4) :: step, iter, potok
        real(8) :: mu_(par_n_zone + 1), Wt_(par_n_zone + 1), Wp_(par_n_zone + 1), Wr_(par_n_zone + 1), X_(par_n_zone + 1)
		logical :: bb
        integer(4) :: i, cell
        real(8) :: sin_, x, phi, y, z, ksi, Vx, Vy, Vz, r_peregel, no, ksi1, ksi2, ksi3, ksi4, ksi5

		integer(4) :: num, to_i, to_j, j, pp, k
		real(8) :: ll, rr, Vphi, Vr, pui_w2, pui_w1
		real(8), allocatable :: vol_sr(:)                                    ! ��� ���������� � �����
		real(8), allocatable :: M_K_Moment_print(:, :, :)
		real(8) :: PAR(9), MAS_PUI(2)
		logical :: inzone
		integer :: nk

        call M_K_Set(SS)    ! ������� �������
        call M_K_init(SS)   ! �������������� ���� � �.�.

        end_time = 0.0
		start_time = 0.0
        call omp_set_num_threads(par_n_potok)
		start_time = omp_get_wtime()
        step = 1
        cell = 1


        call Get_sensor_sdvig(SS, 0)

		!$omp parallel
		!$omp do private(potok, num, mu_, Wt_, Wp_, Wr_, X_, bb, i, Vx, Vy, &
		!$omp Vz, cell, sin_, x, phi, y, z, ksi, r_peregel, no, to_i, to_j, ksi1, &
		!$omp ksi2, ksi3, ksi4, ksi5, ll, rr, Vphi, Vr, inzone)
        do iter = 1, par_n_potok * par_n_parallel

            potok = (omp_get_thread_num() + 1) 
			cell = 1

            !$omp critical
                print*, "start potok = ", potok, " iter = ", iter, "   step = ", step, "from = ", par_n_potok * par_n_parallel
                step = step + 1
			!$omp end critical

            ! ��������� ������� ������� ���� (� ���������)
            do num = 1, MK_N1
                call MK_Init_Parametrs(SS, potok, mu_, Wt_, Wp_, Wr_, X_, bb)

                do i = 1, par_n_zone + 1
                    sin_ = sqrt(1.0 - (X_(i)**2))
					x = (par_Rmax) * X_(i)
					call M_K_rand(SS%sensor(1, 1, potok), SS%sensor(2, 1, potok), SS%sensor(3, 1, potok), ksi)
					phi = 2.0 * par_pi * ksi
					y = (par_Rmax) * sin_ * cos(phi)
					z = (par_Rmax) * sin_ * sin(phi)
					
                    call Geo_Find_Cell(SS, x, sqrt(y**2 + z**2), cell, inzone = inzone)
					if(inzone == .False.) then
						print*, x, sqrt(y**2 + z**2), cell
						STOP "ERROR 98y7gwfubevevgte"
					end if
					call dekard_skorost(x, y, z, Wr_(i), Wp_(i), Wt_(i), Vx, Vy, Vz)
				
					if(cell < 1) then
						STOP "Error cell < 1 MK 0lhy976yihko"
					end if
					
					
					if(i /= par_n_zone + 1 .or. bb == .True.) then
						! ��������� ������� � ����
						SS%stek(potok) = SS%stek(potok) + 1
						SS%M_K_particle(1:7, SS%stek(potok), potok) = (/ x, y, z, Vx, Vy, Vz, mu_(i) * SS%MK_mu1 * SS%MK_Mu_mult /)
						SS%M_K_particle_2(1, SS%stek(potok), potok) = cell       ! � ����� ������ ���������
						SS%M_K_particle_2(2, SS%stek(potok), potok) = SS%gl_all_Cell_zone(cell) ! ����
						call MK_Distination(SS, SS%M_K_particle(1:3, SS%stek(potok), potok), SS%M_K_particle(4:6, SS%stek(potok), potok),&
							to_i, to_j, r_peregel)
						SS%M_K_particle(8, SS%stek(potok), potok) = r_peregel
						SS%M_K_particle_2(3, SS%stek(potok), potok) = to_i  ! ���� ����������
						SS%M_K_particle_2(4, SS%stek(potok), potok) = to_j  ! ���� ����������
					end if

                    call M_K_Fly(SS, SS_int, potok)

                end do


            end do

			! ��������� ������� ������� ���� (����� ������)
			do num = 1, MK_N2
				call M_K_rand(SS%sensor(1, 1, potok), SS%sensor(2, 1, potok), SS%sensor(3, 1, potok), ksi1)
				call M_K_rand(SS%sensor(1, 1, potok), SS%sensor(2, 1, potok), SS%sensor(3, 1, potok), ksi2)
				call M_K_rand(SS%sensor(1, 1, potok), SS%sensor(2, 1, potok), SS%sensor(3, 1, potok), ksi3)
				call M_K_rand(SS%sensor(1, 1, potok), SS%sensor(2, 1, potok), SS%sensor(3, 1, potok), ksi4)
				call M_K_rand(SS%sensor(1, 1, potok), SS%sensor(2, 1, potok), SS%sensor(3, 1, potok), ksi5)
				
				ll = SS%par_Rleft
				rr = -0.001;
				x = ll + ksi1 * (rr - ll)
				phi = ksi2 * 2.0 * par_pi
				Vphi = cos(2.0 * par_pi * ksi3) * sqrt(-log(1.0 - ksi4))
				Vx = SS%par_Velosity_inf + sin(2.0 * par_pi * ksi3) * sqrt(-log(1.0 - ksi4))
				Vr = -sqrt(-log(ksi5))
				y = SS%par_Rup * cos(phi)
				z = SS%par_Rup * sin(phi)
				
				
				call Geo_Find_Cell(SS, x, sqrt(y**2 + z**2), cell, inzone = inzone)
				
				if(cell < 1) then
					print*, x, y, z, cell
					!$MPI call MPI_Abort(MPI_COMM_WORLD, 104, mpi_ierror)
					STOP "Error 0lhy976yihkoqwewqdqwd  "
				end if
				
				
				SS%stek(potok) = SS%stek(potok) + 1
				SS%M_K_particle(1:7, SS%stek(potok), potok) = (/ x, y, z, Vx, cos(phi) * Vr - sin(phi) * Vphi,&
					sin(phi) * Vr + cos(phi) * Vphi,  SS%MK_mu2 * SS%MK_Mu_mult /)
				SS%M_K_particle_2(1, SS%stek(potok), potok) = cell       ! � ����� ������ ���������
				SS%M_K_particle_2(2, SS%stek(potok), potok) = SS%gl_all_Cell_zone(cell) ! ����
				call MK_Distination(SS, SS%M_K_particle(1:3, SS%stek(potok), potok), SS%M_K_particle(4:6, SS%stek(potok), potok),&
					to_i, to_j, r_peregel)
				SS%M_K_particle(8, SS%stek(potok), potok) = r_peregel
				SS%M_K_particle_2(3, SS%stek(potok), potok) = to_i  ! ���� ����������
				SS%M_K_particle_2(4, SS%stek(potok), potok) = to_j  ! ���� ����������

				call M_K_Fly(SS, SS_int, potok)

			end do

			! ��������� ������� �������� ���� (����� �����)
			do num = 1, MK_N3
				call MK_Velosity_initial2(SS, potok, Vx, Vy, Vz)
				call M_K_rand(SS%sensor(1, 1, potok), SS%sensor(2, 1, potok), SS%sensor(3, 1, potok), ksi1)
				call M_K_rand(SS%sensor(1, 1, potok), SS%sensor(2, 1, potok), SS%sensor(3, 1, potok), ksi2)
				rr = sqrt(ksi1 * (SS%par_Rup) * (SS%par_Rup))
				phi = ksi2 * 2.0 * par_pi
				y = rr * cos(phi)
				z = rr * sin(phi)
				
				call Geo_Find_Cell(SS, SS%par_Rleft, sqrt(y**2 + z**2), cell, inzone = inzone)
				
				if(cell < 1) then
					!$MPI call MPI_Abort(MPI_COMM_WORLD, 105, mpi_ierror)
					STOP "Error 0lhy976yihkodfresdfre  "
				end if
				
				SS%stek(potok) = SS%stek(potok) + 1
				SS%M_K_particle(1:7, SS%stek(potok), potok) = (/ SS%par_Rleft, y, z, Vx, Vy, Vz,  SS%MK_mu3 * SS%MK_Mu_mult /)
				SS%M_K_particle_2(1, SS%stek(potok), potok) = cell       ! � ����� ������ ���������
				SS%M_K_particle_2(2, SS%stek(potok), potok) = SS%gl_all_Cell_zone(cell) ! ����
				call MK_Distination(SS, SS%M_K_particle(1:3, SS%stek(potok), potok), SS%M_K_particle(4:6, SS%stek(potok), potok),&
					to_i, to_j, r_peregel)
				SS%M_K_particle(8, SS%stek(potok), potok) = r_peregel
				SS%M_K_particle_2(3, SS%stek(potok), potok) = to_i  ! ���� ����������
				SS%M_K_particle_2(4, SS%stek(potok), potok) = to_j  ! ���� ����������

				call M_K_Fly(SS, SS_int, potok)
			end do
			
			! ��������� ������� ��������� ���� (����� ������� � ����� ���������)
			do num = 1, MK_N4
				call MK_Velosity_initial(SS, potok, Vx, Vy, Vz)
				call M_K_rand(SS%sensor(1, 1, potok), SS%sensor(2, 1, potok), SS%sensor(3, 1, potok), ksi1)
				call M_K_rand(SS%sensor(1, 1, potok), SS%sensor(2, 1, potok), SS%sensor(3, 1, potok), ksi2)
				rr = sqrt(ksi1 * ((SS%par_Rup**2) - (par_Rmax**2)) + (par_Rmax**2))
				phi = ksi2 * 2.0 * par_pi
				y = rr * cos(phi)
				z = rr * sin(phi)
				
				call Geo_Find_Cell(SS, -0.001_8, sqrt(y**2 + z**2), cell, inzone = inzone)
				
				if(cell < 1) then
					!$MPI call MPI_Abort(MPI_COMM_WORLD, 106, mpi_ierror)
					STOP "Error 0lhy976yihko133131312  "
				end if
				
				SS%stek(potok) = SS%stek(potok) + 1
				SS%M_K_particle(1:7, SS%stek(potok), potok) = (/ -0.001_8, y, z, Vx, Vy, Vz,  SS%MK_mu4 * SS%MK_Mu_mult /)
				SS%M_K_particle_2(1, SS%stek(potok), potok) = cell       ! � ����� ������ ���������
				SS%M_K_particle_2(2, SS%stek(potok), potok) = SS%gl_all_Cell_zone(cell) ! ����
				call MK_Distination(SS, SS%M_K_particle(1:3, SS%stek(potok), potok), SS%M_K_particle(4:6, SS%stek(potok), potok),&
					to_i, to_j, r_peregel)
				SS%M_K_particle(8, SS%stek(potok), potok) = r_peregel
				SS%M_K_particle_2(3, SS%stek(potok), potok) = to_i  ! ���� ����������
				SS%M_K_particle_2(4, SS%stek(potok), potok) = to_j  ! ���� ����������

				call M_K_Fly(SS, SS_int, potok)
			end do

        end do
		!$omp end do

		!$omp end parallel

		print*, "END ALL particle"
		
		end_time = omp_get_wtime()
		print *, "Time work: ", (end_time-start_time)/60.0, "   in minutes"

		no = SS%MK_Mu_mult * SS%MK_N * par_n_claster
		SS%M_K_Moment(:, :, :, :) = SS%M_K_Moment(:, :, :, :) / no  ! ����� ���� ��� ��������� ������ �������� ��� ��������

		do i = 2, par_n_potok
			SS%M_K_Moment(:, :, :, 1) = SS%M_K_Moment(:, :, :, 1) + SS%M_K_Moment(:, :, :, i)
		end do

		! ����� �� ���� ������� � ��������� �������
		do i = 1, size(SS%M_K_Moment(1, 1, :, 1))

			no = Geo_Get_Volume_Rotate(SS, i, 360.0_8)
			
			SS%M_K_Moment(:, :, i, 1) = SS%sqv * SS%M_K_Moment(:, :, i, 1) / no

			do j = 1, par_n_sort
				if(SS%M_K_Moment(1, j, i, 1) > 0.000001) then
					SS%M_K_Moment(2:3, j, i, 1) = SS%M_K_Moment(2:3, j, i, 1)/SS%M_K_Moment(1, j, i, 1)  ! ��������
					SS%M_K_Moment(4, j, i, 1) = (1.0/3.0) * ( SS%M_K_Moment(4, j, i, 1)/SS%M_K_Moment(1, j, i, 1) - &
						kvv(SS%M_K_Moment(2, j, i, 1), SS%M_K_Moment(3, j, i, 1), 0.0_8) )  ! Temp
					
					SS%M_K_Moment(9, j, i, 1) = SS%M_K_Moment(9, j, i, 1) / SS%M_K_Moment(1, j, i, 1) - &
						SS%M_K_Moment(2, j, i, 1)**2
					SS%M_K_Moment(10, j, i, 1) = SS%M_K_Moment(10, j, i, 1) / SS%M_K_Moment(1, j, i, 1) - &
						SS%M_K_Moment(2, j, i, 1)*SS%M_K_Moment(3, j, i, 1)
					SS%M_K_Moment(11, j, i, 1) = SS%M_K_Moment(11, j, i, 1) / SS%M_K_Moment(1, j, i, 1) - &
						SS%M_K_Moment(3, j, i, 1)**2
					SS%M_K_Moment(12, j, i, 1) = 2.0 * SS%M_K_Moment(2, j, i, 1)**3 + 3.0 * SS%M_K_Moment(2, j, i, 1) * SS%M_K_Moment(9, j, i, 1) - &
						SS%M_K_Moment(12, j, i, 1) / SS%M_K_Moment(1, j, i, 1) 
					SS%M_K_Moment(13, j, i, 1) = 2.0 * SS%M_K_Moment(3, j, i, 1)**3 + 3.0 * SS%M_K_Moment(3, j, i, 1) * SS%M_K_Moment(11, j, i, 1) - &
						SS%M_K_Moment(13, j, i, 1) / SS%M_K_Moment(1, j, i, 1) 	
				end if
			end do
			
			SS%M_K_Moment(5:8, :, i, 1) = SS%M_K_Moment(5:8, :, i, 1) * SS%par_n_H_LISM
		end do

		! ��� ��������� ������� ����� ��������� �������� � ����������
		loop2: do i = 1, size(SS%M_K_Moment(1, 1, :, 1))
			do j = 1, 4
				pp = SS%gl_all_Cell(j, i)
				if(pp < 1) CYCLE
				
				if((SS%gl_yzel(1, pp, 1) >= 0.0 .and. norm2(SS%gl_yzel(:, pp, 1)) >= par_Rmax)) then
					SS%M_K_Moment(:, :, i, 1) = 0.0
					SS%M_K_Moment(1, 4, i, 1) = 1.0
					SS%M_K_Moment(5, 4, i, 1) = 1.0
					SS%M_K_Moment(10, 4, i, 1) = 0.5
					SS%M_K_Moment(13, 4, i, 1) = 0.5
					SS%M_K_Moment(15, 4, i, 1) = 0.5
					SS%M_K_Moment(2, 4, i, 1) = SS%par_Velosity_inf
					CYCLE loop2
				end if
			end do
		end do loop2



		! ��������� ��������� �������� � ������� ����������

		if(par_n_sort /= SS%n_Hidrogen) print*, "ERROR 9iur987e8y7bn24243289vwrtbryd"

		SS%hydrogen(1, :, :, 1) = SS%M_K_Moment(1, :, :, 1)
		SS%hydrogen(3:5, :, :, 1) = SS%M_K_Moment(2:4, :, :, 1)
		SS%hydrogen(:, :, :, 2) = SS%hydrogen(:, :, :, 1)

    end subroutine M_K_start

    subroutine M_K_Fly(SS, SS_int, n_potok)
        TYPE (Setka), intent(in out) :: SS
        TYPE (Inter_Setka), intent(in) :: SS_int
        integer(4), intent(in) :: n_potok

		real(8) :: particle(8)
		integer(4):: particle_2(5), i, ijk
		logical :: particle_3(par_n_zone + 1, par_m_zone + 1)
		
		integer(4) :: num  ! ����� �������, ������� � �����
		integer(4) :: cell ! ����� ������, � ������� ��������� �������
		integer(4) :: next ! ����� ������, � ������� ������ ������� � ��������� ���
		integer(4) :: next2 ! ����� ������, � ������� ������ ������� � ��������� ���
		integer(4) :: area2  ! ����, � ������� ������ ��������� ������
		integer(4) :: II  ! �� ������� ������ ������������ ���� ��� �����������
		integer(4) :: to_i, to_j, from_i, from_j
		logical :: bb, bb2, inzone
		
		real(8) :: time ! ��������� ����� �� ������ ������� �� ������
		
		real(8) :: cp, vx, vy, vz, ro, PAR(5)  ! ��������� ������ � ������
		real(8) :: uz, nu_ex, kappa, ksi, t_ex, t2, mu_ex, mu2, r_ex(3), r, mu, u, V(3), mu3, r_exit(3)
		real(8) :: uz_M, uz_E, k1, k2, k3, u1, u2, u3, skalar
		real(8) :: Ur, Uthe, Uphi, Vr, Vthe, Vphi, phi
		real(8) :: v1, v2, v3, r_peregel, ddt
		
		real(8) :: nu_ph, kappa_ph, kappa_all, mu_ph, mu_perez
		
		real(8) :: Wr(par_n_zone + 1), Wthe(par_n_zone + 1), Wphi(par_n_zone + 1), mu_(par_n_zone + 1)
		
		integer(4) :: step
	
		step = 0
		
		do while (SS%stek(n_potok) >= 1)
			step = step + 1
			
			
			if(SS%stek(n_potok) > par_stek * 0.9) then
				print*, "1234543fj976r  Perepolnen stek", SS%stek(n_potok) 
				pause
				STOP
			end if
			
				
			num = SS%stek(n_potok)
			SS%stek(n_potok) = SS%stek(n_potok) - 1
			! ���� ��� ��������� �������
			particle = SS%M_K_particle(:, num, n_potok)
			particle_2 = SS%M_K_particle_2(:, num, n_potok)
			if(MK_Mu_stat) particle_3 = SS%M_K_particle_3(:, :, num, n_potok)

			
			
			
			loop1: do  ! ���� ������� �� ������� �� ������� ��� �� ��������� ��������
				!print*, particle(1), norm2(particle(2:3))
				!print*, particle
				!print*, "__________________"
				cell = particle_2(1)
				mu = particle(7)                                ! ��� �������
			
				!print*, "AA"
				call Geo_Time_fly(SS, particle(1:3), particle(4:6), time, cell, next)  ! ������� ����� time �� ������ �� ������
				!print*, "BB"
				time = max(0.0000001_8, time * 1.0001) ! �������� �����, ����� ������� ����� ����� �� ������

				r_exit = particle(1:3) + time * particle(4:6)
				next2 = cell
				call Geo_Find_Cell(SS, r_exit(1), sqrt(r_exit(2)**2 + r_exit(3)**2), next2, inzone)
				do while(next2 == cell .and. inzone == .True.) 
					! print*, "Error 171 next2 == cell  time = ", time
					! print*, particle(1), norm2(particle(2:3))
					! print*, "----------------------"
					! print*, particle(4:6)
					! print*, "----------------------"
					! print*, r_exit(1), norm2(r_exit(2:3))
					! print*, "----------------------"
					! STOP
					time = time * 1.05
					r_exit = particle(1:3) + time * particle(4:6)
					call Geo_Find_Cell(SS, r_exit(1), sqrt(r_exit(2)**2 + r_exit(3)**2), next2, inzone)
				end do
				!print*, "CC"
				
				kappa = 0.0
				do ijk = 1, 3
					
					select case (ijk)
						case(1)
							ddt = 1.0/6.0                                      
						case(2)
							ddt = 1.0/2.0                                    
						case(3)
							ddt = 5.0/6.0                                  
						case default
							print*, "Error uijkhjgfbnbnn hbuhuefw"
							STOP
					end select

					phi = polar_angle(particle(2) + ddt * time * particle(5), particle(3) + ddt * time * particle(6))
					
					PAR = SS%gd(:, cell, 1)  !! ���� ��� ������������
				
					cp = sqrt(PAR(2)/PAR(1))
					vx = PAR(3)
					vy = PAR(4) * cos(phi)
					vz = PAR(4) * sin(phi)
					ro = PAR(1)

					if(ieee_is_nan(cp)) then
                    print*, "error cp nan bntjumki,lo;.piuyntbvtcrwrwfv"
                    print*, "cell = ", cell
                    print*, "______________ 0"
                    print*, PAR
                    print*, "______________ 1"
                    print*, particle
                    print*, "______________ 2"
                    print*, ddt, time
                    print*, "______________ 3"
					pause
                    STOP
                end if
				
					if(ro <= 0.0 .or. ro > 1000.0) then
						print*, PAR
						print*, "___"
						print*, cell
						pause "ERROR ro MK 157 6787yutr4dfghhghjuhj0089"
					end if
				
					! ����� ����� �� ����������� � ���� ������  ****************************************************************************************
					u = sqrt(kvv(particle(4) - vx, particle(5) - vy, particle(6) - vz))
					u1 =  vx - particle(4)
					u2 =  vy - particle(5)
					u3 =  vz - particle(6)
					skalar = particle(4) * u1 + particle(5) * u2 + particle(6) * u3
				
					if (u / cp > 7.0) then
						uz = MK_Velosity_1(u, cp);
						nu_ex = (ro * uz * MK_sigma(uz)) / SS%par_Kn
					else
						nu_ex = (ro * MK_int_1(u, cp)) / SS%par_Kn        ! ������� ��������� ���������� ��������
					end if
			
					kappa = kappa + (nu_ex * time/3.0)  ! �� �����������
				end do
				!print*, "DD"
				area2 = SS%gl_all_Cell_zone(cell) ! ���� ��������
				
				r = norm2(particle(1:3) + time/2.0 * particle(4:6))
				
				if(MK_photoionization) then
					nu_ph = SS%par_nu_ph * (SS%par_R0/r)**2
					kappa_ph = (nu_ph * time)     ! �� �������������
				end if
				
				kappa_all = kappa
				if(MK_photoionization) kappa_all = kappa_all + kappa_ph
				
				call M_K_rand(SS%sensor(1, 2, n_potok), SS%sensor(2, 2, n_potok), SS%sensor(3, 2, n_potok), ksi)
				
				t_ex = -(time / kappa_all) * log(1.0 - ksi * (1.0 - exp(-kappa_all)))  ! ����� �� �����������
				t2 = time - t_ex  ! ����� ������� ������ ����� ����, ��� ���� �������������
				mu_perez = mu * (1.0 - exp(-kappa_all)) ! ��� ��������������� ����� �� ���� ���������
				mu2 = mu * exp(-kappa_all)  ! ��� ����������� ����������������� �����
				mu_ph = 0.0
				if(MK_photoionization) mu_ph = (kappa_ph / kappa_all) * mu_perez  ! ��� ��������������� �����
				mu_ex = mu_perez - mu_ph      ! ��� ��������������� �� �������� �����
				!print*, "EE"
				
				if(mu2 < 0.0) then
					print*, mu2, mu, mu_ex, kappa
					STOP "Eror oiuyyuiojhu987uio9i"
				end if
				

				r_ex = particle(1:3) + t_ex * particle(4:6)   ! ���������� �����������
				
				!! ��������, ���� ����������� ��������� �� ��������� ������ (����� ������� �������� � � ���� ������)

				next2 = cell
				call Geo_Find_Cell(SS, r_ex(1), sqrt(r_ex(2)**2 + r_ex(3)**2), next2) 

				if (cell /= next2) then
					t_ex = time * 0.998
					t2 = time - t_ex
					r_ex = particle(1:3) + t_ex * particle(4:6)  
				end if
				
				r = norm2(r_ex)  ! ���������� �� ����� ����������� �� ������
				
				from_i = MK_geo_zones(SS, r, 1.0_8)     ! ���� �� r � ����� �����������
				from_j = MK_alpha_zones(SS, polar_angle( r_ex(1), sqrt(r_ex(2)**2 + r_ex(3)**2) ) ) ! ���� �� ���� � ����� �����������
				
				!print*, from_i, from_j, r_ex, polar_angle( r_ex(1), sqrt(r_ex(2)**2 + r_ex(3)**2) )
				!pause
				
				if(MK_Mu_stat .and. particle_3(from_i, from_j) == .False.) then
					particle_3(from_i, from_j) = .True.
					
					!$omp critical
					SS%MK_Mu_statistic(from_i, from_j, particle_2(2)) = SS%MK_Mu_statistic(from_i, from_j, particle_2(2)) + &
						mu/max(0.3 * SS%MK_SINKR(from_j), sin(polar_angle( r_ex(1), sqrt(r_ex(2)**2 + r_ex(3)**2) )))
					!$omp end critical
				end if

				!print*, "FF"
				! ���� �������� �������
				call MK_ADD_MOMENT(SS, n_potok, particle_2(2), cell, t_ex, t2, mu, mu2, mu_ex, mu_ph, &
					 particle(4:6), particle(1:3), u, cp, uz, u1, u2, u3, skalar)
				!print*, "GG"
				
				call spherical_skorost(r_ex(1), r_ex(2), r_ex(3), vx, vy, vz, Ur, Uphi, Uthe)
				call spherical_skorost(r_ex(1), r_ex(2), r_ex(3), particle(4), particle(5), particle(6), Vr, Vphi, Vthe)
			
				! ������������
				if (area2 == 1 .or. Ur / cp > 1.8) then   ! ��� ��������������� �����������
					II = 0  ! II - ��� ������� �������������� ������ ����������� (������ ���������)
				else
					II = MK_geo_zones(SS, r, 1.2_8) - 1
				end if
				
				!print*, "GGG"
				call M_K_Change_Velosity4(SS, n_potok, Ur/cp, Uthe/cp, Uphi/cp, Vr/cp, Vthe/cp, Vphi/cp, Wr, Wthe, Wphi, mu_, &
					cp, r, II, r_ex(1), r_ex(2), r_ex(3), bb)
				!print*, "GGGG"
				
				Wr = Wr * cp
				Wthe = Wthe * cp
				Wphi = Wphi * cp
				
				!print*, "A"
				do i = 1, II + 1
					if(i == II + 1 .and. bb == .False.) CYCLE  ! ���� �� ����������� �������� ����
					call dekard_skorost(r_ex(1), r_ex(2), r_ex(3), Wr(i), Wphi(i), Wthe(i), v1, v2, v3)
					V = (/ v1, v2, v3 /)
					
					call MK_Distination(SS, r_ex, V, to_i, to_j, r_peregel)  ! ������� ���� ���������� � ����� ��������� ��� �������
					
					mu3 = mu_ex * mu_(i)
					call MK_ruletka(SS, n_potok, to_i, to_j, from_i, from_j, area2, r, r_peregel, mu3, bb2)
					
					if(bb2 == .False.) CYCLE  ! �� ��������� ��� �������, ��� ����������
					
					! ������� ����� ������� � ����
					SS%stek(n_potok) = SS%stek(n_potok) + 1
					SS%M_K_particle(1:3, SS%stek(n_potok), n_potok) = r_ex
					SS%M_K_particle(4:6, SS%stek(n_potok), n_potok) = V
					SS%M_K_particle(7, SS%stek(n_potok), n_potok) = mu3
					SS%M_K_particle(8, SS%stek(n_potok), n_potok) = r_peregel

					SS%M_K_particle_2(:,SS%stek(n_potok), n_potok) = (/ cell, area2, to_i, to_j /)
					
					if(MK_Mu_stat == .True.) then
						if(particle_2(2) == area2) then
							SS%M_K_particle_3(:, :, SS%stek(n_potok), n_potok) = particle_3
						else 
							SS%M_K_particle_3(:, :, SS%stek(n_potok), n_potok) = .False.
						end if
					end if
					
				end do
				!print*, "B"
				! ���������, ����� �� �������� ���������������� ����� �����
				call MK_ruletka(SS, n_potok, particle_2(3), particle_2(4), from_i, from_j, particle_2(2), &
						r, particle(8), mu2, bb2)
				if(bb2 == .False.) EXIT loop1  ! �� ��������� ��� �������, ��� ����������
				particle(7) = mu2
			
				! ������� ��������� ������
				
				particle(1:3) = particle(1:3) + time * particle(4:6)
				if(next <= 0 .or. inzone == .False.) then
					if (norm2(particle(1:3)) <= 100.0) then
						print*, "Error 23e2323 ", particle(1:3)
						pause
					end if
					EXIT loop1  ! ������� �������� �� ���� �������
				end if
				
				11	continue 
				
				!print*, "C"
				call Geo_Find_Cell(SS, particle(1), sqrt(particle(2)**2 + particle(3)**2), next, inzone = inzone) ! ������� ����� ��������� �������� 
				!print*, "D"
				
				if(inzone == .False.) then
					print*, particle
					print*, "____________________"
					print*, cell, next, inzone
					STOP "Error inzone 09247trfyvwuueicvtghrbertvbyj kjbhvgcfxf"
				end if

				if (next < 1) then ! �����
					next = 3
					particle(1:3) = particle(1:3) + (time/100.0) * particle(4:6)
					print*, "Zatik"
					pause
					GO TO 11
				end if
				
				particle_2(1) = next
				if(next == 0) then
					if (norm2(particle(1:3)) <= 100.0) then
						print*, "Error wqer2r42  ", particle(1:3)
						pause
					end if
					
					EXIT loop1  ! ������� �������� �� ���� �������
				end if
				
				if(particle(1) > 0.00001 .and. norm2(particle(1:3)) > par_Rmax + 0.001) then
					EXIT loop1  ! ������� �������� �� ���� �������
				end if
				
			
				
			
			end do loop1
			
			
		end do

	end subroutine M_K_Fly

	subroutine MK_ADD_MOMENT(SS, n_potok, sort, cell, t_ex, t2, mu, mu2, mu_ex, mu_ph, VV, XX, u, cp, uz, u1, u2, u3, skalar)
		TYPE (Setka), intent(in out) :: SS
		integer(4), intent(in) :: n_potok, sort, cell
		real(8), intent(in) :: t_ex, t2, mu, mu2, VV(3), XX(3), u, cp, uz, u1, u2, u3, mu_ex, mu_ph, skalar

		real(8) :: alpha, v, uz_M, uz_E, k1, k2, k3

		alpha = polar_angle(XX(2), XX(3))
		v = VV(2) * cos(alpha) + VV(3) * sin(alpha)


		SS%M_K_Moment(1, sort, cell, n_potok) = SS%M_K_Moment(1, sort, cell, n_potok) + t_ex * mu + t2 * mu2
		SS%M_K_Moment(2, sort, cell, n_potok) = SS%M_K_Moment(2, sort, cell, n_potok) + (t_ex * mu + t2 * mu2) * VV(1)
		SS%M_K_Moment(3, sort, cell, n_potok) = SS%M_K_Moment(3, sort, cell, n_potok) + (t_ex * mu + t2 * mu2) * v
		SS%M_K_Moment(4, sort, cell, n_potok) = SS%M_K_Moment(4, sort, cell, n_potok) + &
			(t_ex * mu + t2 * mu2) * kvv(VV(1), VV(2), VV(3))

		SS%M_K_Moment(9, sort, cell, n_potok) = SS%M_K_Moment(9, sort, cell, n_potok) + &
			(t_ex * mu + t2 * mu2) * VV(1)**2
		SS%M_K_Moment(10, sort, cell, n_potok) = SS%M_K_Moment(10, sort, cell, n_potok) + &
			(t_ex * mu + t2 * mu2) * VV(1) * v
		SS%M_K_Moment(11, sort, cell, n_potok) = SS%M_K_Moment(11, sort, cell, n_potok) + &
			(t_ex * mu + t2 * mu2) * v * v
		SS%M_K_Moment(12, sort, cell, n_potok) = SS%M_K_Moment(12, sort, cell, n_potok) + &
			(t_ex * mu + t2 * mu2) * VV(1)**3
		SS%M_K_Moment(13, sort, cell, n_potok) = SS%M_K_Moment(13, sort, cell, n_potok) + &
			(t_ex * mu + t2 * mu2) * v**3

		if (u / cp > 7.0) then
			uz_M = MK_Velosity_2(u, cp)/ (uz * (cp**2) * cp * par_pi * par_sqrtpi)
			uz_E = MK_Velosity_3(u, cp)
			
			SS%M_K_Moment(6, sort, cell, n_potok) = SS%M_K_Moment(6, sort, cell, n_potok) - mu_ex * uz_M * u1 / u
			SS%M_K_Moment(7, sort, cell, n_potok) = SS%M_K_Moment(7, sort, cell, n_potok) - mu_ex * uz_M * (u2 * cos(alpha) + u3 * sin(alpha)) / u
			SS%M_K_Moment(8, sort, cell, n_potok) = SS%M_K_Moment(8, sort, cell, n_potok) + &
				mu_ex * (-0.25 * (3.0 * cp**2 + 2.0 * u**2) * (uz_E / uz) - uz_M * skalar / u)
		else
			k1 = MK_int_1(u, cp)
			k2 = MK_int_2(u, cp)
			k3 = MK_int_3(u, cp)
			
			SS%M_K_Moment(6, sort, cell, n_potok) = SS%M_K_Moment(6, sort, cell, n_potok) + mu_ex * (k2/k1) * u1 / u
			SS%M_K_Moment(7, sort, cell, n_potok) = SS%M_K_Moment(7, sort, cell, n_potok) + mu_ex * (k2/k1) * (u2 * cos(alpha) + u3 * sin(alpha)) / u
			SS%M_K_Moment(8, sort, cell, n_potok) = SS%M_K_Moment(8, sort, cell, n_potok) + &
				mu_ex * (-0.5 * k3/k1 + k2/k1 * skalar / u)
		end if

		! ������� ���������� �� �������������  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		SS%M_K_Moment(5, sort, cell, n_potok) = SS%M_K_Moment(5, sort, cell, n_potok) + mu_ph 
		SS%M_K_Moment(6, sort, cell, n_potok) = SS%M_K_Moment(6, sort, cell, n_potok) + mu_ph * VV(1)
		SS%M_K_Moment(7, sort, cell, n_potok) = SS%M_K_Moment(7, sort, cell, n_potok) + mu_ph * v
		SS%M_K_Moment(8, sort, cell, n_potok) = SS%M_K_Moment(8, sort, cell, n_potok) + &
			mu_ph * (0.5 * norm2(VV) + SS%par_E_ph)

	end subroutine MK_ADD_MOMENT


    subroutine M_K_Set(SS)
		TYPE (Setka), intent(in out) :: SS
		integer(4) :: i, j, n, k
		real(8) :: Yr
		logical :: exists
		
		allocate(SS%M_K_particle(8, par_stek, par_n_potok))
		allocate(SS%M_K_particle_2(5, par_stek, par_n_potok))
		
		if(MK_Mu_stat) then
			allocate(SS%M_K_particle_3(par_n_zone + 1, par_m_zone + 1, par_stek, par_n_potok))
			allocate(SS%MK_Mu_statistic(par_n_zone + 1, par_m_zone + 1, par_n_sort))
			SS%M_K_particle_3 = .False.
			SS%MK_Mu_statistic = 0.0
		end if
		
		allocate(SS%sensor(3, 2, par_n_potok))
		allocate(SS%stek(par_n_potok))
		allocate(SS%MK_Mu(par_n_zone + 1, par_m_zone + 1, par_n_sort))
		allocate(SS%M_K_Moment(SS%par_n_moment, par_n_sort, size(SS%gl_all_Cell(1, :)), par_n_potok))
		
	
		SS%MK_Mu = 1.0
		SS%M_K_particle = 0.0
		SS%M_K_particle_2 = 0
		SS%stek = 0
		SS%sensor = 1
		SS%M_K_Moment = 0.0
		
		! ����� ������� ���
		!MK_R_zone(1) = 1.0
		SS%MK_R_zone(1) = 2.0
		SS%MK_R_zone(2) = 4.6
		SS%MK_R_zone(3) = 10.0
		SS%MK_R_zone(4) = 25.0
		SS%MK_R_zone(5) = 60.0
		SS%MK_R_zone(6) = 130.0

		
		! ����� ���� ���
		do i = 1, par_m_zone
			SS%MK_al_zone(i) = i * par_pi/(par_m_zone + 1)
		end do
		
		! ����� ����������� ����
		
		do j = 1, par_m_zone + 1
			do i = 1, par_n_zone
				SS%MK_Mu(i, j, :) = 1.0 !(MK_R_zone(i)/par_Rmax)**(1.3_8)
			end do
		end do
		
		inquire(file="MK_Mu_statistic.txt", exist=exists)
		if (exists == .False.) then
			pause "net faila!!!  cvbdfgertmkopl"
			STOP "net faila!!!"
		end if
		open(2, file = "MK_Mu_statistic.txt", status = 'old')
		
		do k = 1, par_n_sort
			do j = 1, par_m_zone + 1
				do i = 1, par_n_zone + 1
					read(2, *) n, n, n, SS%MK_Mu(i, j, k)
				end do
			end do
		end do
		
		close(2)
		
		SS%MK_Mu(:, :, 1) = SS%MK_Mu(:, :, 1) * 1.5
		SS%MK_Mu(:, :, 2) = SS%MK_Mu(:, :, 2) * 0.2
		SS%MK_Mu(:, :, 3) = SS%MK_Mu(:, :, 3) * 0.5
		SS%MK_Mu(:, :, 4) = SS%MK_Mu(:, :, 4) * 0.5
		
		
		SS%MK_SINKR(1) = sin(SS%MK_al_zone(1))
		SS%MK_SINKR(par_m_zone + 1) = sin(SS%MK_al_zone(par_m_zone))
		
		do i = 2, par_m_zone
			SS%MK_SINKR(i) = max( dabs(sin(SS%MK_al_zone(i))), dabs(sin(SS%MK_al_zone(i - 1))) )
		end do
		
		do i = 1, par_n_zone
			SS%MK_gam_zone(i) = 1.0 / ((par_Rmax / SS%MK_R_zone(i))**2 - 1.0)
			if (SS%MK_gam_zone(i) < 0.0) pause "ERROR gamma 56789oihgfr6uijt6789"
		end do
		
		Yr = dabs(SS%par_Velosity_inf)
		SS%MK_A0_ = (Yr + 1.0 / (2.0 * Yr)) * erf(Yr) + exp(-(Yr**2)) / par_sqrtpi + Yr
		SS%MK_A1_ = 1.0 + (1.0 + 1.0 / (2.0 * (Yr)**2 )) * erf(Yr) + exp(-(Yr**2) ) / (par_sqrtpi * Yr)
	end subroutine M_K_Set


    subroutine M_K_init(SS)
		TYPE (Setka), intent(in out) :: SS
		! Variables
		real(8) :: Y
		real(8) :: PAR(9) 
		integer :: cell
		
		! ������������� ��������� ����������
		SS%par_n_moment = 19
		SS%par_Rleft = SS%par_R_LEFT + 0.001
		SS%par_Rup = SS%par_R_END - 2.0
		
		Y = dabs(SS%par_Velosity_inf)
		SS%sqv_1 = (par_Rmax) * (0.5 * (par_Rmax) * par_pi * Y * &
			(erf(Y) * (1.0 + 1.0 / (2.0 * (Y**2))) + 1.0 + exp(-(Y**2)) / (Y * par_sqrtpi)))  ! ������ � ���������
		SS%sqv_2 = par_sqrtpi * (SS%par_Rup) * dabs(SS%par_Rleft) ! ������� ������
		SS%sqv_3 = par_pi * (SS%par_Rup**2) * exp(-(SS%par_Velosity_inf**2)) * &
		(1.0 + exp(SS%par_Velosity_inf**2) * par_sqrtpi * SS%par_Velosity_inf * (1.0 + erf(SS%par_Velosity_inf))) / (2.0 * par_sqrtpi)
		SS%sqv_4 = (par_sqrtpi / 2.0) * ((par_Rmax**2) - (SS%par_Rup**2)) * &
			exp(-(SS%par_Velosity_inf**2)) * (par_sqrtpi * SS%par_Velosity_inf * erfc(SS%par_Velosity_inf) * exp(SS%par_Velosity_inf**2) - 1.0) ! �������� ������ (���� ���������)
		SS%sqv = SS%sqv_1 + SS%sqv_2 + SS%sqv_3 + SS%sqv_4
		
		
		
		SS%MK_N = MK_N1 + MK_N2 + MK_N3 + MK_N4  
		SS%MK_mu1 = (SS%sqv_1/ SS%sqv) * (1.0 * SS%MK_N / MK_N1)
		SS%MK_mu2 = (SS%sqv_2/ SS%sqv) * (1.0 * SS%MK_N / MK_N2)
		SS%MK_mu3 = (SS%sqv_3/ SS%sqv) * (1.0 * SS%MK_N / MK_N3)
		SS%MK_mu4 = (SS%sqv_4/ SS%sqv) * (1.0 * SS%MK_N / MK_N4)
		! Body of M_K_init
		SS%MK_N = SS%MK_N * par_n_potok * par_n_parallel
	end subroutine M_K_init


    subroutine Get_sensor_sdvig(SS, sdvig)
		! ��������� ������� ��������� ����� �� �����
		! Variables
		TYPE (Setka), intent(in out) :: SS
		integer, intent(in) :: sdvig
		logical :: exists
		integer(4) :: i, a, b, c, j
		
		
		inquire(file="rnd_my.txt", exist=exists)
		
		if (exists == .False.) then
			pause "net faila!!!  345434wertew21313edftr3e"
			STOP "net faila!!!"
		end if
		
		if (par_n_claster * par_n_potok * 2 > 1021) then
			print*, "NE XVATAET DATCHIKOV 31 miuhi8789pok9"
			pause
		end if
		
		open(1, file = "rnd_my.txt", status = 'old')
		
		do i = 1, sdvig
			read(1,*) a, b, c
			read(1,*) a, b, c
		end do
		
		do i = 1, par_n_potok
			read(1,*) a, b, c
			SS%sensor(:, 1, i) = (/ a, b, c /)
			read(1,*) a, b, c
			SS%sensor(:, 2, i) = (/ a, b, c /)
		end do
		
		close(1)
	end subroutine Get_sensor_sdvig

	subroutine M_K_Change_Velosity4(SS, potok, Ur, Uthe, Uphi, Vr, Vthe, Vphi, Wr_, Wthe_, Wphi_, mu_, cp, r, I_, x_ex, y_ex, z_ex, bb)
		! ��� ������ �����, �� ��������� ��� ��-������
		! ����������� ���� �������� � I_ �������������� ������
        TYPE (Setka), intent(in out) :: SS
		integer(4), intent(in) :: potok, I_
		real(8), intent(in) ::  Ur, Uthe, Uphi, Vr, Vthe, Vphi, cp, r, x_ex, y_ex, z_ex
		real(8), intent(out) ::   Wr_(I_ + 1), Wthe_(I_ + 1), Wphi_(I_ + 1), mu_(I_ + 1)
		logical, intent(out) :: bb
		
		real(8) :: X, uu
		real(8) :: gamma_(I_), Wa_(I_), Mho_(I_)
		real(8) :: ksi, gam1, gam2, Wr1, Wr2, Wr0, ksi1, ksi2, W1, W2, Wa, pp1, pp2, c, p, u
		real(8) :: p4, om1, om2, om3, lo, y1, y2, y3, v1, v2, v3, u1, u2, u3, uuu, yy, h, ksi3, ksi4, ksi5, ksi6, D, ko, gg
		integer(4) :: ii, met, k
		
		X = sqrt( (Vr - Ur)**2 + (Vthe - Uthe)**2 + (Vphi - Uphi)**2 )
		uu = exp(-(X**2)) / par_sqrtpi + (X + 1.0 / (2.0 * X)) * erf(X)
		Wr0 = -1.0
		met = 0

		do ii = 1, I_
			gamma_(ii) = 1.0 / ( (r / SS%MK_R_zone(ii))**2 - 1.0)
		end do

		

		! ��������� Wr
		
		do ii = 1, I_
			
			if (ii == 1) then
				gam1 = 0.0
				gam2 = gamma_(ii)
			else
				gam1 = gamma_(ii - 1)
				gam2 = gamma_(ii)
			end if
			
			
			
			pp1 = MK_for_Wr_1(0.0_8, gam2, Ur) - MK_for_Wr_1(0.0_8, gam1, Ur)
			pp2 = MK_for_Wr_2(0.0_8, gam2, Ur, Uthe) - MK_for_Wr_2(0.0_8, gam1, Ur, Uthe)
			
			
			call M_K_rand(SS%sensor(1, 2, potok), SS%sensor(2, 2, potok), SS%sensor(3, 2, potok), ksi)

			if (ksi < pp1 / (pp1 + pp2)) then
				call M_K_rand(SS%sensor(1, 2, potok), SS%sensor(2, 2, potok), SS%sensor(3, 2, potok), ksi)
				Wr1 = -3.0;
				Wr2 = 0.0;

				do while (MK_H_Wr_1(gam1, gam2, Wr1, Ur, pp1, ksi) >= 0.0)
					Wr1 = Wr1 - 1.0
				end do
					
				k = 0
				do while (dabs(Wr2 - Wr1) > 0.0001)     ! ������� �������, ����� �������������
					Wr0 = (Wr1 + Wr2) / 2.0
					if (MK_H_Wr_1(gam1, gam2, Wr0, Ur, pp1, ksi) < 0) then
						Wr1 = Wr0
					else
						Wr2 = Wr0
					end if
					k = k + 1
				end do
				Wr_(ii) = Wr1
			else
				call M_K_rand(SS%sensor(1, 2, potok), SS%sensor(2, 2, potok), SS%sensor(3, 2, potok), ksi)
				Wr1 = -3.0
				Wr2 = 0.0

				do while (MK_H_Wr_2(gam1, gam2, Wr1, Ur, Uthe, pp2, ksi) >= 0.0)
					Wr1 = Wr1 - 1.0
				end do
				
				k = 0
				do while (dabs(Wr2 - Wr1) > 0.0001)     ! ������� �������, ����� �������������
					Wr0 = (Wr1 + Wr2) / 2.0
					if (MK_H_Wr_2(gam1, gam2, Wr0, Ur, Uthe, pp2, ksi) < 0) then
						Wr1 = Wr0
					else
						Wr2 = Wr0
					end if
					k = k + 1
				end do
				Wr_(ii) = Wr1
			end if
		end do

		
		! �����������  Wa
		do ii = 1, I_
			
			if (ii == 1) then
				gam1 = 0.0
				gam2 = gamma_(ii)
			else
				gam1 = gamma_(ii - 1)
				gam2 = gamma_(ii)
			end if

			do while(.True.)
				call M_K_rand(SS%sensor(1, 2, potok), SS%sensor(2, 2, potok), SS%sensor(3, 2, potok), ksi1)
				call M_K_rand(SS%sensor(1, 2, potok), SS%sensor(2, 2, potok), SS%sensor(3, 2, potok), ksi2)
				W1 = sqrt(gam1 * (Wr_(ii)**2))
				W2 = sqrt(gam2 * (Wr_(ii)**2))
				Wa = sqrt(-log(exp(-(W1**2)) - ksi1 * (exp(-(W1**2)) - exp(-(W2**2)))))
				if ((1.0 + (Uthe * Wa)**2) / (1.0 + (Uthe * W2)**2) >= ksi2) EXIT
			end do
			Wa_(ii) = Wa
		end do

		
		! ������� ���� � Mho
		do ii = 1, I_
			c = Uthe * Wa_(ii)
			p = MK_norm_mho(c)
			
			Mho_(ii) = MK_play_mho(SS, potok, c)
			

			Wthe_(ii) = Wa_(ii) * cos(Mho_(ii))
			Wphi_(ii) = Wa_(ii) * sin(Mho_(ii))
			
			if (ii == 1) then
				gam1 = 0.0
				gam2 = gamma_(ii)
			else
				gam1 = gamma_(ii - 1)
				gam2 = gamma_(ii)
			end if

			u = sqrt( (Vr - Wr_(ii))**2 + (Vthe - Wthe_(ii))**2 +(Vphi - Wphi_(ii))**2 )

			if (X > 7.0) then
				mu_(ii) = (u * MK_sigma2(u, cp) / (uu * MK_sigma2(uu, cp))) * (MK_f2(0.0_8, gam1, Ur, Uthe) - MK_f2(0.0_8, gam2, Ur, Uthe)) * exp(-(Uthe**2)) * (p) / (1.0 + (c**2))
			else
				mu_(ii) = (u * MK_sigma2(u, cp) / (MK_int_1(X * cp, cp) / cp)) * (MK_f2(0.0_8, gam1, Ur, Uthe) - MK_f2(0.0_8, gam2, Ur, Uthe)) * exp(-(Uthe**2)) * (p) / (1.0 + (c**2))
			end if
		end do


		! �������� ��������� �����
		p4 = 0.5 * par_sqrtpi * X / (1.0 + 0.5 * par_sqrtpi * X)
		

		if (I_ > 0) then
			gg = gamma_(I_)
		else
			gg = 0.0
		end if

		do while(.True.)
			call M_K_rand(SS%sensor(1, 2, potok), SS%sensor(2, 2, potok), SS%sensor(3, 2, potok), ksi1)
			call M_K_rand(SS%sensor(1, 2, potok), SS%sensor(2, 2, potok), SS%sensor(3, 2, potok), ksi2)
			call M_K_rand(SS%sensor(1, 2, potok), SS%sensor(2, 2, potok), SS%sensor(3, 2, potok), ksi3)
			call M_K_rand(SS%sensor(1, 2, potok), SS%sensor(2, 2, potok), SS%sensor(3, 2, potok), ksi4)
			call M_K_rand(SS%sensor(1, 2, potok), SS%sensor(2, 2, potok), SS%sensor(3, 2, potok), ksi5)
			call M_K_rand(SS%sensor(1, 2, potok), SS%sensor(2, 2, potok), SS%sensor(3, 2, potok), ksi6)
			
			if (p4 < ksi1) then
				om1 = 1.0 - 2.0 * ksi4
				om2 = sqrt(1.0 - (om1**2)) * cos(2.0 * par_pi * ksi5)
				om3 = sqrt(1.0 - (om1**2)) * sin(2.0 * par_pi * ksi5)

				lo = sqrt(-log(ksi2 * ksi3))
				y1 = lo * om1
				y2 = lo * om2
				y3 = lo * om3
			else
				y1 = sqrt(-log(ksi2)) * cos(par_pi * ksi3)
				y2 = sqrt(-log(ksi4)) * cos(2.0 * par_pi * ksi5)
				y3 = sqrt(-log(ksi4)) * sin(2.0 * par_pi * ksi5)
			end if
			
			v1 = y1 + Ur
			v2 = y2 + Uthe
			v3 = y3 + Uphi
			u1 = Vr - v1
			u2 = Vthe - v2
			u3 = Vphi - v3
			uuu = sqrt(kvv(u1, u2, u3))
			yy = sqrt(kvv(y1, y2, y3))
			h = ((uuu * MK_sigma2(uuu, cp)) / (MK_sigma2(X, cp) * (X + yy)))
			
			if (h >= ksi6) EXIT
		end do


		Wr_(I_ + 1) = v1
		Wthe_(I_ + 1) = v2
		Wphi_(I_ + 1) = v3


		if (Wr_(I_ + 1) >= 0.0 .or. (Wthe_(I_ + 1)**2) + (Wphi_(I_ + 1)**2) > gg * (Wr_(I_ + 1)**2)) then
			mu_(I_ + 1) = 1.0
			bb = .True.
			return
		else
			mu_(I_ + 1) = 0.0  ! ����� �� ��������� ���� ����
			bb = .False.
			return
		end if

		bb = .True.
		return
	end subroutine M_K_Change_Velosity4

	subroutine MK_Velosity_initial2(SS, potok, Vx, Vy, Vz)
        TYPE (Setka), intent(in out) :: SS
		integer(4), intent(in) :: potok
		real(8), intent(out) :: Vx, Vy, Vz
		
		real(8) :: ksi1, ksi2, a, ksi3, ksi4, ksi5, ksi6, z, p1
		
		call M_K_rand(SS%sensor(1, 1, potok), SS%sensor(2, 1, potok), SS%sensor(3, 1, potok), ksi1)
		call M_K_rand(SS%sensor(1, 1, potok), SS%sensor(2, 1, potok), SS%sensor(3, 1, potok), ksi2)
		
		a = sqrt(-log(1.0 - ksi2))
		Vy = a * cos(2.0 * par_pi * ksi1)
		Vz = a * sin(2.0 * par_pi * ksi1)
		
		z = 0
		p1 = 0.5 * dabs(SS%par_Velosity_inf) * par_sqrtpi / (0.5 + 0.5 * dabs(SS%par_Velosity_inf) * par_sqrtpi);

		do
			call M_K_rand(SS%sensor(1, 1, potok), SS%sensor(2, 1, potok), SS%sensor(3, 1, potok), ksi3)
			call M_K_rand(SS%sensor(1, 1, potok), SS%sensor(2, 1, potok), SS%sensor(3, 1, potok), ksi4)
			call M_K_rand(SS%sensor(1, 1, potok), SS%sensor(2, 1, potok), SS%sensor(3, 1, potok), ksi5)
			call M_K_rand(SS%sensor(1, 1, potok), SS%sensor(2, 1, potok), SS%sensor(3, 1, potok), ksi6)

			if (p1 > ksi3) then
				z = cos(par_pi * ksi5) * sqrt(-log(ksi4))
			else
				z = sqrt(-log(1.0 - ksi4))
			end if
			
			if((dabs(z + SS%par_Velosity_inf) / (dabs(SS%par_Velosity_inf) + dabs(z)) > ksi6 .and. z >= -SS%par_Velosity_inf)) EXIT
		end do

		Vx = z + SS%par_Velosity_inf
		
		if (Vx <= 0.0) then
			print*, "Error iuygvbnmklo9890pljiouytrtyjhg"
		end if
		
		return
	end subroutine MK_Velosity_initial2
	
	subroutine MK_Velosity_initial(SS, potok, Vx, Vy, Vz)
        TYPE (Setka), intent(in out) :: SS
		integer(4), intent(in) :: potok
		real(8), intent(out) :: Vx, Vy, Vz
		
		real(8) :: ksi1, ksi2, a, ksi3, ksi4, ksi5, ksi6, z, p1
		
		call M_K_rand(SS%sensor(1, 1, potok), SS%sensor(2, 1, potok), SS%sensor(3, 1, potok), ksi1)
		call M_K_rand(SS%sensor(1, 1, potok), SS%sensor(2, 1, potok), SS%sensor(3, 1, potok), ksi2)
		
		a = sqrt(-log(ksi2))
		Vy = a * cos(2.0 * par_pi * ksi1)
		Vz = a * sin(2.0 * par_pi * ksi1)
		
		z = 0
		p1 = dabs(SS%par_Velosity_inf) * par_sqrtpi / (1.0 + dabs(SS%par_Velosity_inf) * par_sqrtpi)

		do
			call M_K_rand(SS%sensor(1, 1, potok), SS%sensor(2, 1, potok), SS%sensor(3, 1, potok), ksi3)
			call M_K_rand(SS%sensor(1, 1, potok), SS%sensor(2, 1, potok), SS%sensor(3, 1, potok), ksi4)
			call M_K_rand(SS%sensor(1, 1, potok), SS%sensor(2, 1, potok), SS%sensor(3, 1, potok), ksi5)
			call M_K_rand(SS%sensor(1, 1, potok), SS%sensor(2, 1, potok), SS%sensor(3, 1, potok), ksi6)

			if (p1 > ksi3) then
				z = cos(par_pi * ksi5) * sqrt(-log(ksi4))
			else
				if (ksi4 <= 0.5) then
					z = -sqrt(-log(2.0 * ksi4))
				else
					z = sqrt(-log(2.0 * (1.0 - ksi4)))
				end if
			end if
			
			if(dabs(z + SS%par_Velosity_inf) / (dabs(SS%par_Velosity_inf) + dabs(z)) >= ksi6 .and. z <= -SS%par_Velosity_inf) EXIT
		end do

		Vx = z + SS%par_Velosity_inf
		
		return
	end subroutine MK_Velosity_initial

    subroutine MK_ruletka(SS, n_potok, to_i, to_j, from_i, from_j, area, r, r_per, mu3, bb2)
        TYPE (Setka), intent(in out) :: SS
		real(8), intent(in out) :: mu3
		real(8), intent(in) :: r, r_per
		integer(4), intent(in) :: n_potok, to_i, to_j, from_i, from_j, area
		logical, intent(out) :: bb2
		
		real(8) :: mu, ksi
		
		mu = min( SS%MK_Mu(to_i, to_j, area) * 0.3 * SS%MK_SINKR(to_j) * (r_per/par_Rmax)**2, &
					SS%MK_Mu(from_i, from_j, area) * 0.3 * SS%MK_SINKR(from_j) * (r/par_Rmax)**2) 
		
		if (mu3 >= mu) then
			bb2 = .True.
			return
		else
			call M_K_rand(SS%sensor(1, 2, n_potok), SS%sensor(2, 2, n_potok), SS%sensor(3, 2, n_potok), ksi)
			if (mu3 >= ksi * mu) then
				mu3 = mu
				bb2 = .True.
				return
			else
				mu3 = 0.0
				bb2 = .False.
				return
			end if
		end if
		
		bb2 = .False.
		return
	
	end subroutine MK_ruletka

    subroutine MK_Init_Parametrs(SS, potok, mu_, Wt_, Wp_, Wr_, X_, bb)
		!! ��������� ������ �� ��������� ����� (�������� � ���������)
        TYPE (Setka), intent(in out) :: SS
		integer(4), intent(in) :: potok
		real(8), intent(out) :: mu_(par_n_zone + 1), Wt_(par_n_zone + 1), Wp_(par_n_zone + 1), Wr_(par_n_zone + 1), X_(par_n_zone + 1)
		logical, intent(out) :: bb
		
		real(8) :: Y, Wr1, Wr2, Wr0, W1, W2, Wa, Yt, c, p
		real(8) :: ksi, ksi1, ksi2
		real(8) :: X1
		real(8) :: X2
		real(8) :: X0                ! ��� ������ ����
		real(8) :: split, gam1, gam2
		real(8) :: Wa_(par_n_zone + 1)
		real(8) :: Mho_(par_n_zone + 1)
		integer(4) :: i
		real(8) :: p1, ksi3, ksi4, z, h, gg, p4_, ksi5, ksi6
		
		
		X0 = 1.0
		Y = dabs(SS%par_Velosity_inf)
		p1 = erf(Y) / (SS%MK_A1_ * (Y**2))
		
		

		! �����������  X
		do i = 1, par_n_zone
			
			call M_K_rand(SS%sensor(1, 1, potok), SS%sensor(2, 1, potok), SS%sensor(3, 1, potok), ksi)
			X1 = 0.0
			X2 = 1.0
			
			if (i == 1) then
				gam1 = 0.0
				gam2 = SS%MK_gam_zone(i)
			else
				gam1 = SS%MK_gam_zone(i - 1)
				gam2 = SS%MK_gam_zone(i)
			end if
			

			do while (dabs(X2 - X1) > 0.0000001)     ! ������� �������, ����� �������������
				X0 = (X1 + X2) / 2.0
				if (MK_Hx(gam1, gam2, X0, Y, ksi) < 0.0) then
					X1 = X0
				else
					X2 = X0
				end if
			end do
			X_(i) = X0
		end do

		! �����������  Wr
		do i = 1, par_n_zone
			
			call M_K_rand(SS%sensor(1, 1, potok), SS%sensor(2, 1, potok), SS%sensor(3, 1, potok), ksi)
			
			Wr1 = -5.0
			Wr2 = 0.0
			
			if (i == 1) then
				gam1 = 0.0
				gam2 = SS%MK_gam_zone(i)
			else
				gam1 = SS%MK_gam_zone(i - 1)
				gam2 = SS%MK_gam_zone(i)
			end if
			
			do while (MK_Hwr(gam1, gam2, Wr1, X_(i), Y, ksi) >= 0.0)
				Wr1 = Wr1 - 1.0
			end do
				
			do while (dabs(Wr2 - Wr1) > 0.000001)     ! ������� �������, ����� �������������
				Wr0 = (Wr1 + Wr2) / 2.0
				if (MK_Hwr(gam1, gam2, Wr0, X_(i), Y, ksi) < 0) then
					Wr1 = Wr0
				else
					Wr2 = Wr0
				end if
			end do

			Wr_(i) = Wr0
		end do

		
		! �����������  Wa
		
		do i = 1, par_n_zone
			Yt = Y * sqrt(1.0 - (X_(i)**2))
			
			if (i == 1) then
				gam1 = 0.0
				gam2 = SS%MK_gam_zone(i)
			else
				gam1 = SS%MK_gam_zone(i - 1)
				gam2 = SS%MK_gam_zone(i)
			end if

			do while(.True.)
				call M_K_rand(SS%sensor(1, 1, potok), SS%sensor(2, 1, potok), SS%sensor(3, 1, potok), ksi1)
				call M_K_rand(SS%sensor(1, 1, potok), SS%sensor(2, 1, potok), SS%sensor(3, 1, potok), ksi2)
				W1 = sqrt(gam1 * (Wr_(i)**2))
				W2 = sqrt(gam2 * (Wr_(i)**2))
				Wa = sqrt(-log(exp(-(W1**2)) - ksi1 * (exp(-(W1**2)) - exp(-(W2**2)))))
				if ((1.0 + (Yt * Wa)**2) / (1.0 + (Yt * W2)**2) >= ksi2) EXIT
			end do


			Wa_(i) = Wa
		end do

		! ������� ����
		do i = 1, par_n_zone
			
			c = Y * sqrt(1.0 - (X_(i)**2)) * Wa_(i)
			p = MK_norm_mho(c)
			
			
			Mho_(i) = MK_play_mho(SS, potok, c)
			

			Wt_(i) = Wa_(i) * cos(Mho_(i))
			Wp_(i) = Wa_(i) * sin(Mho_(i))
			
			if (i == 1) then
				gam1 = 0.0
				gam2 = SS%MK_gam_zone(i)
			else
				gam1 = SS%MK_gam_zone(i - 1)
				gam2 = SS%MK_gam_zone(i)
			end if
			
			mu_(i) = (MK_F(1.0_8, gam2, Y) - MK_F(1.0_8, gam1, Y)) * (p) / (SS%MK_A0_ * (1.0 + (c**2)))
		end do

		! ����������� �������� ����

		! �������� X ������� �� ���������

		
		call M_K_rand(SS%sensor(1, 1, potok), SS%sensor(2, 1, potok), SS%sensor(3, 1, potok), ksi1)
		
		
		if (p1 < ksi1) then
			
			do while(.True.)
				call M_K_rand(SS%sensor(1, 1, potok), SS%sensor(2, 1, potok), SS%sensor(3, 1, potok), ksi2)
				call M_K_rand(SS%sensor(1, 1, potok), SS%sensor(2, 1, potok), SS%sensor(3, 1, potok), ksi3)
				X2 = sqrt(ksi2)
				h = (1.0 + erf(X2 * Y)) / (1.0 + erf(Y))
				if (h > ksi3) EXIT
			end do
		else
			do
				call M_K_rand(SS%sensor(1, 1, potok), SS%sensor(2, 1, potok), SS%sensor(3, 1, potok), ksi2)
				call M_K_rand(SS%sensor(1, 1, potok), SS%sensor(2, 1, potok), SS%sensor(3, 1, potok), ksi3)
				X2 = (1.0 / Y) * sqrt(-log(ksi2)) * cos(par_pi * ksi3 / 2.0)
				if(X2 < 1.0) EXIT
			end do
		end if
		

		X_(par_n_zone + 1) = X2

		gg = 0.0
		
		if (par_n_zone > 0) gg = SS%MK_gam_zone(par_n_zone)
		
		p4_ = par_sqrtpi * (Y * X2) / (1.0 + par_sqrtpi * (Y * X2))                 ! ��� ��������� ��������� ����� �� ������� �� ��������� ������
		
		do while(.True.)
			call M_K_rand(SS%sensor(1, 1, potok), SS%sensor(2, 1, potok), SS%sensor(3, 1, potok), ksi1)
			call M_K_rand(SS%sensor(1, 1, potok), SS%sensor(2, 1, potok), SS%sensor(3, 1, potok), ksi2)
			call M_K_rand(SS%sensor(1, 1, potok), SS%sensor(2, 1, potok), SS%sensor(3, 1, potok), ksi3)
			if (p4_ > ksi1) then
				z = sqrt(-log(ksi2)) * cos(par_pi * ksi3)
			else
				if (ksi2 <= 0.5) then
					z = -sqrt(-log(2.0 * ksi2))
				else
					z = sqrt(-log(2.0 * (1.0 - ksi2)))
				end if
			end if

			Wr_(par_n_zone + 1) = z - (Y * X2)
			h = dabs(-(Y * X2) + z) / ((Y * X2) + dabs(z))
			call M_K_rand(SS%sensor(1, 1, potok), SS%sensor(2, 1, potok), SS%sensor(3, 1, potok), ksi4)
			if (h > ksi4 .and. z < (Y * X2)) EXIT
		end do

			
			call M_K_rand(SS%sensor(1, 1, potok), SS%sensor(2, 1, potok), SS%sensor(3, 1, potok), ksi5)
			call M_K_rand(SS%sensor(1, 1, potok), SS%sensor(2, 1, potok), SS%sensor(3, 1, potok), ksi6)

			Wt_(par_n_zone + 1) = Y * sqrt(1.0 - (X_(par_n_zone + 1)**2)) + sqrt(-log(ksi5)) * cos(2.0 * par_pi * ksi6)
			Wp_(par_n_zone + 1) = sqrt(-log(ksi5)) * sin(2.0 * par_pi * ksi6)


		if (par_n_zone == 0) then
			mu_(par_n_zone + 1) = 1.0
			bb = .True.
			return
		else
			if (Wr_(par_n_zone + 1) >= 0.0 .or. Wt_(par_n_zone + 1)**2 + Wp_(par_n_zone + 1)**2 > &
				SS%MK_gam_zone(par_n_zone) * Wr_(par_n_zone + 1)**2) then
				mu_(par_n_zone + 1) = 1.0
			else
				mu_(par_n_zone + 1) = 0.0  ! ����� �� ��������� ���� ����
				bb = .False.
				return
			end if
		end if

		bb = .True.
		return
	end subroutine MK_Init_Parametrs

    subroutine M_K_rand(s1, s2, s3, b)
		integer(4), intent(in out) :: s1, s2, s3
		real(8), intent(out) :: b
		integer(4):: ic15 = 32768, ic10 = 1024
		integer(4):: mz = 710, my = 17784, mx = 11973
		real(8):: xi = 9.0949470177292824E-13, c = 1.073741824E9
		integer(4) :: i13, i12, i11, ii
		i13 = mz * s1 + my * s2 + mx * s3
		i12 = my * s1 + mx * s2
		i11 = mx * s1
		ii = i11 / ic15
		i12 = i12 + ii
		s1 = i11 - ic15 * ii
		ii = i12 / ic15
		i13 = i13 + ii
		s2 = i12 - ic15 * ii
		s3 = mod(i13,ic10)
		b = xi * (c * s3 + ic15 * s2 + s1)
	end subroutine M_K_rand

    subroutine MK_Distination(SS, r, V, to_i, to_j, rr)
		! Variables
        TYPE (Setka), intent(in) :: SS
		real(8), intent(in) :: r(3), V(3)
		real(8), intent(out) :: rr
		integer(4), intent(out) :: to_i, to_j
		
		real(8) :: time_do_peregel, rk(3)
		
		time_do_peregel = -DOT_PRODUCT(r, V) / kvv(V(1), V(2), V(3))
		rk = r + time_do_peregel * V
		rr = norm2(rk)
		
		to_i = MK_geo_zones(SS, rr, 1.0_8)
		to_j = MK_alpha_zones(SS, polar_angle( rk(1), sqrt(rk(2)**2 + rk(3)**2) ) )
	
	end subroutine MK_Distination

    integer(4) pure function MK_alpha_zones(SS, al)
        TYPE (Setka), intent(in) :: SS
		real(8), intent(in) :: al
		integer(4) :: i
		
		do i = 1, par_m_zone
			if(al < SS%MK_al_zone(i)) then
				MK_alpha_zones = i
				return
			end if
		end do
		
		MK_alpha_zones = par_m_zone + 1
		return
		
	end function MK_alpha_zones
	
	integer(4) pure function MK_geo_zones(SS, r, k)
		! � ����� �������������� ���� ������ ��������� ���� (���� ��������� � ����)
        TYPE (Setka), intent(in) :: SS
		real(8), intent (in) :: r, k
		integer(4) :: i

		do i = 1, par_n_zone
			if (r < k * SS%MK_R_zone(i)) then
				MK_geo_zones = i
				return
			end if
		end do
		
		MK_geo_zones = par_n_zone + 1
		
	end function MK_geo_zones

    real(8) pure function MK_F(X, gam, Y)
		real(8), intent (in) :: X, gam, Y
		real(8) :: bb, gg, Y2, Y4, X2, A, B, C, D
		
		bb = sqrt(1.0 + gam)
		gg = 1.0 + gam
		Y2 = Y**2
		Y4 = Y**4
		X2 = X**2
		A = exp(-Y2) / (2.0 * par_sqrtpi * Y * bb**7)
		B = 2 * X * Y * bb**3 * (2.0 + gam * (3.0 + (X2 - 1.0) * Y2 + gam))
		C = par_sqrtpi * gg * (gg * (gg * (4.0 + gam) + Y2 * (2.0 + 3.0 * gam)) + &
			exp(X2 * Y2 / gg) * (2.0 * X2 * (X2 - 1.0) * Y4 * gam - (gg**2) * (4.0 + gam) + &
				Y2 * gg * (-2.0 - 3.0 * gam + X2 * (2.0 + gam))))
		D = par_sqrtpi * gg * exp(X2 * Y2 / gg) * &
			(-2.0 * X2 * (X2 - 1.0) * Y4 * gam + gg**2 * (4.0 + gam) - Y2 * gg * (-2.0 - 3.0 * gam + X2 * (2.0 + gam))) * &
			erf(X * Y / bb)
			
		MK_F =  A * (B + C - D)
	end function MK_F
	
	real(8) pure function MK_FI(Z, X, gam, Y)
		real(8), intent (in) :: Z, X, gam, Y
		real(8) :: bb, gg, X2, Y2, Y4, Z2, A, B, C, D, ee
		
		bb = sqrt(1.0 + gam)
		gg = 1.0 + gam
		X2 = X**2
		Y2 = Y**2
		Y4 = Y**4
		Z2 = Z**2
		A = exp(-Y2 - 2.0 * X * Y * Z - Z2 * gg) / (par_sqrtpi * bb**7)
		ee = exp((X * Y + Z + Z * gam)**2 / gg)
		B = ee * par_sqrtpi * X * Y * (2.0 * X2 * (X2 - 1.0) * Y4 * gam - 2.0 * gg**2 + (X2 - 1.0) * Y2 * gg * (2.0 + 5.0 * gam))
		C = 2.0 * bb * (X2 * (X2 - 1.0) * Y4 * gam - X * (X2 - 1.0) * Y**3 * Z * gam * gg - gg**2 +  &
			(X2 - 1.0) * Y2 * gg * (1 + gam * (2.0 + Z2 * gg)))
		D = ee * par_sqrtpi * X * Y * (2.0 * X2 * (X2 - 1.0) * Y4 * gam - 2.0 * gg**2 + (X2 - 1.0) * Y2 * gg * (2.0 + 5.0 * gam)) * &
			erf((X * Y + Z + Z * gam) / bb)
			
		MK_FI =  A * (B + C + D)
	end function MK_FI
	
	real(8) pure function MK_Hx(gam1, gam2, X, Y, ksi)
		real(8), intent (in) :: gam1, gam2, X, Y, ksi
		MK_Hx = MK_F(X, gam2, Y) - MK_F(X, gam1, Y) - ksi * (MK_F(1.0_8, gam2, Y) - MK_F(1.0_8, gam1, Y))
	end function MK_Hx
	
	real(8) pure function MK_Hwr(gam1, gam2, Z, X, Y, ksi)
		real(8), intent (in) :: gam1, gam2, Z, X, Y, ksi
		MK_Hwr = MK_FI(Z, X, gam2, Y) - MK_FI(Z, X, gam1, Y) - ksi * (MK_FI(0.0_8, X, gam2, Y) - MK_FI(0.0_8, X, gam1, Y))
	end function MK_Hwr
	
	real(8) pure function MK_norm_mho(c)
		real(8), intent (in) :: c
		real(8) :: s, d, cc, nn
		integer(4) :: i
		
		s = 1.0
		d = 0.0
		cc = 1.0
		nn = 1.0
		
		do i = 1, 99
			cc = cc * c
			nn = nn * i
			d = (cc / nn)**2
			s = s + d
			if (d < 0.00001) then
				EXIT
			end if
		end do
		
		MK_norm_mho = s
	end function MK_norm_mho
	
	real(8) function MK_play_mho(SS, potok, c)
        TYPE (Setka), intent(in out) :: SS
		real(8), intent(in) :: c
		integer(4), intent(in) :: potok
		real(8) :: x, ksi
		
		do while(.True.)
			call M_K_rand(SS%sensor(1, 2, potok), SS%sensor(2, 2, potok), SS%sensor(3, 2, potok), ksi)
			x = 2.0 * ksi * par_pi
			call M_K_rand(SS%sensor(1, 2, potok), SS%sensor(2, 2, potok), SS%sensor(3, 2, potok), ksi)
			if(exp(2.0 * c * cos(x)) > ksi * exp(2.0 * dabs(c)) ) EXIT
		end do

		MK_play_mho = x
	end function MK_play_mho

	real(8) pure function MK_Velosity_1(u, cp)
		real(8), intent (in) :: u, cp
		
		if (u < 0.00001) then
			MK_Velosity_1 = 2.0 * cp / par_sqrtpi + 2.0 * u * u / (3.0 * cp * par_sqrtpi) - &
				u * u * u * u / (15.0 * cp * cp * cp * par_sqrtpi)
		else
			MK_Velosity_1 =  exp(-u * u / cp**2) * cp / par_sqrtpi + (u + (cp**2) / (2.0 * u)) * erf(u / cp)
		end if
		
	end function MK_Velosity_1
	
	real(8) pure function MK_Velosity_2(u, cp)
		real(8), intent (in) :: u, cp
		
		if (u < 0.00001) then
			MK_Velosity_2 = (8.0 / 3.0) * cp**4 * par_pi * u + (8.0 / 15.0) * cp**2 * par_pi * u * u * u - &
				(4.0 / 105.0) * par_pi * u**5
		else
			MK_Velosity_2 =  cp**3 * par_pi * (exp(-u * u / cp**2) * cp * u * 2.0 * (cp**2 + 2.0 * u**2) + &
			par_sqrtpi * (4.0 * u**4 + 4.0 * cp * cp * u**2 - cp**4) * erf(u / cp)) / (4.0 * u * u)
		end if
		
	end function MK_Velosity_2
	
	real(8) pure function MK_Velosity_3(u, cp)
		real(8), intent (in) :: u, cp
		
		if (u < 0.00001) then
			MK_Velosity_3 = 8.0 * cp / (3.0 * par_sqrtpi) + 8.0 * u**2 / (9.0 * cp * par_sqrtpi) - &
				44.0 * u**4 / (135.0 * cp * cp * cp * par_sqrtpi)
		else
			MK_Velosity_3 =  exp(-u**2 / cp**2) * cp * (5.0 * cp**2 + 2.0 * u**2) / (par_sqrtpi * (3.0 * cp**2 + 2.0 * u**2)) + &
			(4.0 * u**4 + 12.0 * cp**2 * u**2 + 3.0 * cp**4) * erf(u / cp) / (2.0 * u * (3.0 * cp**2 + 2.0 * u**2))
		end if
		
	end function MK_Velosity_3
	
	real(8) pure function MK_int_1_f1(x)

		real(8), intent (in) :: x
	
		if (x <= 1.0) then
			MK_int_1_f1 =  6.283185155644284 + 0.000024846677279866114 * x + &
				2.0934329078277405 * x * x + 0.008055998193903208 * x * x * x - &
				0.2355169235647438 * x * x * x * x + 0.03820480582423355 * x * x * x * x * x + &
				0.006992274370591744 * x * x * x * x * x * x
		else if (x <= 3.0) then
			MK_int_1_f1 =  6.437524091973454 - 0.6331520099380095 * x + &
				3.1348881317268997 * x * x - 0.8454201478027856 * x * x * x + &
				0.1004702004260311 * x * x * x * x + 0.0009895488638964746 * x * x * x * x * x - &
				0.000920750276197054 * x * x * x * x * x * x
		else if (x <= 5) then
			MK_int_1_f1 =  4.4920780630505135 + 2.5133093267020654 * x + &
				1.1327223176567935 * x * x - 0.24648691152318875 * x * x * x + &
				0.031326738629523766 * x * x * x * x - 0.0021366031960331384 * x * x * x * x * x + &
				0.00005954097505746697 * x * x * x * x * x * x
		else if (x <= 7) then
			MK_int_1_f1 =  1.9138683588136232 + 5.350374732905213 * x - &
				0.16380205801427633 * x * x + 0.06765898334856263 * x * x * x - &
				0.011071118267864083 * x * x * x * x + 0.0008673476933852199 * x * x * x * x * x - &
				0.00002691859374483661 * x * x * x * x * x * x
		else if (x <= 50.0) then
			MK_int_1_f1 =  1.3138472469154294 + 5.336877156136497 * x + &
				0.020286308991329983 * x * x - 0.9780973533544137 * (x / 10.0)**3 + &
				0.26354051936651874 * (x / 10.0)**4 - 0.03711733070841712 * (x / 10.0)**5 + &
				0.002120935433043921 * (x / 10.0)**6
		else
			MK_int_1_f1 =  1.3138472469154294 + 5.336877156136497 * x + &
				0.020286308991329983 * x * x - 0.9780973533544137 * (x / 10.0)**3 + &
				0.26354051936651874 * (x / 10.0)**4 - 0.03711733070841712 * (x / 10.0)**5 + &
				0.002120935433043921 * (x / 10.0)**6
		end if
				 
		return
	end function MK_int_1_f1
	
	real(8) pure function MK_int_1_f2(x)
		real(8), intent (in) :: x
		
		if (x <= 1.0) then
			MK_int_1_f2 =  1.328216167939543 - 0.000004545681954848391 * x + &
				2.537368073155103 * x * x - 0.0020584991728545624 * x * x * x - &
				0.03742568018912792 * x * x * x * x - 0.010312136385277346 * x * x * x * x * x + &
				0.002767736179209713 * x * x * x * x * x * x
		else if (x <= 3.0) then
			MK_int_1_f2 =  1.2959616295629246 + 0.1533684067037866 * x + &
				2.2354849981206106 * x * x + 0.3113395567715921 * x * x * x - &
				0.21656309882941488 * x * x * x * x + 0.041957500887605075 * x * x * x * x * x - &
				0.0029978773724628604 * x * x * x * x * x * x
		else if (x <= 5.0) then
			MK_int_1_f2 =  1.903643456971281 - 1.4801836911099535 * x + 3.973958664572268 * x * x - &
				0.6482729779428982 * x * x * x + 0.07665007314658864 * x * x * x * x - &
				0.005369758193338703 * x * x * x * x * x + 0.00016605531871992049 * x * x * x * x * x * x
		else if (x <= 7.0) then
			MK_int_1_f2 =  -4.484415105552316 + 5.3747429756690766 * x + &
				0.8892806582308143 * x * x + 0.09767316152573671 * x * x * x - &
				0.025704749778475783 * x * x * x * x + 0.0021937998296249206 * x * x * x * x * x - &
				0.00006928845984076111 * x * x * x * x * x * x
		end if
				 
		return
	end function MK_int_1_f2
	
	real(8) pure function MK_int_1_f3(x)
		real(8), intent (in) :: x
 
		if (x <= 1.0) then
			MK_int_1_f3 = 1.2938345594193854 - 0.000031719847351174835 * x + &
				1.3183710041280094 * x * x - 0.014150512069488197 * x * x * x + &
				0.4226114681928129 * x * x * x * x - 0.06985750969880078 * x * x * x * x * x - &
				0.015347864048406958 * x * x * x * x * x * x
		else if (x <= 3.0) then
			MK_int_1_f3 = 0.9667460440956788 + 1.336271810704016 * x - &
				0.8687355257991665 * x * x + 1.7676868273627229 * x * x * x - &
				0.2731222764016417 * x * x * x * x + 0.004801770033831665 * x * x * x * x * x + &
				0.001780776080720323 * x * x * x * x * x * x
		else if (x <= 5.0) then
			MK_int_1_f3 = 4.760566734123174 - 5.048204299463048 * x + 3.332342585228025 * x * x + &
				0.47584339615235993 * x * x * x - 0.12072786272726124 * x * x * x * x + &
				0.011870955604980658 * x * x * x * x * x - 0.0004580199652402304 * x * x * x * x * x * x
		else if (x <= 7.0) then
			MK_int_1_f3 = 9.370493362469261 - 10.848615619383615 * x + 6.423326878282571 * x * x - &
				0.4148977656870439 * x * x * x + 0.025300923044176957 * x * x * x * x - &
				0.0010108688120876522 * x * x * x * x * x + 0.00001864423130429156 * x * x * x * x * x * x
		end if
			 
		return
	end function MK_int_1_f3
	
	real(8) pure function MK_int_1(x, cp)
		real(8), intent (in) :: x, cp
		real(8) :: b
		
		b = 1.0 - par_a_2 * log(cp)
		MK_int_1 = (cp / (par_sqrtpi**3)) * (b * b * MK_int_1_f1(x / cp) - &
			2.0 * par_a_2 * b * MK_int_1_f2(x / cp) + par_a_2**2 * MK_int_1_f3(x / cp))
	end function MK_int_1
	
	real(8) pure function MK_int_2_f1(x)
		real(8), intent (in) :: x
	
		if (x <= 1.0) then
				MK_int_2_f1 =   8.377571213788123 * x + 0.00047608508679086725 * x * x + &
					1.6710478320575737 * x * x * x + 0.016857530811432916 * x * x * x * x - &
					0.15132474125724812 * x * x * x * x * x + 0.030723378194358945 * x * x * x * x * x * x
				return
		else if (x <= 3.0) then
				MK_int_2_f1 =   -0.11788367995598747 + 8.937936705157014 * x - &
						1.119886471634001 * x * x + 2.8831031948885917 * x * x * x - &
						0.735250146386514 * x * x * x * x + 0.10356311378423572 * x * x * x * x * x - &
						0.006231417172309398 * x * x * x * x * x * x
				return
		else if (x <= 5.0) then
				MK_int_2_f1 =   2.9044497739429858 + 2.757712415967557 * x + 4.239161941189675 * x * x + &
						0.36198838294786784 * x * x * x - 0.05737777787138304 * x * x * x * x + &
						0.004956250079677106 * x * x * x * x * x - 0.0001809238236975877 * x * x * x * x * x * x
				return
		else if (x <= 7.0) then
				MK_int_2_f1 =  41.6323689028892 - 38.118317864344135 * x + 22.211189528076645 * x * x -&
					3.8547348524931246 * x * x * x + 0.5000517174807501 * x * x * x * x -&
					0.03446294709493891 * x * x * x * x * x + 0.0009860204070962582 * x * x * x * x * x * x
				return
		end if
	
		MK_int_2_f1 =  0.0
	end  function MK_int_2_f1
	
	real(8) pure function MK_int_2_f2(x)

		real(8), intent (in) :: x
	
		if (x <= 1.0) then
			MK_int_2_f2 =  3.8653461103376667 * x + 0.0001975300512691014 * x * x + &
				2.4468141895384012 * x * x * x + 0.005987984681429616 * x * x * x * x - &
				0.06453987836713967 * x * x * x * x * x + 0.0066920981111229004 * x * x * x * x * x * x
			return
		else if (x <= 3.0) then
			MK_int_2_f2 =  -0.10983889480661446 + 4.321087890898017 * x - &
				0.7707850845797619 * x * x + 3.1237901158486583 * x * x * x - &
				0.31485222316123385 * x * x * x * x + 0.010270249760261791 * x * x * x * x * x + &
				0.0008259803934338584 * x * x * x * x * x * x
			return
		else if (x <= 5.0) then
			MK_int_2_f2 =  1.8468011847509729 + 1.1986396743254275 * x + &
				1.1421489029509448 * x * x + 2.606316149569781 * x * x * x - &
				0.2788783468089509 * x * x * x * x + 0.019815317035281846 * x * x * x * x * x - &
				0.0006379970557448899 * x * x * x * x * x * x
			return
		else if (x <= 7.0) then
			MK_int_2_f2 =  9.480707804348185 - 8.022988228952784 * x + 5.823555900242488 * x * x + &
				1.3277220473440972 * x * x * x - 0.08074921612981413 * x * x * x * x + &
				0.0033058587723954185 * x * x * x * x * x - 0.00006041810279926061 * x * x * x * x * x * x
			return
		end if
			 
		MK_int_2_f2 = 0.0
	end function MK_int_2_f2
		
	real(8) pure function MK_int_2_f3(x)
        real(8), intent (in) :: x
		if (x <= 1.0) then
			MK_int_2_f3 =  2.6106039258326 * x - 0.0008357997793049243 * x * x + &
				2.0764571907368174 * x * x * x - 0.03182275644841273 * x * x * x * x + &
				0.26310521962808975 * x * x * x * x * x - 0.06034325471643871 * x * x * x * x * x * x
			return
		else if (x <= 3.0) then
			MK_int_2_f3 =   0.20784760901369737 + 1.5920325291316857 * x + &
				2.0985329535259014 * x * x - 0.26286255221171206 * x * x * x + &
				1.4610329096120886 * x * x * x * x - 0.25626862029131897 * x * x * x * x * x + &
				0.01684969647300594 * x * x * x * x * x * x
			return
		else if (x <= 5.0) then
			MK_int_2_f3 =   -6.284115352064703 + 15.665162343948523 * x &
				- 10.766105772158252 * x**2 + 6.0821074614870465 * x**3 - &
				0.3181501403196319 * x**4 + 0.01232319194701587 * x**5 &
				- 0.00017890661550876597 * x**6
			return
		else if (x <= 7.0) then
			MK_int_2_f3 =   -4.355962170454177 + 13.835332665069274 * x &
				- 10.12766071646978 * x**2 + 5.999392227482686 * x**3 - &
				0.32171318647989955 * x**4 + 0.014181987261856027 * x**5 &
				- 0.00030579035551447497 * x**6
			return
		end if
				
		MK_int_2_f3 = 0.0
		
	end function MK_int_2_f3
	
	real(8) pure function MK_int_2(x, cp)
		real(8), intent (in) :: x, cp
		real(8) :: b
		
		b = 1.0 - par_a_2 * log(cp)
		MK_int_2 = -(cp * cp / (par_sqrtpi * par_sqrtpi * par_sqrtpi)) * (b * b * MK_int_2_f1(x / cp) - &
			2.0 * par_a_2 * b * MK_int_2_f2(x / cp) + (par_a_2**2) * MK_int_2_f3(x / cp))
	end function MK_int_2
	
	real(8) pure function MK_int_3(x, cp)
		real(8), intent (in) :: x, cp
		real(8) :: b
		
		b = 1.0 - par_a_2 * log(cp)
		MK_int_3 = (cp**3 / (par_sqrtpi**3)) * (b * b * MK_int_3_f1(x / cp) - 2.0 * par_a_2 * b * MK_int_3_f2(x / cp) + (par_a_2**2) * MK_int_3_f3(x / cp))
	end function MK_int_3

	real(8) pure function MK_int_3_f1(x)
		real(8), intent (in) :: x
		
		if (x <= 1.0) then
			MK_int_3_f1 = 12.566370586001975 - 0.00001944816384202852 * x + 12.567558607381049 * x**2 - &
				0.010507444068349692 * x**3 + 1.2911398125420694 * x**4 - 0.05048482363937502 * x**5 - &
				0.029947937607835744 * x**6
		else if (x <= 3.0) then
			MK_int_3_f1 = 12.451555448799724 + 0.40252674353016715 * x + 12.081033298182223 * x**2 + 0.12939193776331415 * x**3 + &
				1.478876561367302 * x**4 - 0.2237491583356496 * x**5 + 0.014474521138620033 * x**6
		else if (x <= 5.0) then
			MK_int_3_f1 = 8.281962852844913 + 9.783032429527339 * x + 3.165848344887614 * x**2 + 4.711462666119614 * x**3 + &
				0.13739453130827392 * x**4 - 0.012096015315889195 * x**5 + 0.0004514225555943018 * x**6
		else if (x <= 7.0) then
			MK_int_3_f1 = 17.29966669035025 + 2.1457227895928668 * x + 5.572082327920818 * x**2 + 4.4083748449004645 * x**3 + &
				0.13640200155890422 * x**4 - 0.00854302147917508 * x**5 + 0.00022205921430504255 * x**6
		end if
		
	end function MK_int_3_f1

	real(8) pure function MK_int_3_f2(x)
		real(8), intent (in) :: x
		
		if (x <= 1.0) then
			MK_int_3_f2 = 5.798024979296493 - 0.00001772671478406096 * x + 9.98769025405073 * x**2 - 0.0073593516014156535 * x**3 + &
				2.27820901023418 * x**4 - 0.03135086958655956 * x**5 - 0.030403716978821237 * x**6
		else if (x <= 3.0) then
			MK_int_3_f2 = 5.864728705834779 - 0.34799875480550213 * x + 10.760249318358127 * x**2 - 0.9493738978205943 * x**3 + &
				2.948810798239708 * x**4 - 0.29690823082625284 * x**5 + 0.015284639719532296 * x**6
		else if (x <= 5.0) then
			MK_int_3_f2 = -3.21810405152587 + 17.08054504204813 * x - 3.43427364926184 * x**2 + 5.342946243936558 * x**3 + &
				1.3459970589832742 * x**4 - 0.07448103623442687 * x**5 + 0.0021616888853670958 * x**6
		else if (x <= 7.0) then
			MK_int_3_f2 = -59.42860988375287 + 79.24709914283142 * x - 32.304295129821256 * x**2 + 12.558114389999929 * x**3 + &
				0.32138208180756495 * x**4 + 0.003976327853989966 * x**5 - 0.0003703205838879336 * x**6
		end if
		
	end function MK_int_3_f2

	real(8) pure function MK_int_3_f3(x)
		real(8), intent(in) :: x
		
		if (x <= 1.0) then
			MK_int_3_f3 = 3.915885322797866 + 0.000044284982651632276 * x + 7.778946280540595 * x**2 + 0.01990833982643192 * x**3 + &
				2.710531013102172 * x**4 + 0.09480498076917274 * x**5 + 0.026454483470905545 * x**6
		else if (x <= 3.0) then
			MK_int_3_f3 = 4.371710139648201 - 1.787430730965795 * x + 10.530386754306065 * x**2 - 1.9945949141813966 * x**3 + &
				3.310710912254515 * x**4 + 0.1356460090430464 * x**5 - 0.019853464614841342 * x**6
		else if (x <= 5.0) then
			MK_int_3_f3 = 4.6940208896644435 - 6.423081982183987 * x + 18.137713177954524 * x**2 - 7.298291873534797 * x**3 + &
				5.202095367169497 * x**4 - 0.20620743501245353 * x**5 + 0.005094092519200813 * x**6
		else if (x <= 7.0) then
			MK_int_3_f3 = 194.2826492092263 - 195.77327124311523 * x + 95.83253215419927 * x**2 - 23.99281397751645 * x**3 + &
				7.170549820159753 * x**4 - 0.32568903981020303 * x**5 + 0.007955090180048264 * x**6
		end if
		
	end function MK_int_3_f3

	real(8) pure function MK_for_Wr_1(Z, gam, ur)
		!! ��� ��������� Wr � ����������� ��-������ (������ �����)
		real(8), intent (in) :: Z, gam, ur
		
		MK_for_Wr_1 = -exp(-gam * ur**2 / (gam + 1.0)) * (1.0 + erf((-ur + gam * Z + Z) / sqrt(gam + 1.0))) / (2.0 * sqrt(gam + 1.0))
	end function MK_for_Wr_1

	real(8) pure function MK_for_Wr_2(Z, gam, ur, ut)
		!! ��� ��������� Wr � ����������� ��-������ 
		real(8), intent (in) :: Z, gam, ur, ut
		real(8) :: gam1, ur2
		
		gam1 = gam + 1.0  ! gam + 1
		ur2 = ur**2
		
		MK_for_Wr_2 =  1.0 / (4.0 * gam1 * gam1 * gam1 * par_sqrtpi) * ut**2 * exp((-1.0 - gam / gam1) * ur2 - gam1 * Z * Z) * &
			(2.0 * exp(ur * (gam * ur / gam1 + 2.0 * Z)) * gam * gam1 * (ur + Z + gam * Z) + &
				exp(ur2 + gam1 * Z * Z) * sqrt(gam1) * par_sqrtpi * (2.0 + gam * (5.0 + 3.0 * gam + 2.0 * ur2)) * &
				(-1.0 - erf((-ur + Z + gam * Z) / sqrt(gam1))))
	end function MK_for_Wr_2

	real(8) pure function MK_H_Wr_1(gam1, gam2,	V, ur, p, ksi)
		real(8), intent (in) :: gam1, gam2, V, ur, p, ksi
		MK_H_Wr_1 = MK_for_Wr_1(V, gam2, ur) - MK_for_Wr_1(V, gam1, ur) - ksi * p
	end function MK_H_Wr_1

	real(8) pure function MK_H_Wr_2(gam1, gam2, V, ur, ut, p, ksi)
		! ��� ��������� Wr � ����������� ��-������ 
		real(8), intent (in) :: gam1, gam2, V, ur, ut, p, ksi
		MK_H_Wr_2 =  MK_for_Wr_2(V, gam2, ur, ut) - MK_for_Wr_2(V, gam1, ur, ut) - ksi * p
	end function MK_H_Wr_2

	real(8) pure function MK_f2(V, gam, ur, ut)
		real(8), intent (in) :: V, gam, ur, ut
		real(8) :: b, b4, b2, A
		
		b = sqrt(1.0 + gam)
		b4 = b**4
		b2 = b**2

		A = (exp(2.0 * V * ur - (ur**2) - (V**2) * b2) / (2.0 * par_sqrtpi * b4))
		
		if (A <= 0.0000000000001) then
			MK_f2 = 0.0
			return
		end if
			
		MK_f2 = A * (-gam * (ut**2) * (ur + gam * V + V) + &
			(par_sqrtpi / (2.0 * b)) * exp( (V + gam * V - ur)**2 / b2) * (2.0 * b4 + (ut**2) * (b2 * (3.0 * gam + 2.0) + 2.0 * gam * (ur**2))) * &
			(1.0 + erf((-ur + gam * V + V) / b)))
	end function MK_f2


end module Monte_Karlo