    program IWNS
        use zpl5
        use tool
        use mpi
        use iso_fortran_env,only:wp=>real64
    	implicit none
        integer::z_solve_number
        real(wp)::left_edge,right_edge,down_edge,up_edge
        real(wp)::kc_square,epsilon_i,epsilon_accuracy_limit,epsilon_min,epsilon_max,epsilon_0
        complex(wp),allocatable::ans_z_solve(:)
        integer,allocatable::ans_mul_solve(:)
        real(wp),allocatable::ans_z_error(:)
        complex(wp),allocatable::ans_f_solve(:)
        integer::n_circle,n_line,n_error
		real(wp)::ti_div_te
        real(wp)::c_div_v_para_input,omega_pe_div_omega_ce_input,k_para_rho_i_para_input,k_para_rho_e_para_input,k_para_rho_e_per_input,k_per_rho_i_para_input,k_per_rho_i_per_input,k_per_rho_e_para_input,k_per_rho_e_per_input
		real(wp)::k_para_c_div_omega_ci
        integer::fid_process,fid,n,k,region_i
        integer::ierr,my_id,num_procs
        real(wp)::start_cpu_time,finish_cpu_time
        real(wp)::eigen
        complex(wp)::polar(3)
		real(wp),allocatable::x_wave_cma_x(:),x_wave_cma_y(:)
		integer::direction,kmax,ti_number
		real(wp)::k_para_c_div_omega_ci_input,k_per_c_div_omega_ci_input
        call mpi_init(ierr)
        call mpi_comm_rank(mpi_comm_world,my_id,ierr)
        call mpi_comm_size(mpi_comm_world,num_procs,ierr)
        call cpu_time(start_cpu_time)

        fid=10
        n_error=1000
		kc_square=128.0_wp
		epsilon_i=1d-7
		epsilon_accuracy_limit=1d-6
		n_circle=400
		n_line=400
		epsilon_0=0.1_wp
        if (my_id==0) then
			open(fid, file='x_wave_boundary.csv')
        end if
		k_para_c_div_omega_ci_input=0.0_wp
		k_per_c_div_omega_ci_input=0.0_wp
		omega_pe_div_omega_ce_input=1.0_wp
		do k=1,200
			k_per_c_div_omega_ci_input=10.0_wp*k+1000.0_wp
			
			call set_parameter_cold(omega_pe_div_omega_ce_input,k_para_c_div_omega_ci_input,k_per_c_div_omega_ci_input)
			do region_i=1,300
				left_edge=0.01_wp+10.0_wp*(region_i-1)
				right_edge=left_edge+10.0_wp
				up_edge=0.6_wp 
				down_edge=-0.4_wp
				
				allocate(ans_z_solve(n_error))
				allocate(ans_mul_solve(n_error))
				allocate(ans_z_error(n_error))
				allocate(ans_f_solve(n_error))

				call zero_pole_location(cold_plasma_dispersion_function_two_ion_species,ierr,left_edge,right_edge,down_edge,up_edge,kc_square,epsilon_i,epsilon_accuracy_limit,n_circle,n_line,epsilon_0,z_solve_number,ans_z_solve,ans_mul_solve,ans_z_error,ans_f_solve)
				do n=1,z_solve_number
					if (my_id==0) then
						if (ans_mul_solve(n)>0) then
							call polarization(ans_z_solve(n),eigen,polar)
						else 
							eigen=10000.0_wp
							polar(:)=0.0_wp
						end if
						write(fid,'(*(G30.7,:,",",X))') omega_pe_div_omega_ce_input,k_para_c_div_omega_ci_input, k_per_c_div_omega_ci_input,real(ans_z_solve(n)),imag(ans_z_solve(n)),ans_mul_solve(n),eigen,polar(1),polar(2),polar(3)
					end if
				end do
				deallocate(ans_z_solve)
				deallocate(ans_mul_solve)
				deallocate(ans_z_error)
				deallocate(ans_f_solve)
			end do
		end do
        call cpu_time(finish_cpu_time)
		
        if (my_id==0) then
			write(*,*),'running time is',finish_cpu_time-start_cpu_time
			close(fid)
		end if
		call mpi_finalize(ierr)
       end program IWNS

