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
		real(wp)::k_para_rho_d_para_input,k_per_rho_d_para_input,k_per_rho_d_per_input
		integer::fid_1,fid_2,fid_3,kind
		real(wp)::omega_ce_div_x,omega_pe_div_x_square_1,omega_pe_div_x_square_2
        call mpi_init(ierr)
        call mpi_comm_rank(mpi_comm_world,my_id,ierr)
        call mpi_comm_size(mpi_comm_world,num_procs,ierr)
        call cpu_time(start_cpu_time)

		fid_1=1
		fid_2=2

        if (my_id==0) then
			open(fid_1, file='left cutoff.csv')
			open(fid_2, file='resonance.csv')
        end if
		do k=1,2000
			omega_ce_div_x=1836.0_wp+1.0_wp*k
			omega_pe_div_x_square_1=left_cut_off_two_ion_species(omega_ce_div_x)
			omega_pe_div_x_square_2=resonance_two_ion_species(omega_ce_div_x)
			if (my_id==0) then
				write(fid_1,'(*(G30.7,:,",",X))') omega_pe_div_x_square_1,omega_ce_div_x
				write(fid_2,'(*(G30.7,:,",",X))') omega_pe_div_x_square_2,omega_ce_div_x
			end if
		end do
        call cpu_time(finish_cpu_time)
		
        if (my_id==0) then
			write(*,*),'running time is',finish_cpu_time-start_cpu_time
			close(fid_1)
			close(fid_2)
		end if
		call mpi_finalize(ierr)
       end program IWNS

