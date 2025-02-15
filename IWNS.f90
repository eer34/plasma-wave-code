

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
        integer::fid,n,k,region_i
        integer::ierr,my_id,num_procs
        real(wp)::start_cpu_time,finish_cpu_time
        real(wp)::eigen
        complex(wp)::polar(3)
        call mpi_init(ierr)
        call mpi_comm_rank(mpi_comm_world,my_id,ierr)
        call mpi_comm_size(mpi_comm_world,num_procs,ierr)
        call cpu_time(start_cpu_time)

        fid=10
        n_error=1000
		kc_square=128.0_wp
		epsilon_i=1d-7
		epsilon_accuracy_limit=1d-5
		n_circle=400
		n_line=400
		epsilon_0=0.1_wp
        
        if (my_id==0) then
	  		open(fid,file='output.csv')
        end if
		ti_div_te=1.0_wp
        c_div_v_para_input=470000.0_wp
		do k=1,40
			k_para_c_div_omega_ci=340.0_wp*k
			omega_pe_div_omega_ce_input=1.0_wp
			k_para_rho_i_para_input=k_para_c_div_omega_ci/(c_div_v_para_input)**(0.5)
			k_para_rho_e_para_input=k_para_rho_i_para_input*(1.0_wp/1836.0_wp/ti_div_te)**(0.5)
			k_para_rho_e_per_input=k_para_rho_e_para_input
			k_per_rho_i_para_input=0.1_wp*k_para_rho_i_para_input
			k_per_rho_i_per_input=0.1_wp*k_para_rho_i_para_input
			k_per_rho_e_para_input=0.1_wp*k_para_rho_e_para_input
			k_per_rho_e_per_input=0.1_wp*k_para_rho_e_para_input
			call set_parameter(c_div_v_para_input,omega_pe_div_omega_ce_input,k_para_rho_i_para_input,k_para_rho_e_para_input,k_para_rho_e_per_input,k_per_rho_i_para_input,k_per_rho_i_per_input,k_per_rho_e_para_input,k_per_rho_e_per_input)

			do region_i=1,15
				left_edge=1000.0_wp+100.0_wp*(region_i-1)
				right_edge=left_edge+100.0_wp
				down_edge=-6.0_wp
				up_edge=4.0_wp 

				allocate(ans_z_solve(n_error))
				allocate(ans_mul_solve(n_error))
				allocate(ans_z_error(n_error))
				allocate(ans_f_solve(n_error))

				call zero_pole_location(cold_plasma_dispersion_function,ierr,left_edge,right_edge,down_edge,up_edge,kc_square,epsilon_i,epsilon_accuracy_limit,n_circle,n_line,epsilon_0,z_solve_number,ans_z_solve,ans_mul_solve,ans_z_error,ans_f_solve)
				
				if (my_id==0) then
					do n=1,z_solve_number
						write(*,*),n,':'
						write(*,*),'ans_z_solve are',ans_z_solve(n)
						write(*,*),'ans_mul_solve are',ans_mul_solve(n)
						write(*,*),'ans_z_error are',ans_z_error(n)
						write(*,*),'ans_f_solve are',ans_f_solve(n)
						call polarization(ans_z_solve(n),eigen,polar)
						write(fid,'(*(G30.7,:,",",X))') k_para_c_div_omega_ci,n,real(ans_z_solve(n)),aimag(ans_z_solve(n)),ans_mul_solve(n),ans_z_error(n),abs(ans_f_solve(n)),eigen,polar(1),polar(2),polar(3)
						
					end do
				end if
				deallocate(ans_z_solve)
				deallocate(ans_mul_solve)
				deallocate(ans_z_error)
				deallocate(ans_f_solve)
			end do
		end do
        call cpu_time(finish_cpu_time)
        if (my_id==0) then
			close(fid)
			write(*,*),'running time is',finish_cpu_time-start_cpu_time
		end if
		call mpi_finalize(ierr)
       end program IWNS

