

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
        real(wp)::c_div_v_para_input,omega_pe_div_omega_square,omega_pe_div_omega_ce_input,omega_ce_div_omega,omega_div_omega_ci_input
		real(wp)::k_para_rho_i_input,refractive_para
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
		do k=1,30
			omega_ce_div_omega=1.2_wp
			omega_pe_div_omega_square=0.1_wp*k
			c_div_v_para_input=470000.0_wp
			omega_pe_div_omega_ce_input=omega_pe_div_omega_square/omega_ce_div_omega/omega_ce_div_omega
			omega_div_omega_ci_input=1836.0_wp/omega_ce_div_omega
			refractive_para=2.0_wp
			k_para_rho_i_input=refractive_para*omega_div_omega_ci_input/(c_div_v_para_input)**(0.5)
			
			call set_parameter_variable_k_per(c_div_v_para_input,omega_pe_div_omega_ce_input,omega_div_omega_ci_input,k_para_rho_i_input)
			do region_i=1,1
				if (region_i==1) then
					left_edge=0.01_wp
					right_edge=4.01_wp
					down_edge=-0.99_wp
					up_edge=1.01_wp   
				else
					left_edge=4.01_wp+(region_i-2)*5.0_wp
					right_edge=left_edge+5.0_wp
					down_edge=-0.9_wp
					up_edge=1.1_wp
				end if
		  
				allocate(ans_z_solve(1:n_error))
				allocate(ans_mul_solve(1:n_error))
				allocate(ans_z_error(1:n_error))
				allocate(ans_f_solve(1:n_error))
		  
				call zero_pole_location(dispersion_function_variable_k_per,ierr,left_edge,right_edge,down_edge,up_edge,kc_square,epsilon_i,epsilon_accuracy_limit,n_circle,n_line,epsilon_0,z_solve_number,ans_z_solve,ans_mul_solve,ans_z_error,ans_f_solve)
				!call zero_pole_location(cold_plasma_dispersion_function_k_per,left_edge,right_edge,down_edge,up_edge,kc_square,epsilon_i,epsilon_accuracy_limit,n_circle,n_line,epsilon_0,z_solve_number,ans_z_solve,ans_mul_solve,ans_z_error,ans_f_solve)
			
				do n=1,z_solve_number
					write(*,*),n,':'
					write(*,*),'ans_z_solve are',ans_z_solve(n)
					write(*,*),'ans_mul_solve are',ans_mul_solve(n)
					write(*,*),'ans_z_error are',ans_z_error(n)
					write(*,*),'ans_f_solve are',ans_f_solve(n)
					call polarization(ans_z_solve(n),eigen,polar)
					!call cold_polarization(ans_z_solve(n),eigen,polar)
					write(fid,'(*(G30.7,:,",",X))') omega_pe_div_omega_square,omega_ce_div_omega,refractive_para,real(ans_z_solve(n)),aimag(ans_z_solve(n)),ans_mul_solve(n),ans_z_error(n),abs(ans_f_solve(n)),eigen,polar(1),polar(2),polar(3)
				
				end do
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

