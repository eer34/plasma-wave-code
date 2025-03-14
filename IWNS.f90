

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
       	real(wp)::beta_in,kap_n_in,kap_ti_in,kap_te_in,k_para_rho_i_in,k_para_rho_e_in,k_x_rho_i_in,k_y_rho_i_in
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
		epsilon_accuracy_limit=1d-6
		n_circle=400
		n_line=400
		epsilon_0=0.1_wp
        
        if (my_id==0) then
	  		open(fid,file='output.csv')
        end if
		do k=1,1
			beta_in=0.001
			kap_n_in=0.04
			kap_ti_in=0.1
			kap_te_in=0.0
			k_para_rho_i_in=1.256*1d-2
			k_para_rho_e_in=-k_para_rho_i_in/(1836.0)**(0.5)
			k_x_rho_i_in=0.4
			k_y_rho_i_in=0.4

			call set_parameter_itg(beta_in,kap_n_in,kap_ti_in,kap_te_in,k_para_rho_i_in,k_para_rho_e_in,k_x_rho_i_in,k_y_rho_i_in)
			left_edge=-10.01_wp
			right_edge=-0.01_wp
			down_edge=-4*k_para_rho_i_in*100
			up_edge=6.8_wp 

			allocate(ans_z_solve(n_error))
			allocate(ans_mul_solve(n_error))
			allocate(ans_z_error(n_error))
			allocate(ans_f_solve(n_error))

			call zero_pole_location(dispersion_function_itg,ierr,left_edge,right_edge,down_edge,up_edge,kc_square,epsilon_i,epsilon_accuracy_limit,n_circle,n_line,epsilon_0,z_solve_number,ans_z_solve,ans_mul_solve,ans_z_error,ans_f_solve)
			
			if (my_id==0) then
				do n=1,z_solve_number
					write(*,*),n,':'
					write(*,*),'ans_z_solve are',ans_z_solve(n)
					write(*,*),'ans_mul_solve are',ans_mul_solve(n)
					write(*,*),'ans_z_error are',ans_z_error(n)
					write(*,*),'ans_f_solve are',ans_f_solve(n)
					write(fid,'(*(G30.7,:,",",X))') kap_n_in,kap_ti_in,kap_te_in,k_para_rho_i_in,k_x_rho_i_in,k_y_rho_i_in,real(ans_z_solve(n)),aimag(ans_z_solve(n)),ans_mul_solve(n),ans_z_error(n),abs(ans_f_solve(n))
					
				end do
			end if
			deallocate(ans_z_solve)
			deallocate(ans_mul_solve)
			deallocate(ans_z_error)
			deallocate(ans_f_solve)
        end do
        call cpu_time(finish_cpu_time)
        if (my_id==0) then
			close(fid)
			write(*,*),'running time is',finish_cpu_time-start_cpu_time
		end if
		call mpi_finalize(ierr)
       end program IWNS

