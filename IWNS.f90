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
		real(wp),allocatable::right_wave_cma_x(:),right_wave_cma_y(:)
		integer::direction,kmax
		real(wp)::wave_max_imag, wave_max_real, wave_max_real_last, least_damped_ratio
        call mpi_init(ierr)
        call mpi_comm_rank(mpi_comm_world,my_id,ierr)
        call mpi_comm_size(mpi_comm_world,num_procs,ierr)
        call cpu_time(start_cpu_time)

        fid=10
		fid_process=1
        n_error=1000
		kc_square=128.0_wp
		epsilon_i=1d-7
		epsilon_accuracy_limit=1d-6
		n_circle=400
		n_line=400
		epsilon_0=0.1_wp
        kmax=90
        if (my_id==0) then
	  		open(fid_process,file='output_in_process.csv')
        end if
		ti_div_te=1.0_wp
		c_div_v_para_input=470000.0_wp/(0.1_wp)
		k_per_rho_i_para_input=0.0_wp
		k_per_rho_i_per_input=0.0_wp
		k_per_rho_e_para_input=0.0_wp
		k_per_rho_e_per_input=0.0_wp
		allocate(right_wave_cma_x(110))
        allocate(right_wave_cma_y(110))
		do k=1,kmax
			omega_pe_div_omega_ce_input=10**(-1.5_wp+0.05_wp*k)
			direction=2
			k_para_rho_i_para_input=0.0_wp

			do while (direction/=0)
				if (direction==2) then
					k_para_rho_i_para_input=k_para_rho_i_para_input+1.0_wp
				else if (direction==-1) then
					k_para_rho_i_para_input=k_para_rho_i_para_input-0.1_wp
                else if (direction==1) then
                    k_para_rho_i_para_input=k_para_rho_i_para_input+0.01_wp
				end if 
				if (k_para_rho_i_para_input<=0) then
					wave_max_real=1.0_wp
					exit
				end if
				k_para_rho_e_para_input=-k_para_rho_i_para_input*(1.0_wp/1836.0_wp/ti_div_te)**(0.5)
				k_para_rho_e_per_input=k_para_rho_e_para_input
				
				call set_parameter(c_div_v_para_input,omega_pe_div_omega_ce_input,k_para_rho_i_para_input,k_para_rho_e_para_input,k_para_rho_e_per_input,k_per_rho_i_para_input,k_per_rho_i_per_input,k_per_rho_e_para_input,k_per_rho_e_per_input)

				wave_max_real=1000.0_wp
                wave_max_imag=-1000.0_wp
				least_damped_ratio=100000.0_wp
				do region_i=1,18
					if (region_i==1) then
						left_edge=0.1_wp
						right_edge=100.0_wp
						down_edge=max(-4*k_para_rho_i_para_input,-4.0_wp)
						up_edge=4.0_wp 
					else
						left_edge=100.0_wp+100.0_wp*(region_i-2)
						right_edge=left_edge+100.0_wp
						down_edge=max(-20*k_para_rho_i_para_input,-0.1*left_edge)
						up_edge=4.0_wp 
					end if
					allocate(ans_z_solve(n_error))
					allocate(ans_mul_solve(n_error))
					allocate(ans_z_error(n_error))
					allocate(ans_f_solve(n_error))

					call zero_pole_location(dispersion_function_parallel_1,ierr,left_edge,right_edge,down_edge,up_edge,kc_square,epsilon_i,epsilon_accuracy_limit,n_circle,n_line,epsilon_0,z_solve_number,ans_z_solve,ans_mul_solve,ans_z_error,ans_f_solve)
					
					if (my_id==0) then
						do n=1,z_solve_number
							write(*,*),n,':'
							write(*,*),'ans_z_solve are',ans_z_solve(n)
							write(*,*),'ans_mul_solve are',ans_mul_solve(n)
							write(*,*),'ans_z_error are',ans_z_error(n)
							write(*,*),'ans_f_solve are',ans_f_solve(n)
						end do
					end if
					do n=1,z_solve_number
						if (abs(aimag(ans_z_solve(n)))<least_damped_ratio*real(ans_z_solve(n))) then
							wave_max_imag=aimag(ans_z_solve(n))
							wave_max_real=real(ans_z_solve(n))
							least_damped_ratio=abs(wave_max_imag)/wave_max_real
						end if 
					end do
					deallocate(ans_z_solve)
					deallocate(ans_mul_solve)
					deallocate(ans_z_error)
					deallocate(ans_f_solve)
				end do
				if (my_id==0) then
					write(fid_process,'(*(G30.7,:,",",X))') omega_pe_div_omega_ce_input, direction,k_para_rho_i_para_input,wave_max_real,wave_max_imag
				end if
				if (direction*abs(wave_max_imag)>0.01*direction*abs(wave_max_real)) then
					direction=-(direction+abs(direction)-2)/2
                
					if (direction==0) then
						wave_max_real=(wave_max_real+wave_max_real_last)/2
					end if 
					
				end if
				wave_max_real_last=wave_max_real
			end do
			right_wave_cma_x(k)=1836**2*omega_pe_div_omega_ce_input/wave_max_real**2  
			right_wave_cma_y(k)=1836.0_wp/wave_max_real
        end do
        call cpu_time(finish_cpu_time)
		close(fid_process)
		open(fid, file='omega_ce_cma.csv')
        if (my_id==0) then
			do k=1,90
				write(fid,'(*(G30.7,:,",",X))')  right_wave_cma_x(k) ,right_wave_cma_y(k)
			end do
			write(*,*),'running time is',finish_cpu_time-start_cpu_time
		end if
		close(fid)
		deallocate(right_wave_cma_x)
		deallocate(right_wave_cma_y)
		call mpi_finalize(ierr)
       end program IWNS

