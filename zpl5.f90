!****************************************************************************
!
!   PROGRAM: zplfortran
!
!   PURPOSE:  locate the zeros and poles of a meromorphic function.
!   
!   Author: 
!   Zilong Li     emails:lzl22@mails.tsinghua.edu.cn/lizilong@swip.ac.cn
!   Haotian Chen  emails:chenhaotian@swip.ac.cn
!****************************************************************************

module zpl5
use mpi
use quadpack_generic
use spline_interp
use lapack95
use iso_fortran_env,only:wp=>real64
implicit none

private pi,k_integral,angle_integral_signal,spl_angle_circle_1,spl_angle_circle_2,spl_angle_line_1,spl_angle_line_2,spl_lnf_circle_1,spl_lnf_circle_2,spl_lnf_line_1,spl_lnf_line_2
real(wp)::pi
type(spline_t)::spl_angle_circle_1,spl_angle_circle_2,spl_angle_line_1,spl_angle_line_2,spl_lnf_circle_1,spl_lnf_circle_2,spl_lnf_line_1,spl_lnf_line_2
integer::k_integral
real(wp)::angle_integral_signal
abstract interface

    real(wp) function func_1(x)
    !! interface for integral function.
        import :: wp
        implicit none
        real(wp), intent(in) :: x
    end function func_1
    
    complex(wp) function func_2(x)
    !! interface for user-supplied function.
        import :: wp
        implicit none
        complex(wp), intent(in) :: x
    end function func_2
end interface
    
    contains
    
    !zeros_poles_location:locate all the zeros and poles of a meromorphic function in the input region 
    subroutine zero_pole_location(fun,ierr,left_edge_0,right_edge_0,down_edge_0,up_edge_0,kc_square,epsilon_i,epsilon_accuracy_limit,n_circle,n_line,epsilon_0,z_solve_number,ans_z_solve,ans_mul_solve,ans_z_error,ans_f_solve)
        implicit none
        procedure(func_2)::fun  
        integer,intent(in)::ierr
        real(wp),intent(in)::left_edge_0,right_edge_0,down_edge_0,up_edge_0
        real(wp),intent(in)::kc_square,epsilon_i,epsilon_accuracy_limit,epsilon_0
        integer,intent(in)::n_circle,n_line
        integer,intent(inout)::z_solve_number
        complex(wp),intent(inout)::ans_z_solve(:)
        integer,intent(inout)::ans_mul_solve(:)
        real(wp),intent(inout)::ans_z_error(:)
        complex(wp),intent(inout)::ans_f_solve(:)
        integer::n_integral,repeat_condition,edge_good_condition
        integer::j,j_procs
        real(wp)::angle_max,d,d_remain,region_epsilon_min,condition,epsilon_accuracy
        real(wp)::left_edge,right_edge,down_edge,up_edge
        complex(wp),allocatable::f0_integral(:),f1_integral(:),f_integral_procs(:)
        complex(wp),allocatable::kexi(:),kexi_m(:),z_integral(:),zm_integral(:)
        
        integer::my_id,num_procs,n_integral_procs
        call mpi_comm_rank(mpi_comm_world,my_id,ierr)
        call mpi_comm_size(mpi_comm_world,num_procs,ierr)
        pi=acos(-1.0_wp)
        n_integral=2*n_circle+2*n_line
        n_integral_procs=n_integral/num_procs
        angle_max=pi-epsilon_0/2
        
        allocate(f_integral_procs(n_integral_procs))
        allocate(f0_integral(n_integral))
        allocate(f1_integral(n_integral))
        allocate(kexi(n_integral))
        allocate(kexi_m(n_integral))
        allocate(z_integral(n_integral))
        allocate(zm_integral(n_integral))
        z_solve_number=0
        up_edge=up_edge_0
        d_remain=up_edge-down_edge_0
        repeat_condition=0
        region_epsilon_min=1e-8_wp
        
      
        
        do while(d_remain>region_epsilon_min)
            left_edge=left_edge_0
            right_edge=right_edge_0
            d=min(d_remain,abs(right_edge-left_edge)/10)
            down_edge=up_edge-d
            edge_good_condition=0
            
            do while(edge_good_condition==0)
                call integral_point_compute(n_circle,n_line,left_edge,right_edge,down_edge,up_edge,angle_max,kexi,kexi_m,z_integral,zm_integral) 
                
                do j_procs=1,n_integral_procs
		    j=j_procs+n_integral_procs*my_id
                    if (j>=n_circle+n_line+1 .and. j<=2*n_circle+n_line+1) then   
                        if(up_edge>=up_edge_0 .or. repeat_condition/=1) then
                            f_integral_procs(j_procs)=fun(z_integral(j));   
                        else
			    f_integral_procs(j_procs)=f0_integral(2*n_circle+n_line+2-j); 
                        end if
                    else 
                        f_integral_procs(j_procs)=fun(z_integral(j));
                    end if
                    
                   
                end do
                
                call mpi_allgather(f_integral_procs,2*n_integral_procs,mpi_complex,f1_integral,2*n_integral_procs,mpi_complex,mpi_comm_world,ierr)
                
                do j=1,n_integral
		        if(isnan(abs(f1_integral(j)))) then
                        edge_good_condition=0
                        if(j<=n_circle+1) then
                            down_edge=down_edge-0.1_wp*d
                        else if (j<=n_circle+n_line+1) then
                            right_edge=right_edge+0.01_wp*abs(right_edge-left_edge)
                            repeat_condition=0
                        else if (j<=2*n_circle+n_line+1) then
                            up_edge=up_edge+0.1_wp*d
                        else if (j<=n_integral) then
                            left_edge=left_edge-0.01_wp*abs(right_edge-left_edge)
                            repeat_condition=0
                        end if
                        exit
                  else if (j==n_integral) then
                        edge_good_condition=1
                        if(right_edge==right_edge_0 .and. left_edge==left_edge_0) then
                            repeat_condition=1
                        end if 
                   end if
		    end do
		
            end do
            call argument_solve(ierr,left_edge,right_edge,down_edge,up_edge,epsilon_0,epsilon_i,epsilon_accuracy_limit,n_circle,n_line,kc_square,f1_integral,condition,epsilon_accuracy,z_solve_number,ans_z_solve,ans_mul_solve,ans_z_error)
            if (condition>kc_square .or. epsilon_accuracy>epsilon_accuracy_limit) then
                call subdivision(fun,ierr,left_edge,right_edge,down_edge,up_edge,epsilon_0,epsilon_i,epsilon_accuracy_limit,n_circle,n_line,kc_square,f1_integral,z_solve_number,ans_z_solve,ans_mul_solve,ans_z_error)
            end if
            
            do j=1,n_integral
                f0_integral(j)=f1_integral(j)
            end do
            if (my_id==0) then
	            write(*,'(A,F10.3,A,F10.3,A,F10.3,A,F10.3,A)') 'the region with left_edge:',left_edge,' ,right_edge:',right_edge,' ,down_edge:',down_edge,' ,up_edge:',up_edge,'has been solved'
            end if
            up_edge=down_edge
            d_remain=up_edge-down_edge_0
            
        end do
            
	do j=1,z_solve_number
           ans_f_solve(j)=fun(ans_z_solve(j))    
	end do
        deallocate(f_integral_procs)
        deallocate(f0_integral)
        deallocate(f1_integral)
        deallocate(kexi)
        deallocate(kexi_m)
        deallocate(z_integral)
        deallocate(zm_integral)  
    end subroutine zero_pole_location
   
    !subdivision: subdivide the region into four smaller subregions in z-space and locate the zeros and poles in these subregions
    recursive subroutine subdivision(fun,ierr,left_edge_0,right_edge_0,down_edge_0,up_edge_0,epsilon_0,epsilon_i,epsilon_accuracy_limit,n_circle,n_line,kc_square,f0_integral,z_solve_number,ans_z_solve,ans_mul_solve,ans_z_error)
        implicit none
        procedure(func_2)::fun
        integer,intent(in)::ierr
        real(wp),intent(in)::left_edge_0,right_edge_0,down_edge_0,up_edge_0
        real(wp),intent(in)::epsilon_0,epsilon_i,epsilon_accuracy_limit
        integer,intent(in)::n_circle,n_line
        real(wp),intent(in)::kc_square
        complex(wp),intent(in)::f0_integral(:)
        integer,intent(inout)::z_solve_number
        complex(wp),intent(inout)::ans_z_solve(:)
        integer,intent(inout)::ans_mul_solve(:)
        real(wp),intent(inout)::ans_z_error(:)
        integer::edge_good_condition,j,j_procs,n_circle_1,n_circle_2,n_line_1,n_line_2,n_integral
        integer::region_i
        real(wp)::angle_max,x_center_edge,y_center_edge,left_edge,right_edge,down_edge,up_edge
        real(wp)::condition,epsilon_accuracy
        complex(wp),allocatable::f_xy_integral(:),f_edge_integral(:),f1_integral(:),f2_integral(:),f3_integral(:),f4_integral(:)
        complex(wp),allocatable::kexi(:),kexi_m(:),z_integral(:),zm_integral(:)
        complex(wp),allocatable::f_xy_integral_procs(:),f_edge_integral_procs(:)
        integer::my_id,num_procs
        integer::n_integral_procs
        call mpi_comm_rank(mpi_comm_world,my_id,ierr)
        call mpi_comm_size(mpi_comm_world,num_procs,ierr)
        
        edge_good_condition=0
        angle_max=pi-epsilon_0/2
        x_center_edge=(right_edge_0+left_edge_0)/2
        y_center_edge=(down_edge_0+up_edge_0)/2
        n_circle_1=n_circle
        n_circle_2=n_circle
        n_line_1=n_line
        n_line_2=n_line
        n_integral=2*n_line+2*n_circle
        n_integral_procs=n_integral/num_procs
        
        allocate(f_edge_integral_procs(2*n_integral_procs))
        allocate(f_edge_integral(2*n_integral))
	    allocate(kexi(2*n_integral))
	    allocate(kexi_m(2*n_integral))
	    allocate(z_integral(2*n_integral))
	    allocate(zm_integral(2*n_integral))
        call integral_point_compute(2*n_circle,2*n_line,left_edge_0,right_edge_0,down_edge_0,up_edge_0,angle_max,kexi,kexi_m,z_integral,zm_integral) 
        
        do j_procs=1,2*n_integral_procs
	    j=j_procs+2*n_integral_procs*my_id
	    if (mod(j,2)==0) then
	      f_edge_integral_procs(j_procs)=fun(z_integral(j))
	    else
	      f_edge_integral_procs(j_procs)=f0_integral((j+1)/2)
	    end if
        end do
        call mpi_allgather(f_edge_integral_procs,4*n_integral_procs,mpi_complex,f_edge_integral,4*n_integral_procs,mpi_complex,mpi_comm_world,ierr)
         
        !compute the function value of the discrete point on the horizontal and vertical lines.If any function value is too small or too large,move the boundary in 10 steps
        allocate(f_xy_integral_procs(n_integral_procs))
        allocate(f_xy_integral(n_integral))
    
        do while(edge_good_condition==0)
            
            do j_procs=1,n_integral_procs
		        j=j_procs+n_integral_procs*my_id
		        if (j<=2*n_circle) then
		          f_xy_integral_procs(j_procs)=fun(cmplx(left_edge_0+(right_edge_0-left_edge_0)/(2*n_circle)*(j-1),y_center_edge,wp))
		        else
		          f_xy_integral_procs(j_procs)=fun(cmplx(x_center_edge,down_edge_0+(up_edge_0-down_edge_0)/(2*n_line)*(j-2*n_circle-1),wp))
		        end if
            end do
            
            call mpi_allgather(f_xy_integral_procs,2*n_integral_procs,mpi_complex,f_xy_integral,2*n_integral_procs,mpi_complex,mpi_comm_world,ierr)
                       
            do j=1,n_integral              
                if (isnan(abs(f_xy_integral(j)))) then
		    if(j<=2*n_circle) then
		      edge_good_condition=0
		      n_line_1=n_line_1+20
		      n_line_2=n_line_2-20
		      y_center_edge=y_center_edge-10*abs(up_edge_0-down_edge_0)/n_line
		    else
		      edge_good_condition=0
		      n_circle_1=n_circle_1+20
		      n_circle_2=n_circle_2-20
		      x_center_edge=x_center_edge+10*abs(right_edge_0-left_edge_0)/n_circle
		    end if
                    exit
                else if (j==n_integral) then
                    edge_good_condition=1
                end if    
            end do
        end do
        
        allocate(f1_integral(2*n_circle_1+2*n_line_1))
        allocate(f2_integral(2*n_circle_2+2*n_line_1))
        allocate(f3_integral(2*n_circle_1+2*n_line_2))
        allocate(f4_integral(2*n_circle_2+2*n_line_2))
        
        do region_i=1,4
            select case(region_i)
            !Calculate the function values of discrete points on the boundaries of each subregion,
            !taking into account that the repeated boundary function values in adjacent subregion are only calculated once
            case(1)
                left_edge=left_edge_0
                right_edge=x_center_edge
                down_edge=y_center_edge
                up_edge=up_edge_0
               
		        do j=1,2*n_circle_1+2*n_line_1
		          if (j<=n_circle_1) then
		            f1_integral(j)=f_xy_integral(j)
		          else if (j<=(n_circle_1+n_line_1)) then
		            f1_integral(j)=f_xy_integral(j+n_line_2+n_circle_2)
		          else 
		            f1_integral(j)=f_edge_integral(j+n_line_2+2*n_circle_2)
		           end if
		        end do
		
                call argument_solve(ierr,left_edge,right_edge,down_edge,up_edge,epsilon_0,epsilon_i,epsilon_accuracy_limit,n_circle_1,n_line_1,kc_square,f1_integral,condition,epsilon_accuracy,z_solve_number,ans_z_solve,ans_mul_solve,ans_z_error)
             
                if (condition>kc_square .or. epsilon_accuracy>epsilon_accuracy_limit) then
                    call subdivision(fun,ierr,left_edge,right_edge,down_edge,up_edge,epsilon_0,epsilon_i,epsilon_accuracy_limit,n_circle_1,n_line_1,kc_square,f1_integral,z_solve_number,ans_z_solve,ans_mul_solve,ans_z_error)
                end if
            case(2)   
                left_edge=x_center_edge
                right_edge=right_edge_0
                down_edge=y_center_edge
                up_edge=up_edge_0
                
		        do j=1,2*n_circle_2+2*n_line_1
		          if (j<=n_circle_2) then
		            f2_integral(j)=f_xy_integral(n_circle_1+j)
		          else if (j<=(2*n_circle_2+n_line_1+1)) then
		            f2_integral(j)=f_edge_integral(j+n_line_2+n_circle_1)
		          else 
		            f2_integral(j)=f_xy_integral(n_integral+2*n_circle_2+n_line_1+2-j)
		           end if
		        end do
		
                call argument_solve(ierr,left_edge,right_edge,down_edge,up_edge,epsilon_0,epsilon_i,epsilon_accuracy_limit,n_circle_2,n_line_1,kc_square,f2_integral,condition,epsilon_accuracy,z_solve_number,ans_z_solve,ans_mul_solve,ans_z_error)
                if (condition>kc_square .or. epsilon_accuracy>epsilon_accuracy_limit) then
                    call subdivision(fun,ierr,left_edge,right_edge,down_edge,up_edge,epsilon_0,epsilon_i,epsilon_accuracy_limit,n_circle_2,n_line_1,kc_square,f2_integral,z_solve_number,ans_z_solve,ans_mul_solve,ans_z_error)
                end if
            case(3)   
                left_edge=left_edge_0
                right_edge=x_center_edge
                down_edge=down_edge_0
                up_edge=y_center_edge
		        do j=1,2*n_circle_1+2*n_line_2
		          if (j<=n_circle_1+1) then
		            f3_integral(j)=f_edge_integral(j)
		          else if (j<=(n_circle_1+n_line_2+1)) then
		            f3_integral(j)=f_xy_integral(j+n_circle_2)
		          else if (j<=(2*n_circle_1+n_line_2+1)) then
		            f3_integral(j)=f_xy_integral(n_line_2+2*n_circle_1+2-j)
		          else  
		            f3_integral(j)=f_edge_integral(j+2*n_circle_2+2*n_line_1)
		          end if
		        end do
		
		
                call argument_solve(ierr,left_edge,right_edge,down_edge,up_edge,epsilon_0,epsilon_i,epsilon_accuracy_limit,n_circle_1,n_line_2,kc_square,f3_integral,condition,epsilon_accuracy,z_solve_number,ans_z_solve,ans_mul_solve,ans_z_error)

                if (condition>kc_square .or. epsilon_accuracy>epsilon_accuracy_limit) then
                    call subdivision(fun,ierr,left_edge,right_edge,down_edge,up_edge,epsilon_0,epsilon_i,epsilon_accuracy_limit,n_circle_1,n_line_2,kc_square,f3_integral,z_solve_number,ans_z_solve,ans_mul_solve,ans_z_error)
                end if
                case(4)   
                    left_edge=x_center_edge
                    right_edge=right_edge_0
                    down_edge=down_edge_0
                    up_edge=y_center_edge
                    do j=1,2*n_circle_2+2*n_line_2
		              if (j<=n_circle_2+n_line_2+1) then
		                f4_integral(j)=f_edge_integral(j+n_circle_1)
		              else if (j<=(2*n_circle_2+n_line_2+1)) then
		                f4_integral(j)=f_xy_integral(n_circle_1+2*n_circle_2+n_line_2+2-j)
		              else 
		                f4_integral(j)=f_xy_integral(n_circle_1+3*n_circle_2+2*n_line_2+2-j)
		              end if
		            end do
		
                    call argument_solve(ierr,left_edge,right_edge,down_edge,up_edge,epsilon_0,epsilon_i,epsilon_accuracy_limit,n_circle_2,n_line_2,kc_square,f4_integral,condition,epsilon_accuracy,z_solve_number,ans_z_solve,ans_mul_solve,ans_z_error)
                
                    if (condition>kc_square .or. epsilon_accuracy>epsilon_accuracy_limit) then
                        call subdivision(fun,ierr,left_edge,right_edge,down_edge,up_edge,epsilon_0,epsilon_i,epsilon_accuracy_limit,n_circle_2,n_line_2,kc_square,f4_integral,z_solve_number,ans_z_solve,ans_mul_solve,ans_z_error)
                    end if
                
            end select
           
        end do
	    deallocate(kexi)
	    deallocate(kexi_m)
	    deallocate(z_integral)
	    deallocate(zm_integral)
	    deallocate(f_edge_integral_procs) 
        deallocate(f_edge_integral)
        deallocate(f_xy_integral_procs)
        deallocate(f_xy_integral)
        deallocate(f1_integral)
        deallocate(f2_integral)
        deallocate(f3_integral)
        deallocate(f4_integral)
    end subroutine subdivision
    
    !argumnet_solve:Given the boundary of the region and the function values on the boundary, use the argument principle and solve the general eigenvaules problem to locate the zeros and poles in this region 
    subroutine argument_solve(ierr,left_edge,right_edge,down_edge,up_edge,epsilon_0,epsilon_i,epsilon_accuracy_limit,n_circle,n_line,kc_square,f_integral,condition,epsilon_accuracy,z_solve_number,ans_z_solve,ans_mul_solve,ans_z_error)
    
        implicit none
        integer,intent(in)::ierr
        real(wp),intent(in)::left_edge,right_edge,down_edge,up_edge
        real(wp),intent(in)::epsilon_0,epsilon_i,epsilon_accuracy_limit
        integer,intent(in)::n_circle,n_line
        real(wp),intent(in)::kc_square
        complex(wp),intent(in)::f_integral(:)
        real(wp),intent(inout)::condition,epsilon_accuracy
        integer,intent(inout)::z_solve_number
        complex(wp),intent(inout)::ans_z_solve(:)
        integer,intent(inout)::ans_mul_solve(:)
        real(wp),intent(inout)::ans_z_error(:)
        integer::N_error=50
        integer::N=0
        real(wp)::epsilon_svd=1e-7_wp
        real(wp)::ratio_svd=1000.0_wp
        integer::n_integral,i,j,k
        complex(wp),allocatable::G(:)
        complex(wp),allocatable::kexi(:),kexi_m(:),z_integral(:),zm_integral(:)
        integer,allocatable::rank(:)
        real(wp)::angle_max,r_min,L
        complex(wp)::z0
        real(wp),allocatable::angle_integral(:)     
        complex(wp)::integral_result,integral_result_ref
        real(wp),allocatable::angle_spline(:),r_spline(:)
        real(wp),allocatable::angle_integral_circle_spline(:),lnf_integral_circle_spline(:),angle_integral_line_spline(:),lnf_integral_line_spline(:)
        real(wp)::S_circle_1_real,S_circle_1_imag,S_circle_2_real,S_circle_2_imag,S_line_1_real,S_line_1_imag,S_line_2_real,S_line_2_imag
        complex(wp),allocatable::hankel_0(:,:),hankel_1(:,:),hankel_0_ref(:,:),hankel_1_ref(:,:)
        complex(wp),allocatable::ans_kexi(:),ans_kexi_effective(:),ans_z(:),ans_kexi_ref(:),ans_kexi_effective_ref(:),ans_z_ref(:)
        integer,allocatable::ans_mul(:)
        integer::N_effective,N_effective_ref,my_id
        real(wp)::f_abs_max,f_abs_min
        
        call mpi_comm_rank(mpi_comm_world,my_id,ierr)
        n_integral=2*n_circle+2*n_line
        angle_max=pi-epsilon_0/2
        r_min=exp(-2*angle_max*abs(up_edge-down_edge)/abs(right_edge-left_edge))
        L=abs(right_edge-left_edge)
        z0=cmplx((right_edge+left_edge)/2,up_edge,wp)
        allocate(kexi(n_integral))
        allocate(kexi_m(n_integral))
        allocate(z_integral(n_integral))
        allocate(zm_integral(n_integral))
        allocate(G(2*N_error+1))
        allocate(rank(N_error+1))
        allocate(angle_integral(n_integral+1))
    
        allocate(angle_spline(n_circle+1))
        allocate(angle_integral_circle_spline(n_circle+1))
        allocate(lnf_integral_circle_spline(n_circle+1))
        allocate(r_spline(n_line+1))
        allocate(angle_integral_line_spline(n_line+1))
        allocate(lnf_integral_line_spline(n_line+1))
       
        
        !keep track of Arg(f) to make Arg(f) continuous
        call integral_point_compute(n_circle,n_line,left_edge,right_edge,down_edge,up_edge,angle_max,kexi,kexi_m,z_integral,zm_integral)
        call arg_trace(angle_integral,f_integral,n_integral)
        
        do i=1,n_circle+1
            angle_spline(i)=-angle_max+2*angle_max/n_circle*(i-1)
        end do
        do i=1,n_line+1
            r_spline(i)=r_min**(1-(i-1)*1.0_wp/n_line)
        end do
        
        f_abs_max=abs(f_integral(1))
        f_abs_min=abs(f_integral(1))
        
        do i=1,n_integral
            f_abs_max=max(f_abs_max,abs(f_integral(i)))
            f_abs_min=min(f_abs_min,abs(f_integral(i)))
            
        end do
        
        do i=1,n_circle+1
            angle_integral_circle_spline(i)=angle_integral(n_circle+2-i)
            lnf_integral_circle_spline(i)=log(abs(f_integral(n_circle+2-i)))-0.5*log(f_abs_max*f_abs_min)
        end do
                
        call spline_set_coeffs(angle_spline, angle_integral_circle_spline, n_circle+1, spl_angle_circle_1)
        call spline_set_coeffs(angle_spline, lnf_integral_circle_spline, n_circle+1, spl_lnf_circle_1)
        
        do i=1,n_circle+1
            angle_integral_circle_spline(i)=angle_integral(n_circle+n_line+i)
            lnf_integral_circle_spline(i)=log(abs(f_integral(n_circle+n_line+i)))-0.5*log(f_abs_max*f_abs_min)
        end do
        call spline_set_coeffs(angle_spline, angle_integral_circle_spline, n_circle+1, spl_angle_circle_2)
        call spline_set_coeffs(angle_spline, lnf_integral_circle_spline, n_circle+1, spl_lnf_circle_2)
        
        do i=1,n_line+1
            angle_integral_line_spline(i)=angle_integral(n_circle+i)
            lnf_integral_line_spline(i)=log(abs(f_integral(n_circle+i)))-0.5*log(f_abs_max*f_abs_min)
        end do
        call spline_set_coeffs(r_spline, angle_integral_line_spline, n_line+1, spl_angle_line_1)
        call spline_set_coeffs(r_spline, lnf_integral_line_spline, n_line+1, spl_lnf_line_1)
        
        do i=1,n_line+1
            if(i==1) then
                lnf_integral_line_spline(i)=log(abs(f_integral(1)))-0.5*log(f_abs_max*f_abs_min)
            else
                lnf_integral_line_spline(i)=log(abs(f_integral(n_integral+2-i)))-0.5*log(f_abs_max*f_abs_min)
            end if
            angle_integral_line_spline(i)=angle_integral(n_integral+2-i)
        end do
        call spline_set_coeffs(r_spline, angle_integral_line_spline, n_line+1, spl_angle_line_2)
        call spline_set_coeffs(r_spline, lnf_integral_line_spline, n_line+1, spl_lnf_line_2)
        !calculate the contour integral, construct Hankel matrix and rank,determine the total number of the zeros and poles
        integral_result=cmplx((angle_integral(n_integral+1)-angle_integral(1))/2/pi,0,wp)
    
        if (abs(integral_result)<epsilon_svd) then
            integral_result=0
            rank(1)=0
        else
            rank(1)=1
        end if 
        if(my_id==0)then
            write(*,*),'left_edge',left_edge,'right_edge',right_edge,'up_edge',up_edge,'down_edge',down_edge
            write(*,"('rank(',i2,'):',i2)"),1,rank(1)
        end if
        G(1)=integral_result
    
        do k=1,N_error
            do j=0,1
                integral_result_ref=(0,0)
                !The integration results of the trapezoidal integration method serve as a reference
                do i=1,n_integral-1
                    integral_result_ref=integral_result_ref-(2*k-1+j)*kexi(i+1)**(2*k-2+j)*(kexi_m(i+1)-kexi_m(i))*cmplx(log(abs(f_integral(i+1)))-0.5*log(f_abs_max*f_abs_min),angle_integral(i+1),wp)/cmplx(0,2*pi,wp)
                
                end do
            
                integral_result_ref=integral_result_ref-(2*k-1+j)*kexi(1)**(2*k-2+j)*(kexi_m(1)-kexi(1))*cmplx(log(abs(f_integral(1)))-0.5*log(f_abs_max*f_abs_min),angle_integral(1),wp)/cmplx(0,2*pi,wp)
                integral_result_ref=integral_result_ref-(2*k-1+j)*kexi(1)**(2*k-2+j)*(kexi(1)-kexi_m(n_integral))*cmplx(log(abs(f_integral(1)))-0.5*log(f_abs_max*f_abs_min),angle_integral(n_integral),wp)/cmplx(0,2*pi,wp)
                integral_result_ref=integral_result_ref+kexi(1)**(2*k-1+j)*cmplx((angle_integral(n_integral+1)-angle_integral(1))/2/pi,0,wp)
                !Use cubic spline curve to interpolate continuously arg (f) and lnf, and then use quadpack functions for integration
                k_integral=2*k-1+j
                S_circle_1_real=integral(fun_circle_real_1,-angle_max,angle_max,epsilon_i)
                S_circle_1_real=r_min**k_integral*S_circle_1_real
                S_circle_1_imag=integral(fun_circle_imag_1,-angle_max,angle_max,epsilon_i)
                S_circle_1_imag=r_min**k_integral*S_circle_1_imag
                
                S_circle_2_real=integral(fun_circle_real_2,-angle_max,angle_max,epsilon_i)
                S_circle_2_imag=integral(fun_circle_imag_2,-angle_max,angle_max,epsilon_i)
                
                angle_integral_signal=-angle_max
                S_line_1_real=integral(fun_line_real_1,r_min,1.0_wp,epsilon_i)
                S_line_1_imag=integral(fun_line_imag_1,r_min,1.0_wp,epsilon_i)
              
                angle_integral_signal=angle_max
                S_line_2_real=integral(fun_line_real_2,r_min,1.0_wp,epsilon_i)
                S_line_2_imag=integral(fun_line_imag_2,r_min,1.0_wp,epsilon_i)
            
                integral_result=cmplx(S_circle_1_real-S_circle_2_real-S_line_1_real+S_line_2_real,S_circle_1_imag-S_circle_2_imag-S_line_1_imag+S_line_2_imag,wp)
                integral_result=integral_result/2/pi*(0,-1)
                integral_result=integral_result+kexi(1)**(2*k-1+j)*(angle_integral(n_integral+1)-angle_integral(1))/2/pi
                epsilon_svd=max(epsilon_svd,abs(integral_result-integral_result_ref))
                
                G(2*k+j)=integral_result
            
            end do
                allocate(hankel_0(k+1,k+1))
                call Hankel_construct(hankel_0,G,k+1,0)
                rank(k+1)=svd_rank(hankel_0,k+1,epsilon_svd,ratio_svd)
                deallocate(hankel_0)
                if(my_id==0)then
                    write(*,"('rank(',i2,'):',i2)"),k+1,rank(k+1)
                end if
                if (k>=2) then
                    if (rank(k+1)==rank(k) .and. rank(k)==rank(k-1))then
                        N=rank(k+1)
                        exit
                    else if (k==N_error) then
                        N=N_error
                    end if
                end if
            
        end do
        deallocate(rank)
        deallocate(angle_integral)
        deallocate(angle_spline)
        deallocate(angle_integral_circle_spline)
        deallocate(lnf_integral_circle_spline)
        deallocate(r_spline)
        deallocate(angle_integral_line_spline)
        deallocate(lnf_integral_line_spline)
        deallocate(kexi)
        deallocate(kexi_m)
        deallocate(z_integral)
        deallocate(zm_integral)
        !When N is not 0 and N,using the QZ algorithm to solve generalized eigenvalue problems for finding zeros and poles 
        if (N==0) then
            condition=0.0_wp
            epsilon_accuracy=0.0_wp
        else if (N==N_error) then
            condition=10*kc_square
            epsilon_accuracy=1.0_wp
        else
            
            allocate(hankel_0(N,N))
            allocate(hankel_1(N,N))
            allocate(ans_kexi(N))
            ans_kexi(:)=(0.0_wp,0.0_wp)
            call Hankel_construct(hankel_0,G,N,0)     
            call Hankel_construct(hankel_1,G,N,1)
            call qz_solve(hankel_1,hankel_0,N,ans_kexi)
            N_effective=N
            epsilon_accuracy=0.0_wp
            !Remove solutions that are not within the arc region
            do j=1,N
                if(abs(ans_kexi(j))<0.99*r_min .or. abs(ans_kexi(j))>1.01 .or. isnan(abs(ans_kexi(j)))) then
                    N_effective=N_effective-1
                end if
            end do
        
            allocate(ans_kexi_effective(N_effective))
            allocate(ans_z(N_effective))
        
            k=1
            do j=1,N
                if(abs(ans_kexi(j))>=0.99*r_min .and. abs(ans_kexi(j))<=1.01) then
                    ans_kexi_effective(k)=ans_kexi(j)
                    ans_z(k)=inverse_tran(ans_kexi(j),z0,L,epsilon_0)
                    k=k+1
                end if 
            end do
            !determine the condition number with ans_kexi_effective
            if(N_effective==N) then
                condition=condition_number(ans_kexi_effective,N_effective)
            else
                condition=10*kc_square
            end if
            !Estimating errors by solving the generalized eigenvalue problem of Hankel matrices of order N+1xN+1
            if(condition<=kc_square) then
                allocate(hankel_0_ref(N+1,N+1))
                allocate(hankel_1_ref(N+1,N+1))
                allocate(ans_kexi_ref(N+1))
		        ans_kexi_ref(:)=(0.0_wp,0.0_wp)
                call Hankel_construct(hankel_0_ref,G,N+1,0)     
                call Hankel_construct(hankel_1_ref,G,N+1,1)
                call qz_solve(hankel_1_ref,hankel_0_ref,N+1,ans_kexi_ref)
                N_effective_ref=N+1
                do j=1,N+1
                    if(abs(ans_kexi_ref(j))<0.99*r_min .or. abs(ans_kexi_ref(j))>1.01 .or. isnan(abs(ans_kexi_ref(j))) ) then
                        N_effective_ref=N_effective_ref-1
                    end if
                end do
            
                allocate(ans_kexi_effective_ref(N_effective_ref))
                allocate(ans_z_ref(N_effective_ref))
                k=1
                do j=1,N+1
                    if (abs(ans_kexi_ref(j))>=0.99*r_min .and. abs(ans_kexi_ref(j))<=1.01) then
                        ans_kexi_effective_ref(k)=ans_kexi_ref(j)
                        ans_z_ref(k)=inverse_tran(ans_kexi_ref(j),z0,L,epsilon_0)
                        k=k+1
                    end if
                end do
            
                if (N_effective>N_effective_ref) then
                    if(my_id==0) then
                        write(*,*) "N+1 order hankel matrix has fewer solution than the N order hankel matrix: there may be some zeros or poles on the boundary during the subdivision"
                    end if
                    condition=10*kc_square
                else
                    do j=1,N_effective
                        ans_z_error(z_solve_number+j)=1
                        do k=1,N_effective_ref
                            ans_z_error(z_solve_number+j)=min(ans_z_error(z_solve_number+j),abs(ans_z(j)-ans_z_ref(k))/2)
                            
                        end do
                        epsilon_accuracy=max(ans_z_error(z_solve_number+j),epsilon_accuracy)
                    end do
                    if (abs(right_edge-left_edge)<=epsilon_accuracy_limit) then 
                        epsilon_accuracy=epsilon_accuracy_limit
                    end if
                    if (epsilon_accuracy<=epsilon_accuracy_limit) then
                        allocate(ans_mul(N_effective))
                        call multiplicities(ans_mul,ans_kexi_effective,G,N_effective)
                        do j=1,N_effective
                            ans_z_solve(z_solve_number+j)=ans_z(j)
                            ans_mul_solve(z_solve_number+j)=ans_mul(j)
                            if(my_id==0) then
			                    write(*,'(A,F15.7,A,F15.7,A,I4,E15.7)'),'(',real(ans_z_solve(z_solve_number+j)),',',aimag(ans_z_solve(z_solve_number+j)),')',ans_mul_solve(z_solve_number+j),ans_z_error(z_solve_number+j)
			                end if
                        end do
                        z_solve_number=z_solve_number+N_effective
                        deallocate(ans_mul)
                    end if

                end if
                deallocate(hankel_0_ref)
                deallocate(hankel_1_ref)
                deallocate(ans_kexi_ref)
                deallocate(ans_kexi_effective_ref)
                deallocate(ans_z_ref)
            end if
            deallocate(hankel_0)
            deallocate(hankel_1)
            deallocate(ans_kexi)
            deallocate(ans_kexi_effective)
            deallocate(ans_z)
        end if
        deallocate(G)
    end subroutine argument_solve
    
    subroutine multiplicities(ans_mul,ans_kexi,G,N)
        implicit none    
        integer,intent(inout)::ans_mul(:)
        complex(wp),intent(in)::ans_kexi(:)
        complex(wp),intent(in)::G(:)
        integer,intent(in)::N
        complex(wp),allocatable::a(:,:),b(:),c(:)
        integer::j,k
        allocate(a(N,N))
        allocate(b(N))
   
        
        do j=1,N
            b(j)=G(j)
            do k=1,N
                a(k,j)=ans_kexi(j)**(k-1)
            end do
        end do
        call gesv(a,b)
        
        
        do j=1,N
            ans_mul(j)=nint(real(b(j)))
        end do
        deallocate(a)
        deallocate(b)
    end subroutine multiplicities
    
    subroutine qz_solve(hankel_1,hankel_0,k,ans_kexi)
        implicit none    
        complex(wp),intent(inout)::hankel_1(:,:)
        complex(wp),intent(inout)::hankel_0(:,:)
        integer,intent(in)::k
        complex(wp),intent(inout)::ans_kexi(:)
        complex(wp),allocatable::alpha(:),beta(:)
        integer::i
        allocate(alpha(k))
        allocate(beta(k))
        call ggev(hankel_1,hankel_0,alpha,beta)
        do i=1,k
            if(abs(beta(i))<1d-12) then
                if(abs(alpha(i))<1d-12) then
                    write(*,*) "there are infinite general eigenvalue"
                
                end if 
            else 
                ans_kexi(i)=alpha(i)/beta(i)
            end if
        end do
        deallocate(alpha)
        deallocate(beta)
    end subroutine
    
    
    subroutine Hankel_construct(hankel,G,k,t)
        implicit none    
        complex(wp),intent(inout)::hankel(:,:)
        complex(wp),intent(in)::G(:)
        integer,intent(in)::k
        integer,intent(in)::t
        integer::m,n
        do m=1,k
            do n=1,k
                hankel(m,n)=G(m+n-1+t)
            end do
        end do

    end subroutine Hankel_construct

    !integral_point_compute:Calculate the discrete node coordinates of the enclosure boundary
    subroutine integral_point_compute(n_circle,n_line,left_edge,right_edge,down_edge,up_edge,angle_max,kexi,kexi_m,z_integral,zm_integral) 
    

        implicit none

        integer,intent(in)::n_circle
        integer,intent(in)::n_line   
        real(wp),intent(in)::left_edge,right_edge,down_edge,up_edge
        real(wp),intent(in)::angle_max
        complex(wp),intent(inout)::kexi(:)
        complex(wp),intent(inout)::kexi_m(:)
        complex(wp),intent(inout)::z_integral(:)
        complex(wp),intent(inout)::zm_integral(:)
        real(wp)::L
        real(wp)::r_min
        integer:: n_integral
        integer::k
        L=abs(right_edge-left_edge)
        r_min=exp(-2*angle_max*abs(up_edge-down_edge)/L)
        n_integral=2*n_circle+2*n_line

    
    
        do k=1,n_circle
            kexi(k)=r_min*exp((0,1)*angle_max/n_circle*(n_circle-2*(k-1)))
            kexi_m(k)=r_min*exp((0,1)*angle_max/n_circle*(n_circle-(2*k-1)));
            z_integral(k)=cmplx((right_edge+left_edge)/2-L/n_circle/2*(n_circle-2*(k-1)),down_edge,wp)
            zm_integral(k)=cmplx((right_edge+left_edge)/2-L/n_circle/2*(n_circle-(2*k-1)),down_edge,wp)
    
        end do
    
        do k=1,n_line
            kexi(k+n_circle)=exp(-2*angle_max*(up_edge-down_edge)/n_line/L*(n_line-(k-1)))*exp((0,-1)*angle_max)
            kexi_m(k+n_circle)=exp(-2*angle_max*(up_edge-down_edge)/n_line/L*(n_line-(k-0.5)))*exp((0,-1)*angle_max)
            z_integral(k+n_circle)=cmplx(right_edge,down_edge+(k-1)*(up_edge-down_edge)/n_line,wp)
            zm_integral(k+n_circle)=cmplx(right_edge,down_edge+(k-0.5)*(up_edge-down_edge)/n_line,wp)
        end do
    
        do k=1,n_circle
            kexi(k+n_circle+n_line)=exp((0,-1)*angle_max/n_circle*(n_circle-2*(k-1)))
            kexi_m(k+n_circle+n_line)=exp((0,-1)*angle_max/n_circle*(n_circle-(2*k-1)))
            z_integral(k+n_circle+n_line)=cmplx((right_edge+left_edge)/2+L/n_circle/2*(n_circle-2*(k-1)),up_edge,wp)
            zm_integral(k+n_circle+n_line)=cmplx((right_edge+left_edge)/2+L/n_circle/2*(n_circle-(2*k-1)),up_edge,wp)
    
        end do
   

    
        do k=1,n_line
            kexi(k+2*n_circle+n_line)=exp(-2*angle_max*(k-1)*(up_edge-down_edge)/n_line/L)*exp((0,1)*angle_max)
            kexi_m(k+2*n_circle+n_line)=exp(-2*angle_max*(k-0.5)*(up_edge-down_edge)/n_line/L)*exp((0,1)*angle_max)
            z_integral(k+2*n_circle+n_line)=cmplx(left_edge,up_edge-(k-1)*(up_edge-down_edge)/n_line,wp)
            zm_integral(k+2*n_circle+n_line)=cmplx(left_edge,up_edge-(k-0.5)*(up_edge-down_edge)/n_line,wp)
        end do
    
    end subroutine integral_point_compute
    
    !arg_trace:keep track of arg(f) to make arg(f) continuous
    subroutine arg_trace(angle_integral,f_integral,n_integral)
        implicit none
        real(wp),intent(inout)::angle_integral(:)
        complex(wp),intent(in)::f_integral(:)
        integer,intent(in)::n_integral
        integer,allocatable::m(:)
        integer::i
        allocate(m(n_integral+1))
        m(:)=0
        do i=1,n_integral
            angle_integral(i)=angle(f_integral(i))
        end do
    
        do i=1,n_integral-1
            if (abs(angle_integral(i+1)-angle_integral(i))>pi) then 
                if (angle_integral(i+1)>angle_integral(i)) then
                    m(i+1)=m(i)-1
                else 
                    m(i+1)=m(i)+1
                end if 
            else 
                m(i+1)=m(i)
            end if 
        end do
    
        if(abs(angle_integral(1)-angle_integral(n_integral))>pi) then
            if (angle_integral(1)>angle_integral(n_integral)) then
                m(n_integral+1)=m(n_integral)-1
            else
                m(n_integral+1)=m(n_integral)+1
            end if 
        else
            m(n_integral+1)=m(n_integral)   
        end if
    
        do i=1,n_integral
            angle_integral(i)=angle_integral(i)+2*pi*m(i)
    
        end do
        angle_integral(n_integral+1)=angle_integral(1)+2*pi*m(n_integral+1)
        deallocate(m)
    end subroutine arg_trace
    
    !svd_rank:use gesvd from lapack library to perform svd decomposition of the input matrix and determine the rank of the matrix based on singular values
    function svd_rank(hankel,k,epsilon_svd,ratio_svd) result(ans)
        implicit none
        complex(wp),intent(inout)::hankel(:,:)
        integer,intent(in)::k
        real(wp),intent(in)::epsilon_svd
        real(wp),intent(in)::ratio_svd
        real(wp),allocatable::s(:)
        integer::ans,n
        allocate(s(k))
        call gesvd(hankel,s)
        if (s(1)<epsilon_svd) then
            ans=0
        else
            do n=2,k
                if(s(n-1)/s(n)>ratio_svd .or. s(n)<epsilon_svd) then
                    ans=n-1
                    exit
                else if (n==k) then
                    ans=k
                end if
            end do
        end if
        !do n=1,k
        !    write(*,*),s(n)
        !end do
        deallocate(s)
    end function svd_rank
    
    !inetgral:use dqag function from quadpack library to do the integral of input function over a and b
    function integral(h,a,b,epsabs) result (ans)
        implicit none
        procedure(func_1) :: h 
        real(wp), intent(in):: a 
        real(wp), intent(in) :: b 
        real(wp),intent(in)::epsabs
        real(wp)::epsrel
        integer,parameter:: key = 6
        integer,parameter:: limit = 10000
        integer,parameter:: lenw=limit*4
        real(wp) :: abserr, ans, work(lenw)
        integer :: ier, iwork(limit), last, neval
       
      
        epsrel=epsabs
        call dqag(h, a, b, epsabs, epsrel, key, ans, &
                    abserr, neval, ier, limit, lenw, last, &
                    iwork, work)
            
    end function integral
    !fun_circle_real_1,fun_circle_imag_1...:the real/imag part of the contour integral in the circle/line part
    real(wp) function fun_circle_real_1(x)
        implicit none
        real(wp),intent(in)::x
        fun_circle_real_1=-k_integral*(sin(k_integral*x)*spline_evaluate(x,spl_lnf_circle_1)+cos(k_integral*x)*spline_evaluate(x,spl_angle_circle_1))
        
    end function fun_circle_real_1
    
    
    real(wp) function fun_circle_imag_1(x)
        implicit none
        real(wp),intent(in)::x
        fun_circle_imag_1=k_integral*(cos(k_integral*x)*spline_evaluate(x,spl_lnf_circle_1)-sin(k_integral*x)*spline_evaluate(x,spl_angle_circle_1))
    
    end function fun_circle_imag_1
   
      real(wp) function fun_circle_real_2(x)
        implicit none
        real(wp),intent(in)::x
        fun_circle_real_2=-k_integral*(sin(k_integral*x)*spline_evaluate(x,spl_lnf_circle_2)+cos(k_integral*x)*spline_evaluate(x,spl_angle_circle_2))
        
    end function fun_circle_real_2
    
    
    real(wp) function fun_circle_imag_2(x)
        implicit none
        real(wp),intent(in)::x
        fun_circle_imag_2=k_integral*(cos(k_integral*x)*spline_evaluate(x,spl_lnf_circle_2)-sin(k_integral*x)*spline_evaluate(x,spl_angle_circle_2))
    
    end function fun_circle_imag_2
    
    real(wp) function fun_line_real_1(x)
        implicit none
        real(wp),intent(in)::x
        fun_line_real_1=k_integral*x**(k_integral-1)*(cos(k_integral*angle_integral_signal)*spline_evaluate(x,spl_lnf_line_1)-sin(k_integral*angle_integral_signal)*spline_evaluate(x,spl_angle_line_1))   
    end function fun_line_real_1
    
    real(wp) function fun_line_imag_1(x)
        implicit none
        real(wp),intent(in)::x
        fun_line_imag_1=k_integral*x**(k_integral-1)*(sin(k_integral*angle_integral_signal)*spline_evaluate(x,spl_lnf_line_1)+cos(k_integral*angle_integral_signal)*spline_evaluate(x,spl_angle_line_1))
    end function fun_line_imag_1
    
    real(wp) function fun_line_real_2(x)
        implicit none
        real(wp),intent(in)::x
        fun_line_real_2=k_integral*x**(k_integral-1)*(cos(k_integral*angle_integral_signal)*spline_evaluate(x,spl_lnf_line_2)-sin(k_integral*angle_integral_signal)*spline_evaluate(x,spl_angle_line_2))   
    end function fun_line_real_2
    
    real(wp) function fun_line_imag_2(x)
        implicit none
        real(wp),intent(in)::x
        fun_line_imag_2=k_integral*x**(k_integral-1)*(sin(k_integral*angle_integral_signal)*spline_evaluate(x,spl_lnf_line_2)+cos(k_integral*angle_integral_signal)*spline_evaluate(x,spl_angle_line_2))
    end function fun_line_imag_2
    
    !condition_number: caculate the condition number from ans_kexi and N
    function condition_number(ans_kexi,N) result(condition)
        implicit none
        complex(wp),intent(in)::ans_kexi(:)
        integer,intent(in)::N
        real(wp)::condition
        real(wp)::r_max,r_min,multiply,multiply_max
        integer::j,k
        r_max=abs(ans_kexi(1))
        r_min=abs(ans_kexi(1))
        do j=1,N
            r_max=max(r_max,abs(ans_kexi(j)))
            r_min=min(r_min,abs(ans_kexi(j)))
            multiply=1.0_wp
            do k=1,N
                if (k /=j) then
                    multiply=multiply/abs(ans_kexi(j)/abs(ans_kexi(j))-ans_kexi(k)/abs(ans_kexi(k)))
                end if
            end do
        
            if (j==1) then
                multiply_max=multiply
            else 
                multiply_max=max(multiply_max,multiply)
            end if
        end do
        
        multiply_max=multiply_max**2
        condition=N**2*(r_max/r_min)**(2*N-2)*(1+1.0_wp/r_max)**(2*N-2)*multiply_max
   
    end function condition_number
    
    !angle:calculate the angle of a complex number
    function angle(kexi) result(a)
        implicit none
        complex(wp),intent(in)::kexi
        real(wp)::cos_kexi
        real(wp)::a
        real(wp)::r
   
        r=sqrt(real(kexi)**2+imag(kexi)**2)
        cos_kexi=real(kexi)/r
    
        if (imag(kexi)<0) then
            a=-acos(cos_kexi)
        else
            a=acos(cos_kexi)
        end if 
    
    end function angle
    
    !inverse_tran:inverse transform the coordinate from kexi_space to z_space
    function inverse_tran(kexi,z0,L,epsilon_0) result(ans)
        implicit none
        complex(wp),intent(in)::kexi
        complex(wp),intent(in)::z0
        real(wp),intent(in)::L
        real(wp),intent(in)::epsilon_0
        complex(wp)::ans

        ans=(0,1)*L*log(kexi)/(2*pi-epsilon_0)+z0
    
    end function inverse_tran



end module zpl5
    
