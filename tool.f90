 
 
 
 module tool

 use iso_fortran_env,only:wp=>real64
 use lapack95
 use mathfunction
 
 
 implicit none
 private c_div_v_para,omega_div_omega_ci,omega_div_omega_ce,omega_pi_div_omega_ci,omega_pe_div_omega_ce,ratio_ti,ratio_te,refractive_para
 private k_para_rho_i_para,k_per_rho_i_para,k_per_rho_i_per,k_para_rho_e_para,k_para_rho_e_per,k_per_rho_e_para,k_per_rho_e_per
 private omega_pd_div_omega_cd, k_para_c_div_omega_ci,k_per_c_div_omega_ci
 real(wp)::c_div_v_para,omega_div_omega_ci,omega_div_omega_ce,omega_pi_div_omega_ci,omega_pe_div_omega_ce,k_per_rho_e_per,ratio_ti,ratio_te,refractive_para
 real(wp)::k_para_rho_i_para,k_per_rho_i_para,k_per_rho_i_per,k_para_rho_e_para,k_para_rho_e_per,k_per_rho_e_para
 real(wp):: omega_pd_div_omega_cd, k_para_c_div_omega_ci,k_per_c_div_omega_ci
contains

!-----------------------------------------------------------------------------!
!     test_function_1
!-----------------------------------------------------------------------------!
    complex(wp) function test_fun(x)
	implicit none
	complex(wp),intent(in)::x
	!test_fun=1/(x-(1.0_wp,1.0_wp))/(x-(1.0_wp,1.0_wp))+1/(x-(1.0_wp,1.0_wp))
	real(wp)::lammda,miu,P0
        lammda=1.0_wp
        miu=1.0_wp
        P0=1.0_wp
        test_fun=lammda*exp(-x**2)-0.25*(cmplx(-2,miu,wp)*x+P0)**2-cmplx(1.0_wp,0.5_wp*miu,wp)
        !test_fun=log(test_fun)
    end function test_fun  
    
!-----------------------------------------------------------------------------!
!     test_function_2
!-----------------------------------------------------------------------------!    
      complex(wp) function test_fun2(x)
        implicit none
        complex(wp), intent(in) :: x
        real(wp)::A0(3,3),A1(3,3),A2(3,3)
        complex(wp)::D(3,3)
        complex(wp)::determinant=(0.0_wp,0.0_wp)
  
        integer::i,j
        A2=reshape([17.6_wp,1.28_wp,2.89_wp,1.28_wp,0.824_wp,0.413_wp,2.89_wp,0.413_wp,0.725_wp],[3,3])
        A1=reshape([7.66_wp,2.45_wp,2.1_wp,0.23_wp,1.04_wp,0.223_wp,0.6_wp,0.756_wp,0.658_wp],[3,3])
        A0=reshape([12.1_wp,18.9_wp,15.9_wp,0.0_wp,2.7_wp,0.145_wp,11.9_wp,3.64_wp,15.5_wp],[3,3])
        
        do i=1,3
            do j=1,3
            D(i,j)=A2(i,j)*(exp(x)-1)+x**2*A1(i,j)-A0(i,j)
            end do
        end do
	
        determinant=D(1,1)*(D(2,2)*D(3,3)-D(2,3)*D(3,2))-D(1,2)*(D(2,1)*D(3,3)-D(2,3)*D(3,1))+D(1,3)*(D(2,1)*D(3,2)-D(2,2)*D(3,1))
        
        
        test_fun2=determinant
        !test_fun2=log(test_fun2)
    end function test_fun2
!-----------------------------------------------------------------------------!
!     test_function_3_plasma_dispersion_function
!-----------------------------------------------------------------------------!       
    complex(wp) function test_fun3(x)
	implicit none
	complex(wp),intent(in)::x
	call pdf(x,test_fun3)
    !test_fun3=log(test_fun3)
    end function test_fun3
!-----------------------------------------------------------------------------!
!     set_parameter_old_version: given k to calculate omega (old version)
!-----------------------------------------------------------------------------!    
    subroutine set_parameter_old_version(b,n,ti_para,ti_per,te_para,te_per,k_para,k_per)
    implicit none
    real(wp),intent(in)::b,n,ti_para,ti_per,te_para,te_per,k_para,k_per
    real(wp)::c,e,mi,me,epsilon_0,vi_para,vi_per,ve_para,ve_per
    real(wp)::omega_pi,omega_pe,omega_ci,omega_ce,rho_i_para,rho_i_per,rho_e_para,rho_e_per
    c=3d8
    mi=1.672621777d-27
    e=1.602176565d-19
    me=mi/1836
    epsilon_0=8.854187817d-12
    ratio_ti=ti_per/ti_para
    ratio_te=te_per/te_para
    vi_para=(2*e*ti_para/mi)**0.5
    vi_per=(2*e*ti_per/mi)**0.5
    ve_para=(2*e*te_para/me)**0.5
    ve_per=(2*e*te_per/me)**0.5
    omega_pi=(n*e**2/mi/epsilon_0)**0.5
    omega_pe=(n*e**2/me/epsilon_0)**0.5
    omega_ci=e*b/mi
    omega_ce=-e*b/me
    rho_i_para=vi_para/omega_ci
    rho_i_per=vi_per/omega_ci
    rho_e_para=ve_para/omega_ce
    rho_e_per=ve_per/omega_ce
  
    c_div_v_para=c**2/vi_para**2
    omega_pi_div_omega_ci=omega_pi**2/omega_ci**2
    omega_pe_div_omega_ce=omega_pe**2/omega_ce**2
    k_para_rho_i_para=k_para*rho_i_para
    k_per_rho_i_para=k_per*rho_i_para
    k_per_rho_i_per=k_per*rho_i_per
    k_para_rho_e_para=k_para*rho_e_para
    k_para_rho_e_per=k_para*rho_e_per
    k_per_rho_e_per=k_per*rho_e_per
    k_per_rho_e_para=k_per*rho_e_para
    end subroutine set_parameter_old_version

!-----------------------------------------------------------------------------!
!     set_parameter: given k to calculate omega (new version)
!-----------------------------------------------------------------------------!    
    subroutine set_parameter(c_div_v_para_input,omega_pe_div_omega_ce_input,k_para_rho_i_para_input,k_para_rho_e_para_input,k_para_rho_e_per_input,k_per_rho_i_para_input,k_per_rho_i_per_input,k_per_rho_e_para_input,k_per_rho_e_per_input)
		implicit none
		real(wp),intent(in)::c_div_v_para_input,omega_pe_div_omega_ce_input,k_para_rho_i_para_input,k_para_rho_e_para_input,k_para_rho_e_per_input,k_per_rho_i_para_input,k_per_rho_i_per_input,k_per_rho_e_para_input,k_per_rho_e_per_input
		
		c_div_v_para=c_div_v_para_input
		omega_pe_div_omega_ce=omega_pe_div_omega_ce_input
		omega_pi_div_omega_ci=omega_pe_div_omega_ce_input*1836
		
		k_para_rho_i_para=k_para_rho_i_para_input
		k_per_rho_i_para=k_per_rho_i_para_input
		k_per_rho_i_per=k_per_rho_i_per_input
		k_para_rho_e_para=k_para_rho_e_para_input
		k_para_rho_e_per=k_para_rho_e_per_input
		k_per_rho_e_per=k_per_rho_e_per_input
		k_per_rho_e_para=k_per_rho_e_para_input
		end subroutine set_parameter
		
!-----------------------------------------------------------------------------!
!     set_parameter_cold: given k to calculate omega (cold plasma version)
!-----------------------------------------------------------------------------!    
		subroutine set_parameter_cold(omega_pe_div_omega_ce_input,k_para_c_div_omega_ci_input,k_per_c_div_omega_ci_input)
			implicit none
			real(wp),intent(in)::omega_pe_div_omega_ce_input,k_para_c_div_omega_ci_input,k_per_c_div_omega_ci_input
			
			omega_pe_div_omega_ce=omega_pe_div_omega_ce_input
			omega_pi_div_omega_ci=omega_pe_div_omega_ce_input*1836/2
			omega_pd_div_omega_cd=omega_pe_div_omega_ce_input*1836
			k_para_c_div_omega_ci=k_para_c_div_omega_ci_input
			k_per_c_div_omega_ci=k_per_c_div_omega_ci_input
			end subroutine set_parameter_cold

!-----------------------------------------------------------------------------!
!     set_parameter_variable_k_para:given omega k_perp to calculate k_para
!-----------------------------------------------------------------------------!    
    subroutine set_parameter_variable_k_para(c_div_v_para_input,omega_pe_div_omega_ce_input,omega_div_omega_ci_input,k_per_rho_i_input)
    implicit none
    real(wp),intent(in)::c_div_v_para_input,omega_pe_div_omega_ce_input,omega_div_omega_ci_input,k_per_rho_i_input
	real(wp)::mass_ratio
	mass_ratio=1836
    c_div_v_para=c_div_v_para_input
    omega_div_omega_ci=omega_div_omega_ci_input
    omega_div_omega_ce=-omega_div_omega_ci_input/mass_ratio
    omega_pe_div_omega_ce=omega_pe_div_omega_ce_input
    omega_pi_div_omega_ci=omega_pe_div_omega_ce_input*mass_ratio
	
    k_per_rho_i_para=k_per_rho_i_input
    k_per_rho_i_per=k_per_rho_i_input
    k_per_rho_e_per=-k_per_rho_i_input/(mass_ratio)**(0.5)
    k_per_rho_e_para=-k_per_rho_i_input/(mass_ratio)**(0.5)
    end subroutine set_parameter_variable_k_para
    
!-----------------------------------------------------------------------------!
!     set_parameter_variable_k_per:given omega and k_para to calculate k_perp
!-----------------------------------------------------------------------------!    
    subroutine set_parameter_variable_k_per(c_div_v_para_input,omega_pe_div_omega_ce_input,omega_div_omega_ci_input,k_para_rho_i_input)
    implicit none
    real(wp),intent(in)::c_div_v_para_input,omega_pe_div_omega_ce_input,omega_div_omega_ci_input,k_para_rho_i_input
	real(wp)::mass_ratio
	mass_ratio=1836
    c_div_v_para=c_div_v_para_input
    omega_div_omega_ci=omega_div_omega_ci_input
    omega_div_omega_ce=-omega_div_omega_ci_input/mass_ratio
    omega_pe_div_omega_ce=omega_pe_div_omega_ce_input
    omega_pi_div_omega_ci=omega_pe_div_omega_ce_input*mass_ratio
	k_para_rho_i_para=k_para_rho_i_input
	k_para_rho_e_para=-k_para_rho_i_input/(mass_ratio)**(0.5)

    end subroutine set_parameter_variable_k_per   
    

!-----------------------------------------------------------------------------!
!     set_parameter_cold_variable_k_per:cold plasma dispersion relation, given n_para and omega to calculate n_perp
!-----------------------------------------------------------------------------!    
    subroutine set_parameter_cold_variable_k_per(omega_pe_div_omega_ce_input,omega_div_omega_ci_input,refractive_para_input)
    implicit none
    real(wp),intent(in)::omega_pe_div_omega_ce_input,omega_div_omega_ci_input,refractive_para_input
	real(wp)::mass_ratio
	mass_ratio=1836
    omega_div_omega_ci=omega_div_omega_ci_input
    omega_div_omega_ce=-omega_div_omega_ci_input/mass_ratio
    omega_pe_div_omega_ce=omega_pe_div_omega_ce_input
    omega_pi_div_omega_ci=omega_pe_div_omega_ce_input*mass_ratio
	refractive_para=refractive_para_input

    end subroutine set_parameter_cold_variable_k_per   
    
    
!-----------------------------------------------------------------------------!
!     set_parameter2:given k_perp to calculate omega for perpendicular wave
!-----------------------------------------------------------------------------! 
    subroutine set_parameter2(omega_pe_div_omega_ce_sqrt,k_per_rho_e_per_input)
    implicit none
    real(wp),intent(in)::omega_pe_div_omega_ce_sqrt,k_per_rho_e_per_input
    omega_pe_div_omega_ce=omega_pe_div_omega_ce_sqrt**2
    k_per_rho_e_per=k_per_rho_e_per_input
    omega_pi_div_omega_ci=omega_pe_div_omega_ce*1836
    k_per_rho_i_per=k_per_rho_e_per_input*(1836)**0.5
    end subroutine set_parameter2

!-----------------------------------------------------------------------------!
!     set_parameter_special:given k to calculate omega. The density input is replaced by omega_pe_div_omega_ce_input
!-----------------------------------------------------------------------------!
    
    subroutine set_parameter_special(b,ti_para,ti_per,te_para,te_per,k_para,k_per,omega_pe_div_omega_ce_input)
    implicit none
    real(wp),intent(in)::b,ti_para,ti_per,te_para,te_per,k_para,k_per,omega_pe_div_omega_ce_input
    real(wp)::c,e,mi,me,epsilon_0,vi_para,vi_per,ve_para,ve_per
    real(wp)::omega_ci,omega_ce,rho_i_para,rho_i_per,rho_e_para,rho_e_per
    c=3d8
    mi=1.672621777d-27
    e=1.602176565d-19
    me=mi/1836
    epsilon_0=8.854187817d-12
    ratio_ti=ti_per/ti_para
    ratio_te=te_per/te_para
    vi_para=(2*e*ti_para/mi)**0.5
    vi_per=(2*e*ti_per/mi)**0.5
    ve_para=(2*e*te_para/me)**0.5
    ve_per=(2*e*te_per/me)**0.5
    omega_ci=e*b/mi
    omega_ce=-e*b/me
    rho_i_para=vi_para/omega_ci
    rho_i_per=vi_per/omega_ci
    rho_e_para=ve_para/omega_ce
    rho_e_per=ve_per/omega_ce
  
    c_div_v_para=c**2/vi_para**2
    k_para_rho_i_para=k_para*rho_i_para
    k_per_rho_i_para=k_per*rho_i_para
    k_per_rho_i_per=k_per*rho_i_per
    k_para_rho_e_para=k_para*rho_e_para
    k_para_rho_e_per=k_para*rho_e_per
    k_per_rho_e_per=k_per*rho_e_per
    k_per_rho_e_para=k_per*rho_e_para
    omega_pe_div_omega_ce=omega_pe_div_omega_ce_input
    omega_pi_div_omega_ci=omega_pe_div_omega_ce*1836
    end subroutine set_parameter_special
!-----------------------------------------------------------------------------!
!     dispersion_function_parallel_matrix: dispersion matrix of omega for parallel waves
!-----------------------------------------------------------------------------!
    subroutine dispersion_function_parallel_matrix(x,D)
    implicit none
    complex(wp),intent(in)::x
    complex(wp),intent(out)::D(:,:)
    complex(wp)::xi_pdf_1,yi_pdf_1,xi_pdf_2,yi_pdf_2,x_e
    complex(wp)::xe_pdf_1,ye_pdf_1,xe_pdf_2,ye_pdf_2
	x_e=-x/1836
    ratio_ti=1
	ratio_te=1
	xi_pdf_1=(x-(1,0))/k_para_rho_i_para
	call pdf(xi_pdf_1,yi_pdf_1)
	xi_pdf_2=(x+(1,0))/k_para_rho_i_para
	call pdf(xi_pdf_2,yi_pdf_2)
    
    xe_pdf_1=(x_e-(1,0))/k_para_rho_e_para
	call pdf(xe_pdf_1,ye_pdf_1)
	xe_pdf_2=(x_e+(1,0))/k_para_rho_e_para
	call pdf(xe_pdf_2,ye_pdf_2)
    
    
	D(1,1)=1-c_div_v_para*k_para_rho_i_para**2/x**2+0.5*omega_pi_div_omega_ci/x*(yi_pdf_1/k_para_rho_i_para+(ratio_ti-1)*(1+xi_pdf_1*yi_pdf_1)/x)&
	&+0.5*omega_pi_div_omega_ci/x*(yi_pdf_2/k_para_rho_i_para+(ratio_ti-1)*(1+xi_pdf_2*yi_pdf_2)/x)&
    &+0.5*omega_pe_div_omega_ce/x_e*(ye_pdf_1/k_para_rho_e_para+(ratio_te-1)*(1+xe_pdf_1*ye_pdf_1)/x)&
    &+0.5*omega_pe_div_omega_ce/x_e*(ye_pdf_2/k_para_rho_e_para+(ratio_te-1)*(1+xe_pdf_2*ye_pdf_2)/x)  
	D(1,2)=0.5*(0,1)*omega_pi_div_omega_ci/x*(yi_pdf_1/k_para_rho_i_para+(ratio_ti-1)*(1+xi_pdf_1*yi_pdf_1)/x-yi_pdf_2/k_para_rho_i_para-(ratio_ti-1)*(1+xi_pdf_2*yi_pdf_2)/x)&
    &+0.5*(0,1)*omega_pe_div_omega_ce/x_e*(ye_pdf_1/k_para_rho_e_para+(ratio_te-1)*(1+xe_pdf_1*ye_pdf_1)/x_e-ye_pdf_2/k_para_rho_e_para-(ratio_te-1)*(1+xe_pdf_2*ye_pdf_2)/x_e)
	D(1,3)=0
	D(2,1)=-D(1,2)
	D(2,2)=D(1,1)
	D(2,3)=0
	D(3,1)=0
	D(3,2)=0
	xi_pdf_1=x/k_para_rho_i_para
	xe_pdf_1=x_e/k_para_rho_e_para
	call pdf(xi_pdf_1,yi_pdf_1)
	call pdf(xe_pdf_1,ye_pdf_1) 
	D(3,3)=1+2*omega_pi_div_omega_ci/k_para_rho_i_para**2*(1+xi_pdf_1*yi_pdf_1)+2*omega_pe_div_omega_ce/k_para_rho_e_para**2*(1+xe_pdf_1*ye_pdf_1)


    
    end subroutine dispersion_function_parallel_matrix

!-----------------------------------------------------------------------------!
!     dispersion_function_parallel_matrix_variable_k_para: the dispersion matrix of k_para with known omega for parallel waves.
!-----------------------------------------------------------------------------!
    subroutine dispersion_function_parallel_matrix_variable_k_para(x,D)
    implicit none
    complex(wp),intent(in)::x
    complex(wp),intent(out)::D(:,:)
    complex(wp)::xi_pdf_1,yi_pdf_1,xi_pdf_2,yi_pdf_2,x_e
    complex(wp)::xe_pdf_1,ye_pdf_1,xe_pdf_2,ye_pdf_2
    real(wp)::mass_ratio
    mass_ratio=1836
	x_e=-x/(mass_ratio)**(0.5)
    
	xi_pdf_1=(omega_div_omega_ci-(1,0))/x
	call pdf(xi_pdf_1,yi_pdf_1)
	xi_pdf_2=(omega_div_omega_ci+(1,0))/x
	call pdf(xi_pdf_2,yi_pdf_2)
    
    xe_pdf_1=(omega_div_omega_ce-(1,0))/x_e
	call pdf(xe_pdf_1,ye_pdf_1)
	xe_pdf_2=(omega_div_omega_ce+(1,0))/x_e
	call pdf(xe_pdf_2,ye_pdf_2)
    
    
	D(1,1)=1-c_div_v_para*x**2/omega_div_omega_ci**2+0.5*omega_pi_div_omega_ci/omega_div_omega_ci*(yi_pdf_1/x)&
	&+0.5*omega_pi_div_omega_ci/omega_div_omega_ci*(yi_pdf_2/x)&
    &+0.5*omega_pe_div_omega_ce/omega_div_omega_ce*(ye_pdf_1/x_e)&
    &+0.5*omega_pe_div_omega_ce/omega_div_omega_ce*(ye_pdf_2/x_e)  
	D(1,2)=0.5*(0,1)*omega_pi_div_omega_ci/omega_div_omega_ci*(yi_pdf_1/x-yi_pdf_2/x)&
    &+0.5*(0,1)*omega_pe_div_omega_ce/omega_div_omega_ce*(ye_pdf_1/x_e-ye_pdf_2/x_e)
	D(1,3)=0
	D(2,1)=-D(1,2)
	D(2,2)=D(1,1)
	D(2,3)=0
	D(3,1)=0
	D(3,2)=0
	xi_pdf_1=omega_div_omega_ci/x
	xe_pdf_1=omega_div_omega_ce/x_e
	call pdf(xi_pdf_1,yi_pdf_1)
	call pdf(xe_pdf_1,ye_pdf_1) 
	D(3,3)=1+2*omega_pi_div_omega_ci/x**2*(1+xi_pdf_1*yi_pdf_1)+2*omega_pe_div_omega_ce/x_e**2*(1+xe_pdf_1*ye_pdf_1)
    
    end subroutine dispersion_function_parallel_matrix_variable_k_para
    
!-----------------------------------------------------------------------------!
!     dispersion_function_matrix:the dispersion matrix for omega with known k.
!-----------------------------------------------------------------------------!
    subroutine dispersion_function_matrix(x,D)
    implicit none
    complex(wp),intent(in)::x
    complex(wp),intent(out)::D(:,:)
    complex(wp)::series_sum_1,series_sum_2,x_pdf,y_pdf
    complex(wp)::xe_pdf_1,xe_pdf_2,ye_pdf_1,ye_pdf_2
    integer::n,k
    real(wp)::gamma_n,gamma_n1
    complex(wp)::x_e
	
	x_e=-x/1836
	do k=1,9
	    select case(k)
	    case(1)
		series_sum_1=(0,0)
		series_sum_2=(0,0)
		do n=1,10
		    call bessel_gn(k_per_rho_i_per**2/2,n, gamma_n)
		    call bessel_gn(k_per_rho_e_per**2/2,n, gamma_n1)
		    if(abs(k_para_rho_i_para)<1e-6) then
		      !series_sum_1=series_sum_1+0.5*gamma_n*(n**2/k_per_rho_i_per**2/((cmplx(n,0)-x)))
		      !series_sum_1=series_sum_1+0.5*gamma_n*(n**2/k_per_rho_i_per**2/((cmplx(-n,0)-x)))
		      series_sum_1=series_sum_1+4/k_per_rho_i_per**2*omega_pi_div_omega_ci*gamma_n*(n**2/(cmplx(n**2,0)-x**2))
		      series_sum_1=series_sum_1+4/k_per_rho_e_per**2*omega_pe_div_omega_ce*gamma_n1*(n**2/(cmplx(n**2,0)-x_e**2))
		    else
		      x_pdf=(x-cmplx(n,0))/k_para_rho_i_para
		      call pdf(x_pdf,y_pdf)
		      
		      series_sum_1=series_sum_1+2*omega_pi_div_omega_ci/x*gamma_n*(n**2/k_per_rho_i_per**2/k_para_rho_i_para*y_pdf+n**2/x*(1/k_per_rho_i_para**2-1/k_per_rho_i_per**2)*(1+x_pdf*y_pdf))
		      
		      x_pdf=(x+cmplx(n,0))/k_para_rho_i_para
		      call pdf(x_pdf,y_pdf)
		    
		      series_sum_1=series_sum_1+2*omega_pi_div_omega_ci/x*gamma_n*(n**2/k_per_rho_i_per**2/k_para_rho_i_para*y_pdf+n**2/x*(1/k_per_rho_i_para**2-1/k_per_rho_i_per**2)*(1+x_pdf*y_pdf))
		      
		      
		      x_pdf=(x_e-cmplx(n,0))/k_para_rho_e_para
		      call pdf(x_pdf,y_pdf)
		      
		      series_sum_1=series_sum_1+2*omega_pe_div_omega_ce/x_e*gamma_n1*(n**2/k_per_rho_e_per**2/k_para_rho_e_para*y_pdf+n**2/x_e*(1/k_per_rho_e_para**2-1/k_per_rho_e_per**2)*(1+x_pdf*y_pdf))
		      
		      x_pdf=(x_e+cmplx(n,0))/k_para_rho_e_para
		      call pdf(x_pdf,y_pdf)
		    
		      series_sum_1=series_sum_1+2*omega_pe_div_omega_ce/x_e*gamma_n1*(n**2/k_per_rho_e_per**2/k_para_rho_e_para*y_pdf+n**2/x_e*(1/k_per_rho_e_para**2-1/k_per_rho_e_per**2)*(1+x_pdf*y_pdf))
		      
		    end if
		    
		    if(n>=10 .and. abs(series_sum_1-series_sum_2)<1e-8*abs(series_sum_1+series_sum_2) )then
			
			exit
		    else if(n==50) then
			write(*,*) 'The sum of the first 50 series does not converge'
		    end if
		    series_sum_2=series_sum_1
		end do
		
		if (abs(k_para_rho_i_para)<1e-6) then
		    !D(1,1)=(1,0)+omega_pe_div_omega_ce/(1-x_e**2)+4*omega_pi_div_omega_ci/x*series_sum_1
		    D(1,1)=(1,0)+series_sum_1
		else
		    D(1,1)=(1,0)-c_div_v_para*k_para_rho_i_para**2/x**2+series_sum_1
		end if
        
	  case(2)
		series_sum_1=(0,0)
		series_sum_2=(0,0)
		do n=1,10
		    
		    if(abs(k_para_rho_i_para)<1e-6) then
		      call bessel_gn(k_per_rho_i_per**2/2,n, gamma_n)
		      call bessel_gn(k_per_rho_i_per**2/2,n+1,gamma_n1)
		      !series_sum_1=series_sum_1+0.25*(k_per_rho_i_per**2*(gamma_n1-gamma_n)+2*n*gamma_n)*(n/k_per_rho_i_per**2/(cmplx(n,0)-x))
		      !series_sum_1=series_sum_1+0.25*(k_per_rho_i_per**2*(gamma_n1-gamma_n)+2*n*gamma_n)*(-n/k_per_rho_i_per**2/(cmplx(-n,0)-x))
		     
		      series_sum_1=series_sum_1+(0,1)*omega_pi_div_omega_ci/x*(2*n**2)/(cmplx(n**2,0)-x**2)*(gamma_n1-(1-2*n/k_per_rho_i_per**2)*gamma_n)
		      call bessel_gn(k_per_rho_e_per**2/2,n, gamma_n)
		      call bessel_gn(k_per_rho_e_per**2/2,n+1,gamma_n1)
		      series_sum_1=series_sum_1+(0,1)*omega_pe_div_omega_ce/x_e*(2*n**2)/(cmplx(n**2,0)-x_e**2)*(gamma_n1-(1-2*n/k_per_rho_e_per**2)*gamma_n)
		    else
		      call bessel_gn(k_per_rho_i_per**2/2,n, gamma_n)
		      call bessel_gn(k_per_rho_i_per**2/2,n+1,gamma_n1)
		      x_pdf=(x-cmplx(n,0))/k_para_rho_i_para
		      call pdf(x_pdf,y_pdf)
		      
		      series_sum_1=series_sum_1+cmplx(0,omega_pi_div_omega_ci,wp)/x*(k_per_rho_i_per**2*(gamma_n1-gamma_n)+2*n*gamma_n)*(n/k_per_rho_i_per**2/k_para_rho_i_para*y_pdf+n/x*(1/k_per_rho_i_para**2-1/k_per_rho_i_per**2)*(1+x_pdf*y_pdf))
		    
		    
		      x_pdf=(x+cmplx(n,0))/k_para_rho_i_para
		      call pdf(x_pdf,y_pdf)
		   
		      series_sum_1=series_sum_1-cmplx(0,omega_pi_div_omega_ci,wp)/x*(k_per_rho_i_per**2*(gamma_n1-gamma_n)+2*n*gamma_n)*(n/k_per_rho_i_per**2/k_para_rho_i_para*y_pdf+n/x*(1/k_per_rho_i_para**2-1/k_per_rho_i_per**2)*(1+x_pdf*y_pdf))
		    
		      call bessel_gn(k_per_rho_e_per**2/2,n, gamma_n)
		      call bessel_gn(k_per_rho_e_per**2/2,n+1,gamma_n1)
		      x_pdf=(x_e-cmplx(n,0))/k_para_rho_e_para
		      call pdf(x_pdf,y_pdf)
		      
		      series_sum_1=series_sum_1+cmplx(0,omega_pe_div_omega_ce,wp)/x_e*(k_per_rho_e_per**2*(gamma_n1-gamma_n)+2*n*gamma_n)*(n/k_per_rho_e_per**2/k_para_rho_e_para*y_pdf+n/x_e*(1/k_per_rho_e_para**2-1/k_per_rho_e_per**2)*(1+x_pdf*y_pdf))
		    
		    
		      x_pdf=(x_e+cmplx(n,0))/k_para_rho_e_para
		      call pdf(x_pdf,y_pdf)
		   
		      series_sum_1=series_sum_1-cmplx(0,omega_pe_div_omega_ce,wp)/x_e*(k_per_rho_e_per**2*(gamma_n1-gamma_n)+2*n*gamma_n)*(n/k_per_rho_e_per**2/k_para_rho_e_para*y_pdf+n/x_e*(1/k_per_rho_e_para**2-1/k_per_rho_e_per**2)*(1+x_pdf*y_pdf))
		    
		      
		    end if
		    
		    if(n>=10 .and. abs(series_sum_1-series_sum_2)<1e-8*abs(series_sum_1+series_sum_2) )then
			
			exit
		    else if(n==50) then
			write(*,*) 'The sum of the first 50 series does not converge'
		    end if
		    series_sum_2=series_sum_1
		end do
		if (abs(k_para_rho_i_para)<1e-6) then
			D(1,2)=series_sum_1
		else
		    D(1,2)=series_sum_1
		end if
          
	    case(3)
		series_sum_1=(0,0)
		series_sum_2=(0,0)
		if(abs(k_para_rho_i_para)<1e-6) then
		  D(1,3)=(0,0)
		else
		   do n=1,10
		      call bessel_gn(k_per_rho_i_per**2/2,n, gamma_n)
		    
		      x_pdf=(x-cmplx(n,0))/k_para_rho_i_para
		      call pdf(x_pdf,y_pdf)
		      series_sum_1=series_sum_1+2*omega_pi_div_omega_ci/x*n*gamma_n*k_per_rho_i_para/k_para_rho_i_para*(1/k_per_rho_i_para**2+n/x*(1/k_per_rho_i_per**2-1/k_per_rho_i_para**2))*(1+x_pdf*y_pdf)
			
		      x_pdf=(x+cmplx(n,0))/k_para_rho_i_para
		      call pdf(x_pdf,y_pdf)
		      series_sum_1=series_sum_1-2*omega_pi_div_omega_ci/x*n*gamma_n*k_per_rho_i_para/k_para_rho_i_para*(1/k_per_rho_i_para**2-n/x*(1/k_per_rho_i_per**2-1/k_per_rho_i_para**2))*(1+x_pdf*y_pdf)
		      
		      call bessel_gn(k_per_rho_e_per**2/2,n, gamma_n)
		    
		      x_pdf=(x_e-cmplx(n,0))/k_para_rho_e_para
		      call pdf(x_pdf,y_pdf)
		      series_sum_1=series_sum_1+2*omega_pe_div_omega_ce/x_e*n*gamma_n*k_per_rho_e_para/k_para_rho_e_para*(1/k_per_rho_e_para**2+n/x_e*(1/k_per_rho_e_per**2-1/k_per_rho_e_para**2))*(1+x_pdf*y_pdf)
			
		      x_pdf=(x_e+cmplx(n,0))/k_para_rho_e_para
		      call pdf(x_pdf,y_pdf)
		      series_sum_1=series_sum_1-2*omega_pe_div_omega_ce/x_e*n*gamma_n*k_per_rho_e_para/k_para_rho_e_para*(1/k_per_rho_e_para**2-n/x_e*(1/k_per_rho_e_per**2-1/k_per_rho_e_para**2))*(1+x_pdf*y_pdf)
		      
		      if(n>=10 .and. abs(series_sum_1-series_sum_2)<1e-8*abs(series_sum_1+series_sum_2) )then
			  exit
		      else if(n==50) then
			  write(*,*) 'The sum of the first 50 series does not converge'
		      end if
		      series_sum_2=series_sum_1
		   end do
		      D(1,3)=c_div_v_para*k_para_rho_i_para*k_per_rho_i_para/x**2+series_sum_1
		    
		end if
	    case(4)
		D(2,1)=-D(1,2)
	    case(5)
		series_sum_1=(0,0)
		series_sum_2=(0,0)
		do n=0,10
		  
		  if(abs(k_para_rho_i_para)<1e-6) then
		    !call bessel_gn(k_per_rho_i_per**2/2,n, gamma_n)
		    !call bessel_gn(k_per_rho_i_per**2/2,n+1,gamma_n1)
		    !series_sum_1=series_sum_1+0.25*(gamma_n*(k_per_rho_i_per**4-2*n*k_per_rho_i_per**2+2*n**2)-k_per_rho_i_per**4*gamma_n1)*(1/k_per_rho_i_per**2/(cmplx(n,0)-x))
		    !if(n>0) then
		    !  series_sum_1=series_sum_1+0.25*(gamma_n*(k_per_rho_i_per**4-2*n*k_per_rho_i_per**2+2*n**2)-k_per_rho_i_per**4*gamma_n1)*(1/k_per_rho_i_per**2/(cmplx(-n,0)-x))
		    !end if
		    call bessel_gn(k_per_rho_i_per**2/2,n, gamma_n)
		    call bessel_gn(k_per_rho_i_per**2/2,n+1,gamma_n1)
		    series_sum_1=series_sum_1+omega_pi_div_omega_ci/x*((k_per_rho_i_per**2-2*n+2*n**2/k_per_rho_i_per**2)*gamma_n-k_per_rho_i_per**2*gamma_n1)/(cmplx(n,0)-x)
		    if(n>0) then
		      series_sum_1=series_sum_1+omega_pi_div_omega_ci/x*((k_per_rho_i_per**2-2*n+2*n**2/k_per_rho_i_per**2)*gamma_n-k_per_rho_i_per**2*gamma_n1)/(cmplx(-n,0)-x)
		    end if
		    call bessel_gn(k_per_rho_e_per**2/2,n, gamma_n)
		    call bessel_gn(k_per_rho_e_per**2/2,n+1,gamma_n1)
		    series_sum_1=series_sum_1+omega_pe_div_omega_ce/x_e*((k_per_rho_e_per**2-2*n+2*n**2/k_per_rho_e_per**2)*gamma_n-k_per_rho_e_per**2*gamma_n1)/(cmplx(n,0)-x_e)
		    if(n>0) then
		      series_sum_1=series_sum_1+omega_pe_div_omega_ce/x_e*((k_per_rho_e_per**2-2*n+2*n**2/k_per_rho_e_per**2)*gamma_n-k_per_rho_e_per**2*gamma_n1)/(cmplx(-n,0)-x_e)
		    end if
		    
		   
		  else
		    call bessel_gn(k_per_rho_i_per**2/2,n, gamma_n)
		    call bessel_gn(k_per_rho_i_per**2/2,n+1,gamma_n1)
		    x_pdf=(x-cmplx(n,0))/k_para_rho_i_para
		    call pdf(x_pdf,y_pdf)
		   
		    series_sum_1=series_sum_1+omega_pi_div_omega_ci/x*(gamma_n*(k_per_rho_i_per**4-2*n*k_per_rho_i_per**2+2*n**2)-k_per_rho_i_per**4*gamma_n1)*(1/k_per_rho_i_per**2/k_para_rho_i_para*y_pdf+1/x*(1/k_per_rho_i_para**2-1/k_per_rho_i_per**2)*(1+x_pdf*y_pdf))
		    if(n>0) then
		      x_pdf=(x+cmplx(n,0))/k_para_rho_i_para
		      call pdf(x_pdf,y_pdf)
		      series_sum_1=series_sum_1+omega_pi_div_omega_ci/x*(gamma_n*(k_per_rho_i_per**4-2*n*k_per_rho_i_per**2+2*n**2)-k_per_rho_i_per**4*gamma_n1)*(1/k_per_rho_i_per**2/k_para_rho_i_para*y_pdf+1/x*(1/k_per_rho_i_para**2-1/k_per_rho_i_per**2)*(1+x_pdf*y_pdf))
		   
		    end if
		    
		    call bessel_gn(k_per_rho_e_per**2/2,n, gamma_n)
		    call bessel_gn(k_per_rho_e_per**2/2,n+1,gamma_n1)
		    x_pdf=(x_e-cmplx(n,0))/k_para_rho_e_para
		    call pdf(x_pdf,y_pdf)
		   
		    series_sum_1=series_sum_1+omega_pe_div_omega_ce/x_e*(gamma_n*(k_per_rho_e_per**4-2*n*k_per_rho_e_per**2+2*n**2)-k_per_rho_e_per**4*gamma_n1)*(1/k_per_rho_e_per**2/k_para_rho_e_para*y_pdf+1/x_e*(1/k_per_rho_e_para**2-1/k_per_rho_e_per**2)*(1+x_pdf*y_pdf))
		    if(n>0) then
		      x_pdf=(x_e+cmplx(n,0))/k_para_rho_e_para
		      call pdf(x_pdf,y_pdf)
		      series_sum_1=series_sum_1+omega_pe_div_omega_ce/x_e*(gamma_n*(k_per_rho_e_per**4-2*n*k_per_rho_e_per**2+2*n**2)-k_per_rho_e_per**4*gamma_n1)*(1/k_per_rho_e_per**2/k_para_rho_e_para*y_pdf+1/x_e*(1/k_per_rho_e_para**2-1/k_per_rho_e_per**2)*(1+x_pdf*y_pdf))
		   
		    end if
		    
		   end if
		   
		   if(n>=10 .and. abs(series_sum_1-series_sum_2)<1e-8*abs(series_sum_1+series_sum_2) )then
		      exit
		   else if(n==50) then
		      write(*,*) 'The sum of the first 50 series does not converge'
		   end if
		   series_sum_2=series_sum_1
		end do
		if(abs(k_para_rho_i_para)<1e-6) then
		      D(2,2)=1-c_div_v_para*k_per_rho_i_para**2/x**2+series_sum_1
		else
		      D(2,2)=1-c_div_v_para*(k_para_rho_i_para**2+k_per_rho_i_para**2)/x**2+series_sum_1
		end if
        
	    case(6)
		series_sum_1=(0,0)
		series_sum_2=(0,0)
		if (abs(k_para_rho_i_para)<1e-6) then
		  D(2,3)=(0,0)
		
		else
		  
		  do n=0,10
		    call bessel_gn(k_per_rho_i_per**2/2,n, gamma_n)
		    call bessel_gn(k_per_rho_i_per**2/2,n+1,gamma_n1)
		    x_pdf=(x-cmplx(n,0))/k_para_rho_i_para
		    call pdf(x_pdf,y_pdf)
			      
		    series_sum_1=series_sum_1-cmplx(0,omega_pi_div_omega_ci,wp)/x*k_per_rho_i_para/k_para_rho_i_para*(k_per_rho_i_per**2*(gamma_n1-gamma_n)+2*n*gamma_n)*(1/k_per_rho_i_para**2+n/x*(1/k_per_rho_i_per**2-1/k_per_rho_i_para**2))*(1+x_pdf*y_pdf)
		    if(n>0) then
			  x_pdf=(x+cmplx(n,0))/k_para_rho_i_para
			  call pdf(x_pdf,y_pdf)
			  series_sum_1=series_sum_1-cmplx(0,omega_pi_div_omega_ci,wp)/x*k_per_rho_i_para/k_para_rho_i_para*(k_per_rho_i_per**2*(gamma_n1-gamma_n)+2*n*gamma_n)*(1/k_per_rho_i_para**2-n/x*(1/k_per_rho_i_per**2-1/k_per_rho_i_para**2))*(1+x_pdf*y_pdf)
		      
		    end if
		      
		      
		    call bessel_gn(k_per_rho_e_per**2/2,n, gamma_n)
		    call bessel_gn(k_per_rho_e_per**2/2,n+1,gamma_n1)
		    x_pdf=(x_e-cmplx(n,0))/k_para_rho_e_para
		    call pdf(x_pdf,y_pdf)
			      
		    series_sum_1=series_sum_1-cmplx(0,omega_pe_div_omega_ce,wp)/x_e*k_per_rho_e_para/k_para_rho_e_para*(k_per_rho_e_per**2*(gamma_n1-gamma_n)+2*n*gamma_n)*(1/k_per_rho_e_para**2+n/x_e*(1/k_per_rho_e_per**2-1/k_per_rho_e_para**2))*(1+x_pdf*y_pdf)
		    if(n>0) then
			  x_pdf=(x_e+cmplx(n,0))/k_para_rho_e_para
			  call pdf(x_pdf,y_pdf)
			  series_sum_1=series_sum_1-cmplx(0,omega_pe_div_omega_ce,wp)/x_e*k_per_rho_e_para/k_para_rho_e_para*(k_per_rho_e_per**2*(gamma_n1-gamma_n)+2*n*gamma_n)*(1/k_per_rho_e_para**2-n/x_e*(1/k_per_rho_e_per**2-1/k_per_rho_e_para**2))*(1+x_pdf*y_pdf)
		      
		    end if
		     if(n>=10 .and. abs(series_sum_1-series_sum_2)<1e-8*abs(series_sum_1+series_sum_2) )then
				exit
		     else if(n==50) then
				write(*,*) 'The sum of the first 50 series does not converge'
		     end if
		     series_sum_2=series_sum_1
		  end do
		  D(2,3)=series_sum_1
		end if
        case(7)
            D(3,1)=D(1,3)
        !         
		!series_sum_1=(0,0)
		!series_sum_2=(0,0)
		!if(abs(k_para_rho_i_para)<1e-6) then
		!	D(3,1)=(0,0)
		!else
	    ! 
		!	do n=1,10
		!	call bessel_gn(k_per_rho_i_per**2/2,n,gamma_n)
		!	x_pdf=(x-cmplx(n,0))/k_para_rho_i_para
		!	call pdf(x_pdf,y_pdf)
		!	series_sum_1=series_sum_1+2*omega_pi_div_omega_ci*k_per_rho_i_para/k_para_rho_i_para/x*gamma_n*n*(n/x/k_per_rho_i_per**2+((1,0)-n/x)/k_per_rho_i_para**2)*(1+x_pdf*y_pdf)
	    !   
		!	x_pdf=(x+cmplx(n,0))/k_para_rho_i_para
		!	call pdf(x_pdf,y_pdf)
		!	series_sum_1=series_sum_1+2*omega_pi_div_omega_ci*k_per_rho_i_para/k_para_rho_i_para/x*gamma_n*(-n)*(-n/x/k_per_rho_i_per**2+((1,0)+n/x)/k_per_rho_i_para**2)*(1+x_pdf*y_pdf)
	    !   
		!	call bessel_gn(k_per_rho_e_per**2/2,n,gamma_n)
		!	x_pdf=(x_e-cmplx(n,0))/k_para_rho_e_para
		!	call pdf(x_pdf,y_pdf)
		!	series_sum_1=series_sum_1+2*omega_pe_div_omega_ce*k_per_rho_e_para/k_para_rho_e_para/x_e*gamma_n*n*(n/x_e/k_per_rho_e_per**2+((1,0)-n/x_e)/k_per_rho_e_para**2)*(1+x_pdf*y_pdf)
	    !   
		!	x_pdf=(x_e+cmplx(n,0))/k_para_rho_e_para
		!	call pdf(x_pdf,y_pdf)
		!	series_sum_1=series_sum_1+2*omega_pe_div_omega_ce*k_per_rho_e_para/k_para_rho_e_para/x_e*gamma_n*(-n)*(-n/x_e/k_per_rho_e_per**2+((1,0)+n/x_e)/k_per_rho_e_para**2)*(1+x_pdf*y_pdf)
	    !   
		!	if(n>=10 .and. abs(series_sum_1-series_sum_2)<1e-8*abs(series_sum_1+series_sum_2) )then
		!		exit
		!		else if(n==50) then
		!		write(*,*) 'The sum of the first 50 series does not converge'
		!		end if
		!		series_sum_2=series_sum_1
		!	end do
		!	D(3,1)=c_div_v_para*k_para_rho_i_para*k_per_rho_i_para/x**2+series_sum_1
		!end if
	    
        case(8)
           D(3,2)=-D(2,3)
			!series_sum_1=(0,0)
			!series_sum_2=(0,0)
			!if(abs(k_para_rho_i_para)<1e-6) then
			!  D(3,2)=(0,0)
			!else
	  !
			!  do n=0,10
			!	call bessel_gn(k_per_rho_i_per**2/2,n,gamma_n)
			!	call bessel_gn(k_per_rho_i_per**2/2,n+1,gamma_n1)
			!	x_pdf=(x-cmplx(n,0))/k_para_rho_i_para
			!	call pdf(x_pdf,y_pdf)
			!	series_sum_1=series_sum_1+cmplx(0,omega_pi_div_omega_ci*k_per_rho_i_para/k_para_rho_i_para,wp)/x*(k_per_rho_i_per**2*(gamma_n1-gamma_n)+2*n*gamma_n)*(n/x/k_per_rho_i_per**2+((1,0)-n/x)/k_per_rho_i_para**2)*(1+x_pdf*y_pdf)
	  !  
			!	if (n>0) then
			!	  x_pdf=(x+cmplx(n,0))/k_para_rho_i_para
			!	  call pdf(x_pdf,y_pdf)
			!	  series_sum_1=series_sum_1+cmplx(0,omega_pi_div_omega_ci*k_per_rho_i_para/k_para_rho_i_para,wp)/x*(k_per_rho_i_per**2*(gamma_n1-gamma_n)+2*n*gamma_n)*(-n/x/k_per_rho_i_per**2+((1,0)+n/x)/k_per_rho_i_para**2)*(1+x_pdf*y_pdf)
			!	end if
	  !  
			!	call bessel_gn(k_per_rho_e_per**2/2,n,gamma_n)
			!	call bessel_gn(k_per_rho_e_per**2/2,n+1,gamma_n1)
			!	x_pdf=(x_e-cmplx(n,0))/k_para_rho_e_para
			!	call pdf(x_pdf,y_pdf)
			!	series_sum_1=series_sum_1+cmplx(0,omega_pe_div_omega_ce*k_per_rho_e_para/k_para_rho_e_para,wp)/x_e*(k_per_rho_e_per**2*(gamma_n1-gamma_n)+2*n*gamma_n)*(n/x_e/k_per_rho_e_per**2+((1,0)-n/x_e)/k_per_rho_e_para**2)*(1+x_pdf*y_pdf)
			!	if (n>0) then
			!	  x_pdf=(x_e+cmplx(n,0))/k_para_rho_e_para
			!	  call pdf(x_pdf,y_pdf)
			!	  series_sum_1=series_sum_1+cmplx(0,omega_pe_div_omega_ce*k_per_rho_e_para/k_para_rho_e_para,wp)/x_e*(k_per_rho_e_per**2*(gamma_n1-gamma_n)+2*n*gamma_n)*(-n/x_e/k_per_rho_e_per**2+((1,0)+n/x_e)/k_per_rho_e_para**2)*(1+x_pdf*y_pdf)
			!	end if
			!	if(n>=10 .and. abs(series_sum_1-series_sum_2)<1e-8*abs(series_sum_1+series_sum_2) )then
			!	  exit
			!	else if(n==50) then
			!	  write(*,*) 'The sum of the first 50 series does not converge'
			!	end if
			!	series_sum_2=series_sum_1
			!  end do
			!  D(3,2)=series_sum_1
			!end if
	  !      
	     
        case(9)
		series_sum_1=(0,0)
		series_sum_2=(0,0)
		do n=0,10
		  
		  if (abs(k_para_rho_i_para)<1e-6) then
		    call bessel_gn(k_per_rho_i_per**2/2,n,gamma_n)
		    series_sum_1=series_sum_1+omega_pi_div_omega_ci/x*gamma_n/(cmplx(n,0)-x)
		    
		    if (n>0) then
		      series_sum_1=series_sum_1+omega_pi_div_omega_ci/x*gamma_n/(cmplx(-n,0)-x)
		    end if
		    
		    call bessel_gn(k_per_rho_e_per**2/2,n,gamma_n)
		    series_sum_1=series_sum_1+omega_pe_div_omega_ce/x_e*gamma_n/(cmplx(n,0)-x_e)
		    
		    if (n>0) then
		      series_sum_1=series_sum_1+omega_pe_div_omega_ce/x_e*gamma_n/(cmplx(-n,0)-x_e)
		    end if
		  else
		    call bessel_gn(k_per_rho_i_per**2/2,n,gamma_n)
		    x_pdf=(x-cmplx(n,0))/k_para_rho_i_para
		    call pdf(x_pdf,y_pdf)
		    series_sum_1=series_sum_1+2*omega_pi_div_omega_ci*k_per_rho_i_para**2/x/k_para_rho_i_para*x_pdf*gamma_n*(1/k_per_rho_i_para**2+n/x*(1/k_per_rho_i_per**2-1/k_per_rho_i_para**2))*(1+x_pdf*y_pdf)
		      
		    if (n>0) then
		      x_pdf=(x+cmplx(n,0))/k_para_rho_i_para
		      call pdf(x_pdf,y_pdf)
		      series_sum_1=series_sum_1+2*omega_pi_div_omega_ci*k_per_rho_i_para**2/x/k_para_rho_i_para*x_pdf*gamma_n*(1/k_per_rho_i_para**2-n/x*(1/k_per_rho_i_per**2-1/k_per_rho_i_para**2))*(1+x_pdf*y_pdf)
		    end if
		    
		    
		    call bessel_gn(k_per_rho_e_per**2/2,n,gamma_n)
		    x_pdf=(x_e-cmplx(n,0))/k_para_rho_e_para
		    call pdf(x_pdf,y_pdf)
		    series_sum_1=series_sum_1+2*omega_pe_div_omega_ce*k_per_rho_e_para**2/x_e/k_para_rho_e_para*x_pdf*gamma_n*(1/k_per_rho_e_para**2+n/x_e*(1/k_per_rho_e_per**2-1/k_per_rho_e_para**2))*(1+x_pdf*y_pdf)
		      
		    if (n>0) then
		      x_pdf=(x_e+cmplx(n,0))/k_para_rho_e_para
		      call pdf(x_pdf,y_pdf)
		      series_sum_1=series_sum_1+2*omega_pe_div_omega_ce*k_per_rho_e_para**2/x_e/k_para_rho_e_para*x_pdf*gamma_n*(1/k_per_rho_e_para**2-n/x_e*(1/k_per_rho_e_per**2-1/k_per_rho_e_para**2))*(1+x_pdf*y_pdf)
		    end if
		  end if
		  
		  if(n>=10 .and. abs(series_sum_1-series_sum_2)<1e-8*abs(series_sum_1+series_sum_2) )then
		      exit
		  else if(n==50) then
		      write(*,*) 'The sum of the first 50 series does not converge'
		  end if
		  series_sum_2=series_sum_1
	      end do
	      D(3,3)=1-k_per_rho_i_para**2*c_div_v_para/x**2+series_sum_1
	end select
	    
	    
    end do
    
    end subroutine dispersion_function_matrix
    
    
!-----------------------------------------------------------------------------!
!     dispersion_function_matrix_variable_k_para:dispersion matrix for k_para with known k_perp and omega
!-----------------------------------------------------------------------------!
    subroutine dispersion_function_matrix_variable_k_para(x,D)
    implicit none
    complex(wp),intent(in)::x
    complex(wp),intent(out)::D(:,:)
    complex(wp)::series_sum_1,series_sum_2,x_pdf,y_pdf
    complex(wp)::xe_pdf_1,xe_pdf_2,ye_pdf_1,ye_pdf_2
    integer::n,k
    real(wp)::gamma_n,gamma_n1
    complex(wp)::x_e
	
	x_e=-x/(1836)**(0.5)
	do k=1,9
	    select case(k)
	    case(1)
		series_sum_1=(0,0)
		series_sum_2=(0,0)
		do n=1,10
		    call bessel_gn(k_per_rho_i_per**2/2,n, gamma_n)
		    call bessel_gn(k_per_rho_e_per**2/2,n, gamma_n1)
		    x_pdf=(omega_div_omega_ci-cmplx(n,0))/x
		    call pdf(x_pdf,y_pdf)
		      
		    series_sum_1=series_sum_1+2*omega_pi_div_omega_ci/omega_div_omega_ci*gamma_n*(n**2/k_per_rho_i_per**2/x*y_pdf)
		      
		    x_pdf=(omega_div_omega_ci+cmplx(n,0))/x
		    call pdf(x_pdf,y_pdf)
		    
		    series_sum_1=series_sum_1+2*omega_pi_div_omega_ci/omega_div_omega_ci*gamma_n*(n**2/k_per_rho_i_per**2/x*y_pdf)
		      
		      
		    x_pdf=(omega_div_omega_ce-cmplx(n,0))/x_e
		    call pdf(x_pdf,y_pdf)
		      
		    series_sum_1=series_sum_1+2*omega_pe_div_omega_ce/omega_div_omega_ce*gamma_n1*(n**2/k_per_rho_e_per**2/x_e*y_pdf)
		    
            x_pdf=(omega_div_omega_ce+cmplx(n,0))/x_e
		    call pdf(x_pdf,y_pdf)
		    
		    series_sum_1=series_sum_1+2*omega_pe_div_omega_ce/omega_div_omega_ce*gamma_n1*(n**2/k_per_rho_e_per**2/x_e*y_pdf)
		      
		    if(n>=10 .and. abs(series_sum_1-series_sum_2)<1e-8*abs(series_sum_1+series_sum_2) )then
			
			exit
		    else if(n==50) then
			write(*,*) 'The sum of the first 50 series does not converge'
		    end if
		    series_sum_2=series_sum_1
		end do
        D(1,1)=(1,0)-c_div_v_para*x**2/omega_div_omega_ci**2+series_sum_1

	  case(2)
		series_sum_1=(0,0)
		series_sum_2=(0,0)
		do n=1,10
			call bessel_gn(k_per_rho_i_per**2/2,n, gamma_n)
			call bessel_gn(k_per_rho_i_per**2/2,n+1,gamma_n1)
			x_pdf=(omega_div_omega_ci-cmplx(n,0))/x
			call pdf(x_pdf,y_pdf)
		      
			series_sum_1=series_sum_1+cmplx(0,omega_pi_div_omega_ci,wp)/omega_div_omega_ci*(k_per_rho_i_per**2*(gamma_n1-gamma_n)+2*n*gamma_n)*(n/k_per_rho_i_per**2/x*y_pdf)
		    
		    
			x_pdf=(omega_div_omega_ci+cmplx(n,0))/x
			call pdf(x_pdf,y_pdf)
		   
			series_sum_1=series_sum_1-cmplx(0,omega_pi_div_omega_ci,wp)/omega_div_omega_ci*(k_per_rho_i_per**2*(gamma_n1-gamma_n)+2*n*gamma_n)*(n/k_per_rho_i_per**2/x*y_pdf)
		    
			call bessel_gn(k_per_rho_e_per**2/2,n, gamma_n)
			call bessel_gn(k_per_rho_e_per**2/2,n+1,gamma_n1)
			x_pdf=(omega_div_omega_ce-cmplx(n,0))/x_e
			call pdf(x_pdf,y_pdf)
		      
			series_sum_1=series_sum_1+cmplx(0,omega_pe_div_omega_ce,wp)/omega_div_omega_ce*(k_per_rho_e_per**2*(gamma_n1-gamma_n)+2*n*gamma_n)*(n/k_per_rho_e_per**2/x_e*y_pdf)
		    
		    
			x_pdf=(omega_div_omega_ce+cmplx(n,0))/x_e
			call pdf(x_pdf,y_pdf)
		   
			series_sum_1=series_sum_1-cmplx(0,omega_pe_div_omega_ce,wp)/omega_div_omega_ce*(k_per_rho_e_per**2*(gamma_n1-gamma_n)+2*n*gamma_n)*(n/k_per_rho_e_per**2/x_e*y_pdf)
		    
		    if(n>=10 .and. abs(series_sum_1-series_sum_2)<1e-8*abs(series_sum_1+series_sum_2) )then
				exit
            else if(n==50) then
				write(*,*) 'The sum of the first 50 series does not converge'
		    end if
		    series_sum_2=series_sum_1
		end do
        D(1,2)=series_sum_1
 
	    case(3)
		series_sum_1=(0,0)
		series_sum_2=(0,0)
		   do n=1,10
		      call bessel_gn(k_per_rho_i_per**2/2,n, gamma_n)
		    
		      x_pdf=(omega_div_omega_ci-cmplx(n,0))/x
		      call pdf(x_pdf,y_pdf)
		      series_sum_1=series_sum_1+2*omega_pi_div_omega_ci/omega_div_omega_ci*n*gamma_n*k_per_rho_i_para/x*(1/k_per_rho_i_para**2)*(1+x_pdf*y_pdf)
			
		      x_pdf=(omega_div_omega_ci+cmplx(n,0))/x
		      call pdf(x_pdf,y_pdf)
		      series_sum_1=series_sum_1-2*omega_pi_div_omega_ci/omega_div_omega_ci*n*gamma_n*k_per_rho_i_para/x*(1/k_per_rho_i_para**2)*(1+x_pdf*y_pdf)
		      
		      call bessel_gn(k_per_rho_e_per**2/2,n, gamma_n)
		    
		      x_pdf=(omega_div_omega_ce-cmplx(n,0))/x_e
		      call pdf(x_pdf,y_pdf)
		      series_sum_1=series_sum_1+2*omega_pe_div_omega_ce/omega_div_omega_ce*n*gamma_n*k_per_rho_e_para/x_e*(1/k_per_rho_e_para**2)*(1+x_pdf*y_pdf)
			
		      x_pdf=(omega_div_omega_ce+cmplx(n,0))/x_e
		      call pdf(x_pdf,y_pdf)
		      series_sum_1=series_sum_1-2*omega_pe_div_omega_ce/omega_div_omega_ce*n*gamma_n*k_per_rho_e_para/x_e*(1/k_per_rho_e_para**2)*(1+x_pdf*y_pdf)
		      
		      if(n>=10 .and. abs(series_sum_1-series_sum_2)<1e-8*abs(series_sum_1+series_sum_2) )then
				  exit
				  else if(n==50) then
				  write(*,*) 'The sum of the first 50 series does not converge'
		      end if
		      series_sum_2=series_sum_1
		   end do
           D(1,3)=c_div_v_para*x*k_per_rho_i_para/omega_div_omega_ci**2+series_sum_1
		    
	    case(4)
		D(2,1)=-D(1,2)
	    case(5)
		series_sum_1=(0,0)
		series_sum_2=(0,0)
		do n=0,10
			call bessel_gn(k_per_rho_i_per**2/2,n, gamma_n)
			call bessel_gn(k_per_rho_i_per**2/2,n+1,gamma_n1)
			x_pdf=(omega_div_omega_ci-cmplx(n,0))/x
			call pdf(x_pdf,y_pdf)
			series_sum_1=series_sum_1+omega_pi_div_omega_ci/omega_div_omega_ci*(gamma_n*(k_per_rho_i_per**4-2*n*k_per_rho_i_per**2+2*n**2)-k_per_rho_i_per**4*gamma_n1)*(1/k_per_rho_i_per**2/x*y_pdf)
			if(n>0) then
				x_pdf=(omega_div_omega_ci+cmplx(n,0))/x
				call pdf(x_pdf,y_pdf)
				series_sum_1=series_sum_1+omega_pi_div_omega_ci/omega_div_omega_ci*(gamma_n*(k_per_rho_i_per**4-2*n*k_per_rho_i_per**2+2*n**2)-k_per_rho_i_per**4*gamma_n1)*(1/k_per_rho_i_per**2/x*y_pdf)
		   
			end if
		    
			call bessel_gn(k_per_rho_e_per**2/2,n, gamma_n)
			call bessel_gn(k_per_rho_e_per**2/2,n+1,gamma_n1)
			x_pdf=(omega_div_omega_ce-cmplx(n,0))/x_e
			call pdf(x_pdf,y_pdf)
		   
			series_sum_1=series_sum_1+omega_pe_div_omega_ce/omega_div_omega_ce*(gamma_n*(k_per_rho_e_per**4-2*n*k_per_rho_e_per**2+2*n**2)-k_per_rho_e_per**4*gamma_n1)*(1/k_per_rho_e_per**2/x_e*y_pdf)
			if(n>0) then
				x_pdf=(omega_div_omega_ce+cmplx(n,0))/x_e
				call pdf(x_pdf,y_pdf)
				series_sum_1=series_sum_1+omega_pe_div_omega_ce/omega_div_omega_ce*(gamma_n*(k_per_rho_e_per**4-2*n*k_per_rho_e_per**2+2*n**2)-k_per_rho_e_per**4*gamma_n1)*(1/k_per_rho_e_per**2/x_e*y_pdf)
			end if
		   
			if(n>=10 .and. abs(series_sum_1-series_sum_2)<1e-8*abs(series_sum_1+series_sum_2) )then
				exit
			else if(n==50) then
				write(*,*) 'The sum of the first 50 series does not converge'
			end if
			series_sum_2=series_sum_1
		end do
        D(2,2)=1-c_div_v_para*(x**2+k_per_rho_i_para**2)/omega_div_omega_ci**2+series_sum_1

        
	    case(6)
		series_sum_1=(0,0)
		series_sum_2=(0,0)
		  do n=0,10
		    call bessel_gn(k_per_rho_i_per**2/2,n, gamma_n)
		    call bessel_gn(k_per_rho_i_per**2/2,n+1,gamma_n1)
		    x_pdf=(omega_div_omega_ci-cmplx(n,0))/x
		    call pdf(x_pdf,y_pdf)
			      
		    series_sum_1=series_sum_1-cmplx(0,omega_pi_div_omega_ci,wp)/omega_div_omega_ci*k_per_rho_i_para/x*(k_per_rho_i_per**2*(gamma_n1-gamma_n)+2*n*gamma_n)*(1/k_per_rho_i_para**2)*(1+x_pdf*y_pdf)
		    if(n>0) then
			  x_pdf=(omega_div_omega_ci+cmplx(n,0))/x
			  call pdf(x_pdf,y_pdf)
			  series_sum_1=series_sum_1-cmplx(0,omega_pi_div_omega_ci,wp)/omega_div_omega_ci*k_per_rho_i_para/x*(k_per_rho_i_per**2*(gamma_n1-gamma_n)+2*n*gamma_n)*(1/k_per_rho_i_para**2)*(1+x_pdf*y_pdf)
		    end if
		      
		    call bessel_gn(k_per_rho_e_per**2/2,n, gamma_n)
		    call bessel_gn(k_per_rho_e_per**2/2,n+1,gamma_n1)
		    x_pdf=(omega_div_omega_ce-cmplx(n,0))/x_e
		    call pdf(x_pdf,y_pdf)
			      
		    series_sum_1=series_sum_1-cmplx(0,omega_pe_div_omega_ce,wp)/omega_div_omega_ce*k_per_rho_e_para/x_e*(k_per_rho_e_per**2*(gamma_n1-gamma_n)+2*n*gamma_n)*(1/k_per_rho_e_para**2)*(1+x_pdf*y_pdf)
		    if(n>0) then
			  x_pdf=(omega_div_omega_ce+cmplx(n,0))/x_e
			  call pdf(x_pdf,y_pdf)
			  series_sum_1=series_sum_1-cmplx(0,omega_pe_div_omega_ce,wp)/omega_div_omega_ce*k_per_rho_e_para/x_e*(k_per_rho_e_per**2*(gamma_n1-gamma_n)+2*n*gamma_n)*(1/k_per_rho_e_para**2)*(1+x_pdf*y_pdf)
		      
		    end if
		     if(n>=10 .and. abs(series_sum_1-series_sum_2)<1e-8*abs(series_sum_1+series_sum_2) )then
				exit
		     else if(n==50) then
				write(*,*) 'The sum of the first 50 series does not converge'
		     end if
		     series_sum_2=series_sum_1
		  end do
		  D(2,3)=series_sum_1

        case(7)
			D(3,1)=D(1,3)
        
        case(8)
			D(3,2)=-D(2,3)
	     
        case(9)
		series_sum_1=(0,0)
		series_sum_2=(0,0)
		do n=0,10
			call bessel_gn(k_per_rho_i_per**2/2,n,gamma_n)
			x_pdf=(omega_div_omega_ci-cmplx(n,0))/x
			call pdf(x_pdf,y_pdf)
			series_sum_1=series_sum_1+2*omega_pi_div_omega_ci*k_per_rho_i_para**2/omega_div_omega_ci/x*x_pdf*gamma_n*(1/k_per_rho_i_para**2)*(1+x_pdf*y_pdf)
		      
			if (n>0) then
				x_pdf=(omega_div_omega_ci+cmplx(n,0))/x
				call pdf(x_pdf,y_pdf)
				series_sum_1=series_sum_1+2*omega_pi_div_omega_ci*k_per_rho_i_para**2/omega_div_omega_ci/x*x_pdf*gamma_n*(1/k_per_rho_i_para**2)*(1+x_pdf*y_pdf)
			end if
		    
			call bessel_gn(k_per_rho_e_per**2/2,n,gamma_n)
			x_pdf=(omega_div_omega_ce-cmplx(n,0))/x_e
			call pdf(x_pdf,y_pdf)
			series_sum_1=series_sum_1+2*omega_pe_div_omega_ce*k_per_rho_e_para**2/omega_div_omega_ce/x_e*x_pdf*gamma_n*(1/k_per_rho_e_para**2)*(1+x_pdf*y_pdf)
		      
			if (n>0) then
				x_pdf=(omega_div_omega_ce+cmplx(n,0))/x_e
				call pdf(x_pdf,y_pdf)
				series_sum_1=series_sum_1+2*omega_pe_div_omega_ce*k_per_rho_e_para**2/omega_div_omega_ce/x_e*x_pdf*gamma_n*(1/k_per_rho_e_para**2)*(1+x_pdf*y_pdf)
			end if

			if(n>=10 .and. abs(series_sum_1-series_sum_2)<1e-8*abs(series_sum_1+series_sum_2) )then
				exit
			else if(n==50) then
				write(*,*) 'The sum of the first 50 series does not converge'
			end if
			series_sum_2=series_sum_1
	      end do
	      D(3,3)=1-k_per_rho_i_para**2*c_div_v_para/omega_div_omega_ci**2+series_sum_1
	end select
	    
	    
    end do
    
    end subroutine dispersion_function_matrix_variable_k_para    

!-----------------------------------------------------------------------------!
!     dispersion_function_matrix_k_per:dipsersion matrix for k_perp with given omega and k_para
!-----------------------------------------------------------------------------!
    subroutine dispersion_function_matrix_variable_k_per(x,D)
    implicit none
    complex(wp),intent(in)::x
    complex(wp),intent(out)::D(:,:)
    complex(wp)::series_sum_1,series_sum_2,x_pdf,y_pdf
    complex(wp)::xe_pdf_1,xe_pdf_2,ye_pdf_1,ye_pdf_2
    integer::n,k
    complex(wp)::gamma_n,gamma_n1
    complex(wp)::x_e
	
	x_e=-x/(1836)**(0.5)
	do k=1,9
	    select case(k)
	    case(1)
		series_sum_1=(0,0)
		series_sum_2=(0,0)
		do n=1,10
            gamma_n=bessel_gn_complex(n,x**2/2)
		    gamma_n1=bessel_gn_complex(n,x_e**2/2)
		    if(abs(k_para_rho_i_para)<1e-6) then
		      series_sum_1=series_sum_1+4/x**2*omega_pi_div_omega_ci*gamma_n*(n**2/(cmplx(n**2,0)-omega_div_omega_ci**2))
		      series_sum_1=series_sum_1+4/x_e**2*omega_pe_div_omega_ce*gamma_n1*(n**2/(cmplx(n**2,0)-omega_div_omega_ce**2))
		    else
		      x_pdf=(omega_div_omega_ci-cmplx(n,0))/k_para_rho_i_para
		      call pdf(x_pdf,y_pdf)
		      
		      series_sum_1=series_sum_1+2*omega_pi_div_omega_ci/omega_div_omega_ci*gamma_n*(n**2/x**2/k_para_rho_i_para*y_pdf)
		      
		      x_pdf=(omega_div_omega_ci+cmplx(n,0))/k_para_rho_i_para
		      call pdf(x_pdf,y_pdf)
		    
		      series_sum_1=series_sum_1+2*omega_pi_div_omega_ci/omega_div_omega_ci*gamma_n*(n**2/x**2/k_para_rho_i_para*y_pdf)
		      
		      
		      x_pdf=(omega_div_omega_ce-cmplx(n,0))/k_para_rho_e_para
		      call pdf(x_pdf,y_pdf)
		      
		      series_sum_1=series_sum_1+2*omega_pe_div_omega_ce/omega_div_omega_ce*gamma_n1*(n**2/x_e**2/k_para_rho_e_para*y_pdf)
		      
		      x_pdf=(omega_div_omega_ce+cmplx(n,0))/k_para_rho_e_para
		      call pdf(x_pdf,y_pdf)
		    
		      series_sum_1=series_sum_1+2*omega_pe_div_omega_ce/omega_div_omega_ce*gamma_n1*(n**2/x_e**2/k_para_rho_e_para*y_pdf)
		      
		    end if
		    
		    if(n>=10 .and. abs(series_sum_1-series_sum_2)<1e-8*abs(series_sum_1+series_sum_2) )then
				exit
            else if(n==50) then
				write(*,*) 'The sum of the first 50 series does not converge'
		    end if
		    series_sum_2=series_sum_1
		end do
		
		if (abs(k_para_rho_i_para)<1e-6) then
		    D(1,1)=(1,0)+series_sum_1
		else
		    D(1,1)=(1,0)-c_div_v_para*k_para_rho_i_para**2/omega_div_omega_ci**2+series_sum_1
		end if
        
	  case(2)
		series_sum_1=(0,0)
		series_sum_2=(0,0)
		do n=1,10
		    if(abs(k_para_rho_i_para)<1e-6) then
              gamma_n=bessel_gn_complex(n,x**2/2)
              gamma_n1=bessel_gn_complex(n+1,x**2/2)
              series_sum_1=series_sum_1+(0,1)*omega_pi_div_omega_ci/omega_div_omega_ci*(2*n**2)/(cmplx(n**2,0)-omega_div_omega_ci**2)*(gamma_n1-(1-2*n/x**2)*gamma_n)
              gamma_n=bessel_gn_complex(n,x_e**2/2)
              gamma_n1=bessel_gn_complex(n+1,x_e**2/2)
		      series_sum_1=series_sum_1+(0,1)*omega_pe_div_omega_ce/omega_div_omega_ce*(2*n**2)/(cmplx(n**2,0)-omega_div_omega_ce**2)*(gamma_n1-(1-2*n/x_e**2)*gamma_n)
		    else
              gamma_n=bessel_gn_complex(n,x**2/2)
              gamma_n1=bessel_gn_complex(n+1,x**2/2)
		      x_pdf=(omega_div_omega_ci-cmplx(n,0))/k_para_rho_i_para
		      call pdf(x_pdf,y_pdf)
		      
		      series_sum_1=series_sum_1+cmplx(0,omega_pi_div_omega_ci,wp)/omega_div_omega_ci*(x**2*(gamma_n1-gamma_n)+2*n*gamma_n)*(n/x**2/k_para_rho_i_para*y_pdf)
		    
		      x_pdf=(omega_div_omega_ci+cmplx(n,0))/k_para_rho_i_para
		      call pdf(x_pdf,y_pdf)
		   
		      series_sum_1=series_sum_1-cmplx(0,omega_pi_div_omega_ci,wp)/omega_div_omega_ci*(x**2*(gamma_n1-gamma_n)+2*n*gamma_n)*(n/x**2/k_para_rho_i_para*y_pdf)
		    
              gamma_n=bessel_gn_complex(n,x_e**2/2)
              gamma_n1=bessel_gn_complex(n+1,x_e**2/2)
		      x_pdf=(omega_div_omega_ce-cmplx(n,0))/k_para_rho_e_para
		      call pdf(x_pdf,y_pdf)
		      
		      series_sum_1=series_sum_1+cmplx(0,omega_pe_div_omega_ce,wp)/omega_div_omega_ce*(x_e**2*(gamma_n1-gamma_n)+2*n*gamma_n)*(n/x_e**2/k_para_rho_e_para*y_pdf)
		    
		    
		      x_pdf=(omega_div_omega_ce+cmplx(n,0))/k_para_rho_e_para
		      call pdf(x_pdf,y_pdf)
		   
		      series_sum_1=series_sum_1-cmplx(0,omega_pe_div_omega_ce,wp)/omega_div_omega_ce*(x_e**2*(gamma_n1-gamma_n)+2*n*gamma_n)*(n/x_e**2/k_para_rho_e_para*y_pdf)

		    end if
		    
		    if(n>=10 .and. abs(series_sum_1-series_sum_2)<1e-8*abs(series_sum_1+series_sum_2) )then
				exit
            else if(n==50) then
				write(*,*) 'The sum of the first 50 series does not converge'
		    end if
		    series_sum_2=series_sum_1
		end do
			if (abs(k_para_rho_i_para)<1e-6) then
				D(1,2)=series_sum_1
			else
				D(1,2)=series_sum_1
			end if
          
	    case(3)
		series_sum_1=(0,0)
		series_sum_2=(0,0)
		if(abs(k_para_rho_i_para)<1e-6) then
		  D(1,3)=(0,0)
		else
		   do n=1,10
              gamma_n=bessel_gn_complex(n,x**2/2)
		      x_pdf=(omega_div_omega_ci-cmplx(n,0))/k_para_rho_i_para
		      call pdf(x_pdf,y_pdf)
		      series_sum_1=series_sum_1+2*omega_pi_div_omega_ci/omega_div_omega_ci*n*gamma_n*x/k_para_rho_i_para*(1/x**2)*(1+x_pdf*y_pdf)
			
		      x_pdf=(omega_div_omega_ci+cmplx(n,0))/k_para_rho_i_para
		      call pdf(x_pdf,y_pdf)
		      series_sum_1=series_sum_1-2*omega_pi_div_omega_ci/omega_div_omega_ci*n*gamma_n*x/k_para_rho_i_para*(1/x**2)*(1+x_pdf*y_pdf)

              gamma_n=bessel_gn_complex(n,x_e**2/2)
		      x_pdf=(omega_div_omega_ce-cmplx(n,0))/k_para_rho_e_para
		      call pdf(x_pdf,y_pdf)
		      series_sum_1=series_sum_1+2*omega_pe_div_omega_ce/omega_div_omega_ce*n*gamma_n*x_e/k_para_rho_e_para*(1/x_e**2)*(1+x_pdf*y_pdf)
			
		      x_pdf=(omega_div_omega_ce+cmplx(n,0))/k_para_rho_e_para
		      call pdf(x_pdf,y_pdf)
		      series_sum_1=series_sum_1-2*omega_pe_div_omega_ce/omega_div_omega_ce*n*gamma_n*x_e/k_para_rho_e_para*(1/x_e**2)*(1+x_pdf*y_pdf)
		      
		      if(n>=10 .and. abs(series_sum_1-series_sum_2)<1e-8*abs(series_sum_1+series_sum_2) )then
				  exit
              else if(n==50) then
				  write(*,*) 'The sum of the first 50 series does not converge'
		      end if
		      series_sum_2=series_sum_1
		   end do
		      D(1,3)=c_div_v_para*k_para_rho_i_para*x/omega_div_omega_ci**2+series_sum_1
		    
		end if
	    case(4)
			D(2,1)=-D(1,2)
        case(5)
            
		series_sum_1=(0,0)
		series_sum_2=(0,0)
		do n=0,10
		  
		  if(abs(k_para_rho_i_para)<1e-6) then
            gamma_n=bessel_gn_complex(n,x**2/2)
			gamma_n1=bessel_gn_complex(n+1,x**2/2)
		    series_sum_1=series_sum_1+omega_pi_div_omega_ci/omega_div_omega_ci*((x**2-2*n+2*n**2/x**2)*gamma_n-x**2*gamma_n1)/(cmplx(n,0)-omega_div_omega_ci)
		    if(n>0) then
		      series_sum_1=series_sum_1+omega_pi_div_omega_ci/omega_div_omega_ci*((x**2-2*n+2*n**2/x**2)*gamma_n-x**2*gamma_n1)/(cmplx(-n,0)-omega_div_omega_ci)
            end if
            gamma_n=bessel_gn_complex(n,x_e**2/2)
			gamma_n1=bessel_gn_complex(n+1,x_e**2/2)
		    series_sum_1=series_sum_1+omega_pe_div_omega_ce/omega_div_omega_ce*((x_e**2-2*n+2*n**2/x_e**2)*gamma_n-x_e**2*gamma_n1)/(cmplx(n,0)-omega_div_omega_ce)
		    if(n>0) then
		      series_sum_1=series_sum_1+omega_pe_div_omega_ce/omega_div_omega_ce*((x_e**2-2*n+2*n**2/x_e**2)*gamma_n-x_e**2*gamma_n1)/(cmplx(-n,0)-omega_div_omega_ce)
		    end if
		    
		  else
            gamma_n=bessel_gn_complex(n,x**2/2)
			gamma_n1=bessel_gn_complex(n+1,x**2/2)
		    x_pdf=(omega_div_omega_ci-cmplx(n,0))/k_para_rho_i_para
		    call pdf(x_pdf,y_pdf)
            
		    series_sum_1=series_sum_1+omega_pi_div_omega_ci/omega_div_omega_ci*(gamma_n*(x**4-2*n*x**2+2*n**2)-x**4*gamma_n1)*(1/x**2/k_para_rho_i_para*y_pdf)
		    if(n>0) then
		      x_pdf=(omega_div_omega_ci+cmplx(n,0))/k_para_rho_i_para
		      call pdf(x_pdf,y_pdf)
		      series_sum_1=series_sum_1+omega_pi_div_omega_ci/omega_div_omega_ci*(gamma_n*(x**4-2*n*x**2+2*n**2)-x**4*gamma_n1)*(1/x**2/k_para_rho_i_para*y_pdf)
		   
		    end if
            gamma_n=bessel_gn_complex(n,x_e**2/2)
			gamma_n1=bessel_gn_complex(n+1,x_e**2/2)
		    x_pdf=(omega_div_omega_ce-cmplx(n,0))/k_para_rho_e_para
		    call pdf(x_pdf,y_pdf)
		   
		    series_sum_1=series_sum_1+omega_pe_div_omega_ce/omega_div_omega_ce*(gamma_n*(x_e**4-2*n*x_e**2+2*n**2)-x_e**4*gamma_n1)*(1/x_e**2/k_para_rho_e_para*y_pdf)
		    if(n>0) then
		      x_pdf=(omega_div_omega_ce+cmplx(n,0))/k_para_rho_e_para
		      call pdf(x_pdf,y_pdf)
		      series_sum_1=series_sum_1+omega_pe_div_omega_ce/omega_div_omega_ce*(gamma_n*(x_e**4-2*n*x_e**2+2*n**2)-x_e**4*gamma_n1)*(1/x_e**2/k_para_rho_e_para*y_pdf)
		   
		    end if
		    
		   end if
		   
		   if(n>=10 .and. abs(series_sum_1-series_sum_2)<1e-8*abs(series_sum_1+series_sum_2) )then
		      exit
		   else if(n==50) then
		      write(*,*) 'The sum of the first 50 series does not converge'
		   end if
		   series_sum_2=series_sum_1
		end do
		if(abs(k_para_rho_i_para)<1e-6) then
		      D(2,2)=1-c_div_v_para*x**2/omega_div_omega_ci**2+series_sum_1
		else
		      D(2,2)=1-c_div_v_para*(k_para_rho_i_para**2+x**2)/omega_div_omega_ci**2+series_sum_1
		end if
        
	    case(6)
		series_sum_1=(0,0)
		series_sum_2=(0,0)
		if (abs(k_para_rho_i_para)<1e-6) then
		  D(2,3)=(0,0)
		
		else
		  
		  do n=0,10
            gamma_n=bessel_gn_complex(n,x**2/2)
            gamma_n1=bessel_gn_complex(n+1,x**2/2)
		    x_pdf=(omega_div_omega_ci-cmplx(n,0))/k_para_rho_i_para
		    call pdf(x_pdf,y_pdf)
			      
		    series_sum_1=series_sum_1-cmplx(0,omega_pi_div_omega_ci,wp)/omega_div_omega_ci*x/k_para_rho_i_para*(x**2*(gamma_n1-gamma_n)+2*n*gamma_n)*(1/x**2)*(1+x_pdf*y_pdf)
		    if(n>0) then
			  x_pdf=(omega_div_omega_ci+cmplx(n,0))/k_para_rho_i_para
			  call pdf(x_pdf,y_pdf)
			  series_sum_1=series_sum_1-cmplx(0,omega_pi_div_omega_ci,wp)/omega_div_omega_ci*x/k_para_rho_i_para*(x**2*(gamma_n1-gamma_n)+2*n*gamma_n)*(1/x**2)*(1+x_pdf*y_pdf)
		    end if
            gamma_n=bessel_gn_complex(n,x_e**2/2)
            gamma_n1=bessel_gn_complex(n+1,x_e**2/2)
		    x_pdf=(omega_div_omega_ce-cmplx(n,0))/k_para_rho_e_para
		    call pdf(x_pdf,y_pdf)
			      
		    series_sum_1=series_sum_1-cmplx(0,omega_pe_div_omega_ce,wp)/omega_div_omega_ce*x_e/k_para_rho_e_para*(x_e**2*(gamma_n1-gamma_n)+2*n*gamma_n)*(1/x_e**2)*(1+x_pdf*y_pdf)
		    if(n>0) then
			  x_pdf=(omega_div_omega_ce+cmplx(n,0))/k_para_rho_e_para
			  call pdf(x_pdf,y_pdf)
			  series_sum_1=series_sum_1-cmplx(0,omega_pe_div_omega_ce,wp)/omega_div_omega_ce*x_e/k_para_rho_e_para*(x_e**2*(gamma_n1-gamma_n)+2*n*gamma_n)*(1/x_e**2)*(1+x_pdf*y_pdf)
		      
		    end if
		     if(n>=10 .and. abs(series_sum_1-series_sum_2)<1e-8*abs(series_sum_1+series_sum_2) )then
				exit
		     else if(n==50) then
				write(*,*) 'The sum of the first 50 series does not converge'
		     end if
		     series_sum_2=series_sum_1
		  end do
		  D(2,3)=series_sum_1
		end if
        case(7)
            D(3,1)=D(1,3)
        case(8)
           D(3,2)=-D(2,3)
        case(9)
            
		series_sum_1=(0,0)
		series_sum_2=(0,0)
		do n=0,10
		  
		  if (abs(k_para_rho_i_para)<1e-6) then
            gamma_n=bessel_gn_complex(n,x**2/2)
		    series_sum_1=series_sum_1+omega_pi_div_omega_ci/omega_div_omega_ci*gamma_n/(cmplx(n,0)-omega_div_omega_ci)
		    
		    if (n>0) then
		      series_sum_1=series_sum_1+omega_pi_div_omega_ci/omega_div_omega_ci*gamma_n/(cmplx(-n,0)-omega_div_omega_ci)
		    end if
		    gamma_n=bessel_gn_complex(n,x_e**2/2)
		    series_sum_1=series_sum_1+omega_pe_div_omega_ce/omega_div_omega_ce*gamma_n/(cmplx(n,0)-omega_div_omega_ce)
		    
		    if (n>0) then
		      series_sum_1=series_sum_1+omega_pe_div_omega_ce/omega_div_omega_ce*gamma_n/(cmplx(-n,0)-omega_div_omega_ce)
		    end if
          else
            gamma_n=bessel_gn_complex(n,x**2/2)
		    x_pdf=(omega_div_omega_ci-cmplx(n,0))/k_para_rho_i_para
		    call pdf(x_pdf,y_pdf)
		    series_sum_1=series_sum_1+2*omega_pi_div_omega_ci*x**2/omega_div_omega_ci/k_para_rho_i_para*x_pdf*gamma_n*(1/x**2)*(1+x_pdf*y_pdf)
		      
		    if (n>0) then
		      x_pdf=(omega_div_omega_ci+cmplx(n,0))/k_para_rho_i_para
		      call pdf(x_pdf,y_pdf)
		      series_sum_1=series_sum_1+2*omega_pi_div_omega_ci*x**2/omega_div_omega_ci/k_para_rho_i_para*x_pdf*gamma_n*(1/x**2)*(1+x_pdf*y_pdf)
		    end if
		    gamma_n=bessel_gn_complex(n,x_e**2/2)
		    x_pdf=(omega_div_omega_ce-cmplx(n,0))/k_para_rho_e_para
		    call pdf(x_pdf,y_pdf)
		    series_sum_1=series_sum_1+2*omega_pe_div_omega_ce*x_e**2/omega_div_omega_ce/k_para_rho_e_para*x_pdf*gamma_n*(1/x_e**2)*(1+x_pdf*y_pdf)
		      
		    if (n>0) then
		      x_pdf=(omega_div_omega_ce+cmplx(n,0))/k_para_rho_e_para
		      call pdf(x_pdf,y_pdf)
		      series_sum_1=series_sum_1+2*omega_pe_div_omega_ce*x_e**2/omega_div_omega_ce/k_para_rho_e_para*x_pdf*gamma_n*(1/x_e**2)*(1+x_pdf*y_pdf)
		    end if
		  end if
		  
		  if(n>=10 .and. abs(series_sum_1-series_sum_2)<1e-8*abs(series_sum_1+series_sum_2) )then
		      exit
		  else if(n==50) then
		      write(*,*) 'The sum of the first 50 series does not converge'
		  end if
		  series_sum_2=series_sum_1
	      end do
	      D(3,3)=1-x**2*c_div_v_para/omega_div_omega_ci**2+series_sum_1
	end select
	    
    end do
    
    end subroutine dispersion_function_matrix_variable_k_per
    
    
!-----------------------------------------------------------------------------!
!     electron_dispersion_function_matrix_k_per:the dispersion relation for k_perp only considering electron 
!-----------------------------------------------------------------------------!
    subroutine electron_dispersion_function_matrix_variable_k_per(x,D)
    implicit none
    complex(wp),intent(in)::x
    complex(wp),intent(out)::D(:,:)
    complex(wp)::series_sum_1,series_sum_2,x_pdf,y_pdf
    complex(wp)::xe_pdf_1,xe_pdf_2,ye_pdf_1,ye_pdf_2
    integer::n,k
    complex(wp)::gamma_n,gamma_n1
    complex(wp)::x_e
	
	x_e=-x/(1836)**(0.5)
	do k=1,9
	    select case(k)
	    case(1)
		series_sum_1=(0,0)
		series_sum_2=(0,0)
		do n=1,10
            
		    gamma_n1=bessel_gn_complex(n,x_e**2/2)
		    if(abs(k_para_rho_i_para)<1e-6) then
		      series_sum_1=series_sum_1+4/x_e**2*omega_pe_div_omega_ce*gamma_n1*(n**2/(cmplx(n**2,0)-omega_div_omega_ce**2))
		    else
		      x_pdf=(omega_div_omega_ce-cmplx(n,0))/k_para_rho_e_para
		      call pdf(x_pdf,y_pdf)
		      
		      series_sum_1=series_sum_1+2*omega_pe_div_omega_ce/omega_div_omega_ce*gamma_n1*(n**2/x_e**2/k_para_rho_e_para*y_pdf)
		      
		      x_pdf=(omega_div_omega_ce+cmplx(n,0))/k_para_rho_e_para
		      call pdf(x_pdf,y_pdf)
		    
		      series_sum_1=series_sum_1+2*omega_pe_div_omega_ce/omega_div_omega_ce*gamma_n1*(n**2/x_e**2/k_para_rho_e_para*y_pdf)
		      
		    end if
		    
		    if(n>=10 .and. abs(series_sum_1-series_sum_2)<1e-8*abs(series_sum_1+series_sum_2) )then
				exit
            else if(n==50) then
				write(*,*) 'The sum of the first 50 series does not converge'
		    end if
		    series_sum_2=series_sum_1
		end do
		
		if (abs(k_para_rho_i_para)<1e-6) then
		    D(1,1)=(1,0)+series_sum_1
		else
		    D(1,1)=(1,0)-c_div_v_para*k_para_rho_i_para**2/omega_div_omega_ci**2+series_sum_1
		end if
        
	  case(2)
		series_sum_1=(0,0)
		series_sum_2=(0,0)
		do n=1,10
		    if(abs(k_para_rho_i_para)<1e-6) then
              gamma_n=bessel_gn_complex(n,x_e**2/2)
              gamma_n1=bessel_gn_complex(n+1,x_e**2/2)
		      series_sum_1=series_sum_1+(0,1)*omega_pe_div_omega_ce/omega_div_omega_ce*(2*n**2)/(cmplx(n**2,0)-omega_div_omega_ce**2)*(gamma_n1-(1-2*n/x_e**2)*gamma_n)
		    else
              gamma_n=bessel_gn_complex(n,x_e**2/2)
              gamma_n1=bessel_gn_complex(n+1,x_e**2/2)
		      x_pdf=(omega_div_omega_ce-cmplx(n,0))/k_para_rho_e_para
		      call pdf(x_pdf,y_pdf)
		      
		      series_sum_1=series_sum_1+cmplx(0,omega_pe_div_omega_ce,wp)/omega_div_omega_ce*(x_e**2*(gamma_n1-gamma_n)+2*n*gamma_n)*(n/x_e**2/k_para_rho_e_para*y_pdf)
		    
		    
		      x_pdf=(omega_div_omega_ce+cmplx(n,0))/k_para_rho_e_para
		      call pdf(x_pdf,y_pdf)
		   
		      series_sum_1=series_sum_1-cmplx(0,omega_pe_div_omega_ce,wp)/omega_div_omega_ce*(x_e**2*(gamma_n1-gamma_n)+2*n*gamma_n)*(n/x_e**2/k_para_rho_e_para*y_pdf)

		    end if
		    
		    if(n>=10 .and. abs(series_sum_1-series_sum_2)<1e-8*abs(series_sum_1+series_sum_2) )then
				exit
            else if(n==50) then
				write(*,*) 'The sum of the first 50 series does not converge'
		    end if
		    series_sum_2=series_sum_1
		end do
			if (abs(k_para_rho_i_para)<1e-6) then
				D(1,2)=series_sum_1
			else
				D(1,2)=series_sum_1
			end if
          
	    case(3)
		series_sum_1=(0,0)
		series_sum_2=(0,0)
		if(abs(k_para_rho_i_para)<1e-6) then
		  D(1,3)=(0,0)
		else
		   do n=1,10
              gamma_n=bessel_gn_complex(n,x_e**2/2)
		      x_pdf=(omega_div_omega_ce-cmplx(n,0))/k_para_rho_e_para
		      call pdf(x_pdf,y_pdf)
		      series_sum_1=series_sum_1+2*omega_pe_div_omega_ce/omega_div_omega_ce*n*gamma_n*x_e/k_para_rho_e_para*(1/x_e**2)*(1+x_pdf*y_pdf)
			
		      x_pdf=(omega_div_omega_ce+cmplx(n,0))/k_para_rho_e_para
		      call pdf(x_pdf,y_pdf)
		      series_sum_1=series_sum_1-2*omega_pe_div_omega_ce/omega_div_omega_ce*n*gamma_n*x_e/k_para_rho_e_para*(1/x_e**2)*(1+x_pdf*y_pdf)
		      
		      if(n>=10 .and. abs(series_sum_1-series_sum_2)<1e-8*abs(series_sum_1+series_sum_2) )then
				  exit
              else if(n==50) then
				  write(*,*) 'The sum of the first 50 series does not converge'
		      end if
		      series_sum_2=series_sum_1
		   end do
		      D(1,3)=c_div_v_para*k_para_rho_i_para*x/omega_div_omega_ci**2+series_sum_1
		    
		end if
	    case(4)
			D(2,1)=-D(1,2)
        case(5)
            
		series_sum_1=(0,0)
		series_sum_2=(0,0)
		do n=0,10
		  
		  if(abs(k_para_rho_i_para)<1e-6) then
            gamma_n=bessel_gn_complex(n,x_e**2/2)
			gamma_n1=bessel_gn_complex(n+1,x_e**2/2)
		    series_sum_1=series_sum_1+omega_pe_div_omega_ce/omega_div_omega_ce*((x_e**2-2*n+2*n**2/x_e**2)*gamma_n-x_e**2*gamma_n1)/(cmplx(n,0)-omega_div_omega_ce)
		    if(n>0) then
		      series_sum_1=series_sum_1+omega_pe_div_omega_ce/omega_div_omega_ce*((x_e**2-2*n+2*n**2/x_e**2)*gamma_n-x_e**2*gamma_n1)/(cmplx(-n,0)-omega_div_omega_ce)
		    end if
		    
		  else
            gamma_n=bessel_gn_complex(n,x_e**2/2)
			gamma_n1=bessel_gn_complex(n+1,x_e**2/2)
		    x_pdf=(omega_div_omega_ce-cmplx(n,0))/k_para_rho_e_para
		    call pdf(x_pdf,y_pdf)
		   
		    series_sum_1=series_sum_1+omega_pe_div_omega_ce/omega_div_omega_ce*(gamma_n*(x_e**4-2*n*x_e**2+2*n**2)-x_e**4*gamma_n1)*(1/x_e**2/k_para_rho_e_para*y_pdf)
		    if(n>0) then
		      x_pdf=(omega_div_omega_ce+cmplx(n,0))/k_para_rho_e_para
		      call pdf(x_pdf,y_pdf)
		      series_sum_1=series_sum_1+omega_pe_div_omega_ce/omega_div_omega_ce*(gamma_n*(x_e**4-2*n*x_e**2+2*n**2)-x_e**4*gamma_n1)*(1/x_e**2/k_para_rho_e_para*y_pdf)
		   
		    end if
		    
		   end if
		   
		   if(n>=10 .and. abs(series_sum_1-series_sum_2)<1e-8*abs(series_sum_1+series_sum_2) )then
		      exit
		   else if(n==50) then
		      write(*,*) 'The sum of the first 50 series does not converge'
		   end if
		   series_sum_2=series_sum_1
		end do
		if(abs(k_para_rho_i_para)<1e-6) then
		      D(2,2)=1-c_div_v_para*x**2/omega_div_omega_ci**2+series_sum_1
		else
		      D(2,2)=1-c_div_v_para*(k_para_rho_i_para**2+x**2)/omega_div_omega_ci**2+series_sum_1
		end if
        
	    case(6)
		series_sum_1=(0,0)
		series_sum_2=(0,0)
		if (abs(k_para_rho_i_para)<1e-6) then
		  D(2,3)=(0,0)
		
		else
		  
		  do n=0,10
            gamma_n=bessel_gn_complex(n,x_e**2/2)
            gamma_n1=bessel_gn_complex(n+1,x_e**2/2)
		    x_pdf=(omega_div_omega_ce-cmplx(n,0))/k_para_rho_e_para
		    call pdf(x_pdf,y_pdf)
			      
		    series_sum_1=series_sum_1-cmplx(0,omega_pe_div_omega_ce,wp)/omega_div_omega_ce*x_e/k_para_rho_e_para*(x_e**2*(gamma_n1-gamma_n)+2*n*gamma_n)*(1/x_e**2)*(1+x_pdf*y_pdf)
		    if(n>0) then
			  x_pdf=(omega_div_omega_ce+cmplx(n,0))/k_para_rho_e_para
			  call pdf(x_pdf,y_pdf)
			  series_sum_1=series_sum_1-cmplx(0,omega_pe_div_omega_ce,wp)/omega_div_omega_ce*x_e/k_para_rho_e_para*(x_e**2*(gamma_n1-gamma_n)+2*n*gamma_n)*(1/x_e**2)*(1+x_pdf*y_pdf)
		      
		    end if
		     if(n>=10 .and. abs(series_sum_1-series_sum_2)<1e-8*abs(series_sum_1+series_sum_2) )then
				exit
		     else if(n==50) then
				write(*,*) 'The sum of the first 50 series does not converge'
		     end if
		     series_sum_2=series_sum_1
		  end do
		  D(2,3)=series_sum_1
		end if
        case(7)
            D(3,1)=D(1,3)
        case(8)
           D(3,2)=-D(2,3)
        case(9)
            
		series_sum_1=(0,0)
		series_sum_2=(0,0)
		do n=0,10
		  
		  if (abs(k_para_rho_i_para)<1e-6) then
		    gamma_n=bessel_gn_complex(n,x_e**2/2)
		    series_sum_1=series_sum_1+omega_pe_div_omega_ce/omega_div_omega_ce*gamma_n/(cmplx(n,0)-omega_div_omega_ce)
		    
		    if (n>0) then
		      series_sum_1=series_sum_1+omega_pe_div_omega_ce/omega_div_omega_ce*gamma_n/(cmplx(-n,0)-omega_div_omega_ce)
		    end if
          else
		    gamma_n=bessel_gn_complex(n,x_e**2/2)
		    x_pdf=(omega_div_omega_ce-cmplx(n,0))/k_para_rho_e_para
		    call pdf(x_pdf,y_pdf)
		    series_sum_1=series_sum_1+2*omega_pe_div_omega_ce*x_e**2/omega_div_omega_ce/k_para_rho_e_para*x_pdf*gamma_n*(1/x_e**2)*(1+x_pdf*y_pdf)
		      
		    if (n>0) then
		      x_pdf=(omega_div_omega_ce+cmplx(n,0))/k_para_rho_e_para
		      call pdf(x_pdf,y_pdf)
		      series_sum_1=series_sum_1+2*omega_pe_div_omega_ce*x_e**2/omega_div_omega_ce/k_para_rho_e_para*x_pdf*gamma_n*(1/x_e**2)*(1+x_pdf*y_pdf)
		    end if
		  end if
		  
		  if(n>=10 .and. abs(series_sum_1-series_sum_2)<1e-8*abs(series_sum_1+series_sum_2) )then
		      exit
		  else if(n==50) then
		      write(*,*) 'The sum of the first 50 series does not converge'
		  end if
		  series_sum_2=series_sum_1
	      end do
	      D(3,3)=1-x**2*c_div_v_para/omega_div_omega_ci**2+series_sum_1
	end select
	    
    end do
    
    end subroutine electron_dispersion_function_matrix_variable_k_per
    
!-----------------------------------------------------------------------------!
!     electron_bernstein_wave: the dispersion relation of electron bernstein waves.
!-----------------------------------------------------------------------------!
    complex(wp) function electron_bernstein_wave(x)
	implicit none
	complex(wp),intent(in)::x
        complex(wp)::xi
        real(wp)::gamma_n,gamma_n1
        complex(wp)::series_sum_1,series_sum_2
        integer::n
        xi=x*1836
        series_sum_1=(0,0)
        series_sum_2=(0,0)
        do n=1,10
        
	  call bessel_gn(k_per_rho_i_per**2/2,n, gamma_n)
	  call bessel_gn(k_per_rho_e_per**2/2,n, gamma_n1)
	  !series_sum_1=series_sum_1+4/k_per_rho_i_per**2*omega_pi_div_omega_ci*gamma_n*(n**2/(cmplx(n**2,0)-xi**2))
	  series_sum_1=series_sum_1+4/k_per_rho_e_per**2*omega_pe_div_omega_ce*gamma_n1*(n**2/(cmplx(n**2,0)-x**2))
	  
	  if(n>=10 .and. abs(series_sum_1-series_sum_2)<1e-5*abs(series_sum_1+series_sum_2) )then
	      
	      exit
	  else if(n==50) then
	      write(*,*) 'The sum of the first 50 series does not converge'
	      write(*,*) k_per_rho_e_per
	      write(*,*) gamma_n,gamma_n1
	  end if
	  series_sum_2=series_sum_1
	end do
	electron_bernstein_wave=1+series_sum_1
    end function electron_bernstein_wave
    
!-----------------------------------------------------------------------------!
!     dispersion_function_parallel_1: the dispersion relation of right hand wave with omega as the varable.
!-----------------------------------------------------------------------------!
    complex(wp) function dispersion_function_parallel_1(x)
	implicit none
	complex(wp),intent(in)::x
	complex(wp)::D(3,3)
    complex(wp)::y
    y=x
	call dispersion_function_parallel_matrix(y,D)
	dispersion_function_parallel_1=D(1,1)+(0,1)*D(1,2)
    end function dispersion_function_parallel_1
       
!-----------------------------------------------------------------------------!
!     dispersion_function_parallel_2: the dispersion relation of left hand wave with omega as the varable.
!-----------------------------------------------------------------------------! 
    complex(wp) function dispersion_function_parallel_2(x)
	implicit none
	complex(wp),intent(in)::x
	complex(wp)::D(3,3)
    	complex(wp)::y
   		 y=x
	call dispersion_function_parallel_matrix(y,D)
	dispersion_function_parallel_2=D(1,1)-(0,1)*D(1,2)
    end function dispersion_function_parallel_2
    
!-----------------------------------------------------------------------------!
!     dispersion_function_parallel_3:the dispersion relation of langmuir wave and ion acoustic wave with omega as the varable.
!-----------------------------------------------------------------------------! 
    complex(wp) function dispersion_function_parallel_3(x)
	implicit none
	complex(wp),intent(in)::x
	complex(wp)::D(3,3)
    complex(wp)::y
    y=x
	call dispersion_function_parallel_matrix(y,D)
	dispersion_function_parallel_3=D(3,3)
    end function dispersion_function_parallel_3

!-----------------------------------------------------------------------------!
!     dispersion_function_parallel_variable_k_para_1: the dispersion relation of left hand wave with k_para as the varable.
!-----------------------------------------------------------------------------!
    complex(wp) function dispersion_function_parallel_variable_k_para_1(x)
	implicit none
	complex(wp),intent(in)::x
	complex(wp)::D(3,3)
    complex(wp)::y
    y=(1.0,1.0)*x/(2.0_wp)**(0.5)
	call dispersion_function_parallel_matrix_variable_k_para(y,D)
	dispersion_function_parallel_variable_k_para_1=D(1,1)+(0,1)*D(1,2)
    end function dispersion_function_parallel_variable_k_para_1
       
!-----------------------------------------------------------------------------!
!     dispersion_function_parallel_2:the dispersion relation of right hand wave with k_para as the varable.
!-----------------------------------------------------------------------------! 
    complex(wp) function dispersion_function_parallel_variable_k_para_2(x)
	implicit none
	complex(wp),intent(in)::x
	complex(wp)::D(3,3)
    	complex(wp)::y
		y=(1.0,1.0)*x/(2.0_wp)**(0.5)
		call dispersion_function_parallel_matrix_variable_k_para(y,D)
		dispersion_function_parallel_variable_k_para_2=D(1,1)-(0,1)*D(1,2)
    end function dispersion_function_parallel_variable_k_para_2
    
!-----------------------------------------------------------------------------!
!     dispersion_function_parallel_3:the dispersion relation of langmuir and ion acoustic wave with k_para as the varable.
!-----------------------------------------------------------------------------! 
    complex(wp) function dispersion_function_parallel_variable_k_para_3(x)
	implicit none
	complex(wp),intent(in)::x
	complex(wp)::D(3,3)
    complex(wp)::y
    y=(1.0,1.0)*x/(2.0_wp)**(0.5)
	call dispersion_function_parallel_matrix_variable_k_para(y,D)
	dispersion_function_parallel_variable_k_para_3=D(3,3)
    end function dispersion_function_parallel_variable_k_para_3    

!-----------------------------------------------------------------------------!
!     dispersion_function:dispersion relation with omega as the variable
!-----------------------------------------------------------------------------! 
    complex(wp) function dispersion_function(x)
	implicit none
	complex(wp),intent(in)::x
	complex(wp)::D(3,3)
	integer::n,k
	call dispersion_function_matrix(x,D)
	do n=1,3
	    do k=1,3
		D(n,k)=D(n,k)/1000
	    end do
	end do
	dispersion_function=det(D,3)
  
    end function dispersion_function

!-----------------------------------------------------------------------------!
!     dispersion_function_perpendicular:dispersion relation with omega as the variable for perpendicular waves
!-----------------------------------------------------------------------------! 
    complex(wp) function dispersion_function_perpendicular(x)
	implicit none
	complex(wp),intent(in)::x
	complex(wp)::D(3,3)
	integer::n,k
	call dispersion_function_matrix(x,D)
	do n=1,3
	    do k=1,3
		D(n,k)=D(n,k)/1000
	    end do
	end do
	dispersion_function_perpendicular=D(1,1)*D(2,2)-D(1,2)*D(2,1)
  
    end function dispersion_function_perpendicular
    
!-----------------------------------------------------------------------------!
!     dispersion_function_variable_k_para: dispersion relation with k_para as the variable
!-----------------------------------------------------------------------------! 
    complex(wp) function dispersion_function_variable_k_para(x)
	implicit none
	complex(wp),intent(in)::x
	complex(wp)::D(3,3)
	integer::n,k
    complex(wp)::y
    y=(1.0,1.0)*x/(2.0_wp)**(0.5)
	call dispersion_function_matrix_variable_k_para(y,D)
	do n=1,3
	    do k=1,3
		D(n,k)=D(n,k)
	    end do
	end do
	dispersion_function_variable_k_para=det(D,3)
  
    end function dispersion_function_variable_k_para

!-----------------------------------------------------------------------------!
!     dispersion_function_variable_k_per: dispersion relation with k_perp as the variable
!-----------------------------------------------------------------------------! 
    complex(wp) function dispersion_function_variable_k_per(x)
	implicit none
	complex(wp),intent(in)::x
	complex(wp)::D(3,3)
	integer::n,k
	complex(wp)::y
    y=x*omega_div_omega_ci/(c_div_v_para)**(0.5)
	call dispersion_function_matrix_variable_k_per(y,D)
	do n=1,3
	    do k=1,3
		D(n,k)=D(n,k)
	    end do
	end do
	dispersion_function_variable_k_per=det(D,3)
  
    end function dispersion_function_variable_k_per 
    
!-----------------------------------------------------------------------------!
!     electron_dispersion_function_variable_k_per: dispersion relation only considering electron with k_perp as the variable
!-----------------------------------------------------------------------------! 
    complex(wp) function electron_dispersion_function_variable_k_per(x)
	implicit none
	complex(wp),intent(in)::x
	complex(wp)::D(3,3)
	integer::n,k
	complex(wp)::y
    y=x*omega_div_omega_ci/(c_div_v_para)**(0.5)
	call electron_dispersion_function_matrix_variable_k_per(y,D)
	do n=1,3
	    do k=1,3
		D(n,k)=D(n,k)
	    end do
	end do
	electron_dispersion_function_variable_k_per=det(D,3)
  
    end function electron_dispersion_function_variable_k_per 
    
!-----------------------------------------------------------------------------!
!     dispersion_function_parallel_variable_k_para:dispersion relation for parallel waves with k_para as the variable
!-----------------------------------------------------------------------------! 
    complex(wp) function dispersion_function_parallel_variable_k_para(x)
	implicit none
	complex(wp),intent(in)::x
	complex(wp)::D(3,3)
	integer::n,k
    complex(wp)::y
    y=x*omega_div_omega_ci/(c_div_v_para)**(0.5)
    y=(1.0,1.0)*y/(2.0_wp)**(0.5)
	call dispersion_function_parallel_matrix_variable_k_para(y,D)
	do n=1,3
	    do k=1,3
		D(n,k)=D(n,k)
	    end do
	end do
	dispersion_function_parallel_variable_k_para=det(D,3)
  
    end function dispersion_function_parallel_variable_k_para   
    
!-----------------------------------------------------------------------------!
!     cold_plasma_dispersion_function:dispersion relation for cold plasma with omega as the variable
!-----------------------------------------------------------------------------! 
    subroutine cold_plasma_dispersion_function_matrix(x,D)
    implicit none
    complex(wp),intent(in)::x
    complex(wp),intent(out)::D(:,:)
    complex(wp)::x_e,ka,k_per,k_para

        x_e=-x/1836
        k_per=1-omega_pi_div_omega_ci/(x**2-1)-omega_pe_div_omega_ce/(x_e**2-1)
        ka=omega_pi_div_omega_ci/(x**2-1)/x+omega_pe_div_omega_ce/(x_e**2-1)/x_e
        k_para=1-omega_pi_div_omega_ci/x**2-omega_pe_div_omega_ce/x_e**2
        D(1,1)=c_div_v_para*k_para_rho_i_para**2/x**2-k_per
        D(1,2)=(0,1)*ka
        D(1,3)=-c_div_v_para*k_para_rho_i_para*k_per_rho_i_para/x**2
        D(2,1)=-D(1,2)
        D(2,2)=c_div_v_para*(k_para_rho_i_para**2+k_per_rho_i_para**2)/x**2-k_per
        D(2,3)=(0,0)
        D(3,1)=D(1,3)
        D(3,2)=(0,0)
        D(3,3)=c_div_v_para*k_per_rho_i_para**2/x**2-k_para
    
    end subroutine cold_plasma_dispersion_function_matrix
    
    
    complex(wp) function cold_plasma_dispersion_function(x)
	implicit none
	complex(wp),intent(in)::x
	complex(wp)::x_e,ka,k_per,k_para
	complex(wp)::D(3,3)
	integer::n,k
        x_e=-x/1836
        k_per=1-omega_pi_div_omega_ci/(x**2-1)-omega_pe_div_omega_ce/(x_e**2-1)
        ka=omega_pi_div_omega_ci/(x**2-1)/x+omega_pe_div_omega_ce/(x_e**2-1)/x_e
        k_para=1-omega_pi_div_omega_ci/x**2-omega_pe_div_omega_ce/x_e**2
        D(1,1)=c_div_v_para*k_para_rho_i_para**2/x**2-k_per
        D(1,2)=(0,1)*ka
        D(1,3)=-c_div_v_para*k_para_rho_i_para*k_per_rho_i_para/x**2
        D(2,1)=-D(1,2)
        D(2,2)=c_div_v_para*(k_para_rho_i_para**2+k_per_rho_i_para**2)/x**2-k_per
        D(2,3)=(0,0)
        D(3,1)=D(1,3)
        D(3,2)=0
        D(3,3)=c_div_v_para*k_per_rho_i_para**2/x**2-k_para
	!do n=1,3
	!	do k=1,3
	!		D(n,k)=D(n,k)/1000
!      
	!	end do
	!end do
	cold_plasma_dispersion_function=det(D,3)
	cold_plasma_dispersion_function=log(cold_plasma_dispersion_function)
    end function cold_plasma_dispersion_function

!-----------------------------------------------------------------------------!
!     cold_plasma_dispersion_function_two_ion_species:dispersion relation for cold plasma with omega as the variable with two ion species.
!-----------------------------------------------------------------------------! 
	subroutine cold_plasma_dispersion_function_two_ion_species_matrix(x,D)
		implicit none
		complex(wp),intent(in)::x
		complex(wp),intent(out)::D(:,:)
		complex(wp)::x_e,x_d,ka,k_per,k_para
	
			x_e=-x/1836
			x_d=2*x
			k_per=1-omega_pi_div_omega_ci/(x**2-1)-omega_pe_div_omega_ce/(x_e**2-1)-omega_pd_div_omega_cd/(x_d**2-1)
			ka=omega_pi_div_omega_ci/(x**2-1)/x+omega_pe_div_omega_ce/(x_e**2-1)/x_e+omega_pd_div_omega_cd/(x_d**2-1)/x_d
			k_para=1-omega_pi_div_omega_ci/x**2-omega_pe_div_omega_ce/x_e**2-omega_pd_div_omega_cd/x_d**2
			D(1,1)=k_para_c_div_omega_ci**2/x**2-k_per
			D(1,2)=(0,1)*ka
			D(1,3)=-k_para_c_div_omega_ci*k_per_c_div_omega_ci/x**2
			D(2,1)=-D(1,2)
			D(2,2)=(k_para_c_div_omega_ci**2+k_per_c_div_omega_ci**2)/x**2-k_per
			D(2,3)=(0,0)
			D(3,1)=D(1,3)
			D(3,2)=(0,0)
			D(3,3)=k_per_c_div_omega_ci**2/x**2-k_para
		
	end subroutine cold_plasma_dispersion_function_two_ion_species_matrix
		
		
	complex(wp) function cold_plasma_dispersion_function_two_ion_species(x)
		implicit none
		complex(wp),intent(in)::x
		complex(wp)::x_e,x_d,ka,k_per,k_para
		complex(wp)::D(3,3)
		integer::n,k
		call cold_plasma_dispersion_function_two_ion_species_matrix(x,D)
		cold_plasma_dispersion_function_two_ion_species=det(D,3)
	end function cold_plasma_dispersion_function_two_ion_species
	
!-----------------------------------------------------------------------------!
!     cold_plasma_dispersion_function_k_per:dispersion relation for cold plasma with k_perp as the variable
!-----------------------------------------------------------------------------! 
    subroutine cold_plasma_dispersion_function_matrix_k_per(x,D)
    implicit none
    complex(wp),intent(in)::x
    complex(wp),intent(out)::D(:,:)
    complex(wp)::ka,k_per,k_para

    k_per=1-omega_pi_div_omega_ci/(omega_div_omega_ci**2-1)-omega_pe_div_omega_ce/(omega_div_omega_ce**2-1)
    ka=omega_pi_div_omega_ci/(omega_div_omega_ci**2-1)/omega_div_omega_ci+omega_pe_div_omega_ce/(omega_div_omega_ce**2-1)/omega_div_omega_ce
    k_para=1-omega_pi_div_omega_ci/omega_div_omega_ci**2-omega_pe_div_omega_ce/omega_div_omega_ce**2
    D(1,1)=refractive_para**2-k_per
    D(1,2)=(0,1)*ka
    D(1,3)=-refractive_para*x
    D(2,1)=-D(1,2)
    D(2,2)=refractive_para**2+x**2-k_per
    D(2,3)=(0,0)
    D(3,1)=D(1,3)
    D(3,2)=(0,0)
    D(3,3)=x**2-k_para
    
    end subroutine cold_plasma_dispersion_function_matrix_k_per
    
    
    complex(wp) function cold_plasma_dispersion_function_k_per(x)
	implicit none
	complex(wp),intent(in)::x
	complex(wp)::ka,k_per,k_para
	complex(wp)::D(3,3)
	integer::n,k
    
    k_per=1-omega_pi_div_omega_ci/(omega_div_omega_ci**2-1)-omega_pe_div_omega_ce/(omega_div_omega_ce**2-1)
    ka=omega_pi_div_omega_ci/(omega_div_omega_ci**2-1)/omega_div_omega_ci+omega_pe_div_omega_ce/(omega_div_omega_ce**2-1)/omega_div_omega_ce
    k_para=1-omega_pi_div_omega_ci/omega_div_omega_ci**2-omega_pe_div_omega_ce/omega_div_omega_ce**2
    D(1,1)=refractive_para**2-k_per
    D(1,2)=(0,1)*ka
    D(1,3)=-refractive_para*x
    D(2,1)=-D(1,2)
    D(2,2)=refractive_para**2+x**2-k_per
    D(2,3)=(0,0)
    D(3,1)=D(1,3)
    D(3,2)=0
    D(3,3)=x**2-k_para
	cold_plasma_dispersion_function_k_per=det(D,3)
    end function cold_plasma_dispersion_function_k_per

!-----------------------------------------------------------------------------!
!     left_cut_off:calculate the left cut off line (omega_pe_div_omega_ce_2)
!-----------------------------------------------------------------------------!  
    real(wp) function left_cut_off(omega_ce_div_x)
		implicit none
		real(wp),intent(in)::omega_ce_div_x
		real(wp)::mass_ratio
		mass_ratio=1836.0_wp
		left_cut_off=(1-omega_ce_div_x/mass_ratio)*(1+omega_ce_div_x)/(1+1/mass_ratio)
    end function left_cut_off
 
!-----------------------------------------------------------------------------!
!     right_cut_off:calculate the right cut off line (omega_pe_div_omega_ce_2)
!-----------------------------------------------------------------------------!    
    real(wp) function right_cut_off(omega_ce_div_x)
		implicit none
		real(wp),intent(in)::omega_ce_div_x
		real(wp)::mass_ratio
		mass_ratio=1836.0_wp
		right_cut_off=(1+omega_ce_div_x/mass_ratio)*(1-omega_ce_div_x)/(1+1/mass_ratio)
    end function right_cut_off
 
!-----------------------------------------------------------------------------!
!     resonance:calculate the resonance line (omega_pe_div_omega_ce_2)
!-----------------------------------------------------------------------------!  
    real(wp) function resonance(omega_ce_div_x)
		implicit none
        real(wp),intent(in)::omega_ce_div_x
        real(wp)::mass_ratio
        mass_ratio=1836.0_wp
        resonance=1/(1/(1-omega_ce_div_x**2)+1/mass_ratio/(1-omega_ce_div_x**2/mass_ratio**2))
    end function resonance
    
 
!-----------------------------------------------------------------------------!
!     test_fun5:test function 5
!-----------------------------------------------------------------------------!  
    complex(wp) function test_fun5(x)
	implicit none
	complex(wp),intent(in)::x
        complex(wp)::x_pdf,y_pdf
        x_pdf=10*(x-(1,0))
        call pdf(x_pdf,y_pdf)
	test_fun5=1-4700/x**2+19000*y_pdf/x
    end function test_fun5
     
    complex(wp) function l_wave(x)
	implicit none
	complex(wp),intent(in)::x
        complex(wp)::x_pdf,y_pdf
        complex(wp)::y
        y=x/100
	!l_wave=1-c_div_v_para*k_para_rho_i_para**2/x/x-omega_pi_div_omega_ci/x/x/(1-1/x)
	l_wave=1-c_div_v_para*k_para_rho_i_para**2/y/y-omega_pi_div_omega_ci/y/y/(1-1/y)
    end function l_wave
    
!-----------------------------------------------------------------------------!
!     polarization: caculate the polarization of the plasma wave
!-----------------------------------------------------------------------------!  
    subroutine polarization(x,eigen,polar)
	implicit none
	complex(wp),intent(in)::x
	real(wp),intent(out)::eigen
	complex(wp),intent(out)::polar(:)
	complex(wp)::D(3,3),w(3),vl(3,3),vr(3,3)
	integer::n,k,info,e_min
	complex(wp)::y
	call cold_plasma_dispersion_function_two_ion_species_matrix(x,D)
    do n=1,3
	  do k=1,3
		  D(n,k)=D(n,k)/1000
  
	  end do
	end do
	call geev(D,w,vl,vr,info)
	eigen=abs(w(1))
	e_min=1
	do n=2,3
	  if (abs(w(n))<eigen) then
		  eigen=abs(w(n))
		  e_min=n
	  end if
	end do
	polar=vr(1:3,e_min)
    end subroutine polarization
    
    subroutine polarization_2(x,eigen,polar)
	implicit none
	complex(wp),intent(in)::x
	real(wp),intent(out)::eigen
	complex(wp),intent(out)::polar(:)
	complex(wp)::D(3,3),w(3),vl(3,3),vr(3,3)
	integer::n,k,info,e_min
	complex(wp)::y
    y=x*omega_div_omega_ci/(c_div_v_para)**(0.5)
	y=(1.0,1.0)*y/(2.0_wp)**(0.5)
	call  dispersion_function_parallel_matrix_variable_k_para(y,D)
   
    do n=1,3
	  do k=1,3
		  D(n,k)=D(n,k)
  
	  end do
	end do
	call geev(D,w,vl,vr,info)
	eigen=abs(w(1))
	e_min=1
	do n=2,3
	  if (abs(w(n))<eigen) then
		  eigen=abs(w(n))
		  e_min=n
	  end if
	end do
	polar=vr(1:3,e_min)
    end subroutine polarization_2
    
    
     subroutine cold_polarization(x,eigen,polar)
	implicit none
	complex(wp),intent(in)::x
	real(wp),intent(out)::eigen
	complex(wp),intent(out)::polar(:)
	complex(wp)::D(3,3),w(3),vl(3,3),vr(3,3)
	integer::n,k,info,e_min
	complex(wp)::y
    y=x
	call cold_plasma_dispersion_function_matrix_k_per(y,D)
!         do n=1,3
! 	  do k=1,3
! 		  D(n,k)=D(n,k)/1000
!   
! 	  end do
! 	end do
	call geev(D,w,vl,vr,info)
	eigen=abs(w(1))
	e_min=1
	do n=2,3
	  if (abs(w(n))<eigen) then
	      eigen=abs(w(n))
	      e_min=n
	  end if
	end do
	polar=vr(1:3,e_min)
    end subroutine cold_polarization
!-----------------------------------------------------------------------------!
!     damping_rate:calculate the damping rate with an approximate expression
!-----------------------------------------------------------------------------!     
    complex(wp) function damping_rate(x)
	implicit none
	complex(wp),intent(in)::x
        complex(wp)::x_pdf,y_pdf
        
        x_pdf=(x-(1,0))/k_para_rho_i_para
        call pdf(x_pdf,y_pdf)
        y_pdf=cmplx(real(y_pdf),0,wp)
		damping_rate=-(acos(-1.0_wp))**0.5/k_para_rho_i_para*omega_pi_div_omega_ci*exp(-x_pdf**2)/(1+c_div_v_para*k_para_rho_i_para**2/x**2+2*omega_pi_div_omega_ci/k_para_rho_i_para**2*(1+x_pdf*y_pdf))
    
    end function damping_rate
    
    
end module tool


