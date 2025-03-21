MODULE MATHFUNCTION
      use iso_fortran_env,only:wp=>real64
      use lapack95
      use quadpack_generic
      
      
      private x_integral,n_integral,z_integral,z1_integral
      real(wp)::x_integral
      integer::n_integral
      complex(wp)::z_integral,z1_integral
      CONTAINS
!-----------------------------------------------------------------------------!
! Here, we porvide some useful functions which Fortran doesn't have originally
!-----------------------------------------------------------------------------!
!-----------------------------------------------------------------------------!
!     HERMITE_FUNCTION
!-----------------------------------------------------------------------------!
         SUBROUTINE HERMITE_FUNCTION(alpha, n, x, result)

         IMPLICIT NONE
         REAL*8    ::alpha, x, result
         INTEGER   ::n
         REAL*8    ::h0, h1
         INTEGER   ::j

         h0 = 1.d0
         h1 = 2.d0*alpha*x

         IF(n .EQ. 0)THEN
             result = h0
         ELSE IF(n .EQ. 1)THEN
             result = h1
         ELSE IF(n .GT. 1)THEN
             DO j = 1, n-1
                  result  = 2.d0*alpha*x*h1 - 2.d0*j*h0
                  h0      = h1
                  h1      = result
             END DO
         END IF

         END SUBROUTINE HERMITE_FUNCTION


!-----------------------------------------------------------------------------!
!     FACTORIAL
!-----------------------------------------------------------------------------!
         SUBROUTINE FACTORIAL(n, result)

         IMPLICIT NONE
         INTEGER   ::n
         REAL*8    ::result      
         INTEGER   ::j

         IF(n .EQ. 0)THEN
             result = 1.d0
         ELSE IF (n .EQ. 1)THEN
             result = 1.d0
         ELSE IF (n .GT. 1)THEN
             result = 1.d0
             DO j = 1, n
                  result = result*dfloat(j)
             END DO
         END IF

         END SUBROUTINE FACTORIAL


!-----------------------------------------------------------------------------!
!     HERMITE_GAUSSIAN_FUNCTION
!-----------------------------------------------------------------------------!
         SUBROUTINE HERMITE_GAUSSIAN_FUNCTION(alpha, n, x, result, nmax)

         IMPLICIT NONE
         REAL*8    ::alpha, x, result
         INTEGER   ::n, nmax
         REAL*8    ::y1, y2, y3
         REAL*8    ::pi

         pi     = dacos(-1.d0)

         IF((n .GE. 0) .AND. (n .LE. nmax))THEN
            CALL FACTORIAL(n, y1)
            CALL HERMITE_FUNCTION(alpha, n, x, y2)
            y3     = EXP(- 0.5d0*alpha*alpha*x*x)

            result = (1.d0/SQRT((2**n)*y1*SQRT(pi)))*SQRT(alpha)*y3*y2
         ELSE IF(n .LT. 0) THEN
            result = 0.d0
         ELSE IF(n .GE. nmax) THEN
            result = 0.d0
         END IF

         END SUBROUTINE HERMITE_GAUSSIAN_FUNCTION


!-----------------------------------------------------------------------------!
!     WELL
!-----------------------------------------------------------------------------!
         SUBROUTINE WELL(omega, x, well_p, well_q, well_r)

         IMPLICIT NONE
         REAL*8    ::x
         COMPLEX*16::omega, well_p, well_q, well_r
         REAL*8    ::pi
         COMPLEX*16::cj
         REAL*8    ::DG_FLR, DG_FLRP, DG_FLRPP
         COMPLEX*16::MM_R, MM_Q, MM_P, MM_Zi, MM_Ze
         REAL*8    ::MM_e2i, MM_Tperp_i2e, MM_Tpara_i2e, MM_Ti_perp2para, MM_Te_perp2para, MM_rhoTperp_e2i
         REAL*8    ::MM_a_i, MM_a_e, MM_rhoTiperp2L, MM_rhoTeperp2L, MM_vTpara_i2e, MM_beta_i, MM_beta_e, MM_omega_i0, MM_omega_e0
         REAL*8    ::MM_kx, MM_kpara
         REAL*8    ::MM_C1, MM_C2i, MM_C2e, DG_i, DG_e, DGP_i, DGP_e, DGPP_i, DGPP_e

         cj               = (0.d0, 1.d0)               !Complex unit
         pi               = dacos(-1.d0)               !Pi

         ! ! Weber Equation
         ! well_p = 1.d0
         ! well_q = 0.d0
         ! well_r = - 0.25d0*x*x + omega + 0.5d0

         ! !Harmonic Oscillator
         ! well_p = 1.d0
         ! well_q = 0.d0
         ! well_r = - x*x + omega 

         ! !1D MM

         MM_e2i          = 1.d0/1836.d0
         MM_Tperp_i2e    = 10.d0
         MM_Tpara_i2e    = 10.d0
         MM_Ti_perp2para = 2.4d0
         MM_Te_perp2para = 2.4d0
         MM_beta_i       = 1.21d0
         MM_beta_e       = MM_beta_i*(1.d0/MM_Tperp_i2e)
         MM_rhoTperp_e2i = -dsqrt((1.d0/MM_Tperp_i2e)*MM_e2i)
         MM_rhoTiperp2L  = 1.d0/3.d0
         MM_rhoTeperp2L  = MM_rhoTiperp2L*MM_rhoTperp_e2i
         MM_vTpara_i2e   = dsqrt(MM_Tpara_i2e*MM_e2i)
         MM_kx           = 0.1d0
         MM_kpara        = 1.d0
         MM_a_i          = (1.d0/MM_Ti_perp2para) - 1.d0
         MM_a_e          = (1.d0/MM_Te_perp2para) - 1.d0
         MM_omega_i0     = -(MM_kx/MM_kpara)*(1.d0/(MM_beta_i + MM_beta_e))
         MM_omega_e0     = MM_omega_i0*MM_rhoTperp_e2i
         MM_C1           = (1.d0/(MM_beta_i+MM_beta_e)) + 1.d0
         MM_C2i          = MM_beta_i/(1.d0 + MM_a_i)
         MM_C2e          = MM_beta_e/(1.d0 + MM_a_e)

         CALL PDF(omega, MM_Zi)
         CALL PDF(omega*MM_vTpara_i2e, MM_Ze)

         ! MM_Zi           = cj*sqrt(pi)
         ! MM_Ze           = cj*sqrt(pi)

         CALL BESSEL_GFLR(0.5d0*(MM_kx*MM_kx + x*x), DG_i, DGP_i, DGPP_i)
         CALL BESSEL_GFLR(0.5d0*(MM_kx*MM_kx + x*x)*MM_rhoTperp_e2i*MM_rhoTperp_e2i, DG_e, DGP_e, DGPP_e)
         ! DG_i            = 1.d0/(1.d0 + 0.75d0*(MM_kx*MM_kx + x*x))
         ! DG_e            = 1.d0 - 0.75d0*(MM_kx*MM_kx + x*x)*MM_rhoTperp_e2i*MM_rhoTperp_e2i

         ! MM_R = - 3.d0*MM_C2i*(MM_C1*(omega*MM_Zi - MM_a_i) - MM_omega_i0*MM_Zi)*MM_rhoTiperp2L*MM_rhoTiperp2L*DG_i

         ! MM_Q = + 2.d0*MM_C2i*(MM_C1*(omega*MM_Zi - MM_a_i) - MM_omega_i0*MM_Zi)*cj*MM_rhoTiperp2L*DG_i

         ! MM_P = - 1.d0 + MM_C2i*((omega-MM_omega_i0)*MM_Zi - MM_a_i)*DG_i

         MM_R = - 3.d0*MM_C2i*(MM_C1*(omega*MM_Zi - MM_a_i) - MM_omega_i0*MM_Zi)*MM_rhoTiperp2L*MM_rhoTiperp2L*DG_i&
                - 3.d0*MM_C2e*(MM_C1*(omega*MM_vTpara_i2e*MM_Ze - MM_a_e) - MM_omega_e0*MM_Ze)*MM_rhoTiperp2L*MM_rhoTiperp2L*DG_e

         MM_Q = + 2.d0*MM_C2i*(MM_C1*(omega*MM_Zi - MM_a_i) - MM_omega_i0*MM_Zi)*cj*MM_rhoTiperp2L*DG_i&
                + 2.d0*MM_C2e*(MM_C1*(omega*MM_vTpara_i2e*MM_Ze - MM_a_e) - MM_omega_e0*MM_Ze)*MM_rhoTiperp2L*MM_rhoTiperp2L*DG_e

         MM_P = - 1.d0 + MM_C2i*((omega-MM_omega_i0)*MM_Zi - MM_a_i)*DG_i&
                       + MM_C2e*((omega*MM_vTpara_i2e-MM_omega_e0)*MM_Ze - MM_a_e)*DG_e

         well_p = MM_P
         well_q = MM_Q
         well_r = MM_R

         END SUBROUTINE WELL


!-----------------------------------------------------------------------------!
!     DETERMINANT
!-----------------------------------------------------------------------------!
         RECURSIVE SUBROUTINE GETDET(m, n, d)
        
         IMPLICIT NONE
         COMPLEX*16::m(:,:)
         INTEGER   ::n
         COMPLEX*16::d
         INTEGER   ::k,i,j
         COMPLEX*16,ALLOCATABLE::b(:,:)
         COMPLEX*16,ALLOCATABLE::d1(:)
        
         d=0

         IF (n .EQ. 1) THEN
            d=m(1,1)
         ELSE
            ALLOCATE(d1(n))
            ALLOCATE(b(n-1,n-1))
            DO k=1, n
               DO i=2, n
                  DO j=1, n
                     IF (j .LT. k) THEN
                        b(i-1,j)=m(i,j)
                     ELSE IF (j .GT. k) THEN
                        b(i-1,j-1)=m(i,j)   
                     END IF
                  END DO
               END DO
               CALL GETDET(b, n-1, d1(k))
               d=d+m(1,k)*(-1)**(k+1)*d1(k)
            END DO
            DEALLOCATE(b)
            DEALLOCATE(d1)
         END IF

         END SUBROUTINE GETDET


!=======================================================================
!    SUBROUTINE BESSEL
!
!         (1)calculate the 0th-order modified bessel function \Gamma_0s(b_s) = I_0(b_s)*exp(-b_s)
!=======================================================================
         SUBROUTINE BESSEL_G0(x, DG0)

         IMPLICIT NONE
         INTEGER   ::i, k
         REAL*8    ::a(0:64), b(0:69), c(0:44), t, w, x, y, DG0, DG1, DG_FLR, DG_FLRp

         DATA (a(i), i = 0, 12) / &
         8.5246820682016865877d-11, 2.5966600546497407288d-9, 7.9689994568640180274d-8, 1.9906710409667748239d-6, 4.0312469446528002532d-5, &
         6.4499871606224265421d-4, 7.9012345761930579108d-3, 7.1111111109207045212d-2, 4.4444444444472490900d-1, 1.7777777777777532045d0, &
         4.0000000000000011182d0, 3.9999999999999999800d0, 1.0000000000000000001d0 /
         DATA (a(i), i = 13, 25) / &
         1.1520919130377195927d-10, 2.2287613013610985225d-9, 8.1903951930694585113d-8, 1.9821560631611544984d-6, 4.0335461940910133184d-5, &
         6.4495330974432203401d-4, 7.9013012611467520626d-3, 7.1111038160875566622d-2, 4.4444450319062699316d-1, 1.7777777439146450067d0, &
         4.0000000132337935071d0, 3.9999999968569015366d0, 1.0000000003426703174d0 /
         DATA (a(i), i = 26, 38) / &
         1.5476870780515238488d-10, 1.2685004214732975355d-9, 9.2776861851114223267d-8, 1.9063070109379044378d-6, 4.0698004389917945832d-5, &
         6.4370447244298070713d-4, 7.9044749458444976958d-3, 7.1105052411749363882d-2, 4.4445280640924755082d-1, 1.7777694934432109713d0, &
         4.0000055808824003386d0, 3.9999977081165740932d0, 1.0000004333949319118d0 /
         DATA (a(i), i = 39, 51) / &
         2.0675200625006793075d-10, -6.1689554705125681442d-10, 1.2436765915401571654d-7, 1.5830429403520613423d-6, 4.2947227560776583326d-5, &
         6.3249861665073441312d-4, 7.9454472840953930811d-3, 7.0994327785661860575d-2, 4.4467219586283000332d-1, 1.7774588182255374745d0, &
         4.0003038986252717972d0, 3.9998233869142057195d0, 1.0000472932961288324d0 /
         DATA (a(i), i = 52, 64) / &
         2.7475684794982708655d-10, -3.8991472076521332023d-9, 1.9730170483976049388d-7, 5.9651531561967674521d-7, 5.1992971474748995357d-5, &
         5.7327338675433770752d-4, 8.2293143836530412024d-3, 6.9990934858728039037d-2, 4.4726764292723985087d-1, 1.7726685170014087784d0, &
         4.0062907863712704432d0, 3.9952750700487845355d0, 1.0016354346654179322d0 /

         DATA (b(i), i = 0, 13) / &
         6.7852367144945531383d-8, 4.6266061382821826854d-7, 6.9703135812354071774d-6, 7.6637663462953234134d-5, 7.9113515222612691636d-4, &
         7.3401204731103808981d-3, 6.0677114958668837046d-2, 4.3994941411651569622d-1, 2.7420017097661750609d0, 14.289661921740860534d0, &
         59.820609640320710779d0, 188.78998681199150629d0, 399.87313678256011180d0, 427.56411572180478514d0 /
         DATA (b(i), i = 14, 27) / &
         1.8042097874891098754d-7, 1.2277164312044637357d-6, 1.8484393221474274861d-5, 2.0293995900091309208d-4, 2.0918539850246207459d-3, &
         1.9375315654033949297d-2, 1.5985869016767185908d-1, 1.1565260527420641724d0, 7.1896341224206072113d0, 37.354773811947484532d0, &
         155.80993164266268457d0, 489.52113711585409180d0, 1030.9147225169564806d0, 1093.5883545113746958d0 /
         DATA (b(i), i = 28, 41) / &
         4.8017305613187493564d-7, 3.2613178439123800740d-6, 4.9073137508166159639d-5, 5.3806506676487583755d-4, 5.5387918291051866561d-3, &
         5.1223717488786549025d-2, 4.2190298621367914765d-1, 3.0463625987357355872d0, 18.895299447327733204d0, 97.915189029455461554d0, &
         407.13940115493494659d0, 1274.3088990480582632d0, 2670.9883037012547506d0, 2815.7166284662544712d0 /
         DATA (b(i), i = 42, 55) / &
         1.2789926338424623394d-6, 8.6718263067604918916d-6, 1.3041508821299929489d-4, 1.4282247373727478920d-3, 1.4684070635768789378d-2, &
         1.3561403190404185755d-1, 1.1152592585977393953d0, 8.0387088559465389038d0, 49.761318895895479206d0, 257.26842323135291380d0, &
         1066.8543146269566231d0, 3328.3874581009636362d0, 6948.8586598121634874d0, 7288.4893398212481055d0 /
         DATA (b(i), i = 56, 69) / &
         3.4093503681970328930d-6, 2.3079025203103376076d-5, 3.4691373283901830239d-4, 3.7949949772229085450d-3, 3.8974209677945602145d-2, &
         3.5949483804148783710d-1, 2.9522878893539528226d0, 21.246564609514287056d0, 131.28727387146173141d0, 677.38107093296675421d0, &
         2802.3724744545046518d0, 8718.5731420798254081d0, 18141.348781638832286d0, 18948.925349296308859d0 /

         DATA (c(i), i = 0, 8) / &
         2.5568678676452702768d-15, 3.0393953792305924324d-14, 6.3343751991094840009d-13, 1.5041298011833009649d-11, 4.4569436918556541414d-10, &
         1.7463930514271679510d-8, 1.0059224011079852317d-6, 1.0729838945088577089d-4, 5.1503226936425277380d-2 /
         DATA (c(i), i = 9, 17) / &
         5.2527963991711562216d-15, 7.2021184814210056410d-15, 7.2561421229904797156d-13, 1.4823121466731042510d-11, 4.4602670450376245434d-10, &
         1.7463600061788679671d-8, 1.0059226091322347560d-6, 1.0729838937545111487d-4, 5.1503226936437300716d-2 /
         DATA (c(i), i = 18, 26) / &
         1.3365917359358069908d-14, -1.2932643065888544835d-13, 1.7450199447905602915d-12, 1.0419051209056979788d-11, 4.5804788198059832600d-10, &
         1.7442405450073548966d-8, 1.0059461453281292278d-6, 1.0729837434500161228d-4, 5.1503226940658446941d-2 /
         DATA (c(i), i = 27, 35) / &
         5.3771611477352308649d-14, -1.1396193006413731702d-12, 1.2858641335221653409d-11, -5.9802086004570057703d-11, 7.3666894305929510222d-10, &
         1.6731837150730356448d-8, 1.0070831435812128922d-6, 1.0729733111203704813d-4, 5.1503227360726294675d-2 /
         DATA (c(i), i = 36, 44) / &
         3.7819492084858931093d-14, -4.8600496888588034879d-13, 1.6898350504817224909d-12, 4.5884624327524255865d-11, 1.2521615963377513729d-10, &
         1.8959658437754727957d-8, 1.0020716710561353622d-6, 1.0730371198569275590d-4, 5.1503223833002307750d-2 /

         w = x
         IF(w .LT. 8.5d0) THEN
            t = w*w*0.0625d0
            k = 13*int(t)
            y = (((((((((((a(k)*t + a(k + 1))*t + a(k + 2))*t + a(k + 3))*t + a(k + 4))*t + a(k + 5))*t + a(k + 6))*t + a(k + 7))*t + a(k + 8))*t &
                + a(k + 9))*t + a(k + 10))*t + a(k + 11))*t + a(k + 12)
            y = y*dexp(-w)
         ELSE IF(w .LT. 12.5d0) THEN
            k = int(w)
            t = w - k
            k = 14*(k - 8)
            y = ((((((((((((b(k)*t + b(k + 1))*t + b(k + 2))*t + b(k + 3))*t + b(k + 4))*t + b(k + 5))*t + b(k + 6))*t + b(k + 7))*t + b(k + 8))*t &
                + b(k + 9))*t + b(k + 10))*t + b(k + 11))*t + b(k + 12))*t + b(k + 13)
            y = y*dexp(-w)
         ELSE IF(w .LT. 25.d0) THEN
            t = 60/w
            k = 9*int(t)
            y = ((((((((c(k)*t + c(k + 1))*t + c(k + 2))*t + c(k + 3))*t + c(k + 4))*t + c(k + 5))*t + c(k + 6))*t + c(k + 7))*t &
                + c(k + 8))*dsqrt(t)*dexp(w)
            y = y*dexp(-w)
         ELSE
            y = 1.d0/dsqrt(2.d0*dacos(-1.d0)*w)*(1.d0 + 1.d0/(8.d0*w)*(1.d0 + 9.d0/(2.d0*8.d0*w)*(1.d0 + 25.d0/(6.d0*8.d0*w))))
         END IF

         DG0 = y

         END SUBROUTINE BESSEL_G0


!=======================================================================
!    SUBROUTINE BESSEL_G1
!
!         (1)calculate the 1th-order modified bessel function \Gamma_1s(b_s) = I_1(b_s)*exp(-b_s)
!=======================================================================
         SUBROUTINE BESSEL_G1(x, DG1)

         IMPLICIT NONE
         REAL*8    ::x, DG1
         REAL*8    ::delta, x1, x2, x3, x4
         REAL*8    ::DG0_1, DG0_2, DG0_3, DG0_4
         REAL*8    ::I0_1, I0_2, I0_3, I0_4
         REAL*8    ::I1

         delta = 0.000001d0
         x1    = x - 2.d0*delta
         x2    = x - 1.d0*delta
         x3    = x + 1.d0*delta
         x4    = x + 2.d0*delta

         CALL BESSEL_G0(x1, DG0_1)
         CALL BESSEL_G0(x2, DG0_2)
         CALL BESSEL_G0(x3, DG0_3)
         CALL BESSEL_G0(x4, DG0_4)

         I0_1  = DG0_1*dexp(x1)
         I0_2  = DG0_2*dexp(x2)     
         I0_3  = DG0_3*dexp(x3)
         I0_4  = DG0_4*dexp(x4)

         I1    = (- I0_4 + 8.d0*I0_3 - 8.d0*I0_2 + I0_1)/(12.d0*delta)
         DG1   = I1*dexp(-x)

         END SUBROUTINE BESSEL_G1

   

!=======================================================================
!    SUBROUTINE BESSEL_GN
!
!         (1)calculate the nth-order modified bessel function \Gamma_ns(b_s) = I_n(b_s)*exp(-b_s)
!=======================================================================
    real(wp) function bessel_jn_integral_real(x)
        implicit none
        real(wp),intent(in)::x
        bessel_jn_integral_real=real(exp((0,1.0)*(z_integral*sin(x)-n_integral*x)))
    end function bessel_jn_integral_real
    
    real(wp) function bessel_jn_integral_imag(x)
        implicit none
        real(wp),intent(in)::x
        bessel_jn_integral_imag=aimag(exp((0,1.0)*(z_integral*sin(x)-n_integral*x)))
    end function bessel_jn_integral_imag
     
    
    complex(wp) function bessel_jn_complex(n,z)
        implicit none
        integer,intent(in)::n
        complex(wp),intent(in)::z
        real(wp):: a 
        real(wp):: b 
        real(wp)::epsabs
        real(wp)::epsrel
        integer,parameter:: key = 6
        integer,parameter:: limit = 10000
        integer,parameter:: lenw=limit*4
        real(wp) :: abserr, ans_real,ans_imag, work(lenw)
        integer :: ier, iwork(limit), last, neval
        if (abs(z)>0.1) then 
            a=0.0_wp
            b=2*acos(-1.0_wp)
            epsabs=0
            epsrel=1d-8
            z_integral=z
            n_integral=n
            call dqag(bessel_jn_integral_real, a, b, epsabs, epsrel, key, ans_real, &
                        abserr, neval, ier, limit, lenw, last, &
                        iwork, work)
            call dqag(bessel_jn_integral_imag, a, b, epsabs, epsrel, key, ans_imag, &
                        abserr, neval, ier, limit, lenw, last, &
                        iwork, work)  
            bessel_jn_complex=ans_real/b*(1.0_wp,0)+ans_imag/b*(0,1.0_wp)
        else
            bessel_jn_complex=bessel_jn_complex_2(n,z)
            
        end if
      
    end function bessel_jn_complex
    
    complex(wp) function bessel_jn_complex_2(n,x)
        implicit none
        integer,intent(in)::n
        complex(wp),intent(in)::x
        integer::k,j
        real(wp)::res_1,res_2
        real(wp)::ln_factorial_1,ln_factorial_2
        complex(wp)::ln_sum,series_sum_1,series_sum_2
        series_sum_1=0.0_wp
        series_sum_2=0.0_wp
        
        do k=0,1000
            !call FACTORIAL(k, res_1)
            !call FACTORIAL(k+n,res_2)
            ln_factorial_1=0.0_wp
            ln_factorial_2=0.0_wp
            do j=1,k
                ln_factorial_1=ln_factorial_1+log(j*1.0_wp)
            end do
            
            do j=1,k+n
                ln_factorial_2=ln_factorial_2+log(j*1.0_wp)
            end do
            ln_sum=(2*k+n)*log(x/2)-ln_factorial_1-ln_factorial_2
            !series_sum_1=series_sum_1+(-1)**(k)*(x/2)**(2*k+n)/res_1/res_2
            series_sum_1=series_sum_1+(-1)**(k)*exp(ln_sum)
            if(k>=10 .and. abs(series_sum_1-series_sum_2)<1e-8*abs(series_sum_1+series_sum_2) )then
				exit
            end if
            series_sum_2=series_sum_1
        end do
        bessel_jn_complex_2=series_sum_2
    end function bessel_jn_complex_2
    
    complex(wp) function bessel_gn_complex(n,z)
        implicit none
        integer,intent(in)::n
        complex(wp),intent(in)::z
        if (real(z)>700) then
            bessel_gn_complex=(0.0_wp,0.0_wp)
        else
            bessel_gn_complex=exp(-(0,1.0)*n*acos(-1.0_wp)/2.0)*bessel_jn_complex(n,(0,1.0)*z)*exp(-z)
        end if
    end function bessel_gn_complex
      
    
    real(wp) function bessel_integral_real(x)
        implicit none
        real(wp),intent(in)::x
        bessel_integral_real=real(2*x*bessel_jn_complex(n_integral,(2*z1_integral)**0.5*x)**2*exp(-x**2))
    
    end function bessel_integral_real
    
    real(wp) function bessel_integral_imag(x)
        implicit none
        real(wp),intent(in)::x
        bessel_integral_imag=aimag(2*x*bessel_jn_complex(n_integral,(2*z1_integral)**0.5*x)**2*exp(-x**2))
    
    end function bessel_integral_imag
     
    subroutine bessel_gn_complex_2(z,n,gn)
        implicit none
        complex(wp),intent(in)::z
        integer,intent(in)::n
        complex(wp),intent(out)::gn
        real(wp):: a 
        real(wp):: b 
        real(wp)::epsabs
        real(wp)::epsrel
        integer,parameter:: key = 6
        integer,parameter:: limit = 10000
        integer,parameter:: lenw=limit*4
        real(wp) :: abserr, ans_real,ans_imag, work(lenw)
        integer :: ier, iwork(limit), last, neval

        a=0.0_wp
        b=10.0_wp
        epsabs=1d-7
        epsrel=1d-7
        z1_integral=z
        n_integral=n
        call dqag(bessel_integral_real, a, b, epsabs, epsrel, key, ans_real, &
                    abserr, neval, ier, limit, lenw, last, &
                    iwork, work)
        call dqag(bessel_integral_imag, a, b, epsabs, epsrel, key, ans_imag, &
            abserr, neval, ier, limit, lenw, last, &
            iwork, work)
        gn=cmplx(ans_real,ans_imag)
    end subroutine bessel_gn_complex_2
    
    real(wp) function bessel_integral(x)
        implicit none
        real(wp),intent(in)::x
        bessel_integral=2*x*bessel_jn(n_integral,(2*x_integral)**0.5*x)**2*exp(-x**2)
    
    end function bessel_integral
    
    subroutine bessel_gn(x,n,gn)
        implicit none
        real(wp),intent(in)::x
        integer,intent(in)::n
        real(wp),intent(out)::gn
        real(wp):: a 
        real(wp):: b 
        real(wp)::epsabs
        real(wp)::epsrel
        integer,parameter:: key = 6
        integer,parameter:: limit = 10000
        integer,parameter:: lenw=limit*4
        real(wp) :: abserr, ans, work(lenw)
        integer :: ier, iwork(limit), last, neval
        if (x>=700.0_wp) then
            gn=0.0_wp
        else
            a=0.0_wp
            b=10.0_wp
            epsabs=1d-7
            epsrel=1d-7
            x_integral=x
            n_integral=n
            call dqag(bessel_integral, a, b, epsabs, epsrel, key, ans, &
                        abserr, neval, ier, limit, lenw, last, &
                        iwork, work)
                
            gn=ans
        end if
    end subroutine bessel_gn
    
    
!=======================================================================
!    SUBROUTINE BESSEL_FLR
!
!         (1)calculate \Gamma_0s(b_s) - \Gamma_1s(b_s)
!         (2)calculate D(\Gamma_0s(b_s) - \Gamma_1s(b_s))/D(b_s)
!         (2)calculate D^2(\Gamma_0s(b_s) - \Gamma_1s(b_s))/D(b_s)^2
!=======================================================================
         SUBROUTINE BESSEL_GFLR(x, DG_FLR, DG_FLRP, DG_FLRPP)

         IMPLICIT NONE
         REAL*8    ::x, DG_FLR, DG_FLRP, DG_FLRPP
         REAL*8    ::delta, x1, x2, x3, x4
         REAL*8    ::DG0, DG0_1, DG0_2, DG0_3, DG0_4
         REAL*8    ::DG1, DG1_1, DG1_2, DG1_3, DG1_4
         REAL*8    ::DG_FLR_1, DG_FLR_2, DG_FLR_3, DG_FLR_4

         delta    = 0.01d0
         x1       = x - 2.d0*delta
         x2       = x - 1.d0*delta
         x3       = x + 1.d0*delta
         x4       = x + 2.d0*delta

         !k
         CALL BESSEL_G0(x, DG0)
         CALL BESSEL_G1(x, DG1)
         DG_FLR   = DG0 - DG1

         !k-2
         CALL BESSEL_G0(x1, DG0_1)
         CALL BESSEL_G1(x1, DG1_1)
         DG_FLR_1 = DG0_1 - DG1_1

         !k-1
         CALL BESSEL_G0(x2, DG0_2)
         CALL BESSEL_G1(x2, DG1_2)
         DG_FLR_2 = DG0_2 - DG1_2

         !k+1
         CALL BESSEL_G0(x3, DG0_3)
         CALL BESSEL_G1(x3, DG1_3)
         DG_FLR_3 = DG0_3 - DG1_3

         !k+2
         CALL BESSEL_G1(x4, DG1_4)
         CALL BESSEL_G0(x4, DG0_4)
         DG_FLR_4 = DG0_4 - DG1_4

         !calculate first derivative
         DG_FLRP  = (- DG_FLR_4 + 8.d0*DG_FLR_3 - 8.d0*DG_FLR_2 + DG_FLR_1)/(12.d0*delta)

         !calculate second derivative
         DG_FLRPP = (- DG_FLR_4 + 16.d0*DG_FLR_3 - 30.d0*DG_FLR + 16.d0*DG_FLR_2 - DG_FLR_1)/(12.d0*delta*delta)

         END SUBROUTINE BESSEL_GFLR


!=======================================================================
!    SUBROUTINE PDF
!
!=======================================================================
         SUBROUTINE PDF(omega, z)

         USE ERRFUN

         IMPLICIT NONE
         COMPLEX*16  ::omega, z
         COMPLEX*16  ::cj
         REAL*8      ::pi

         cj = (0.d0, 1.d0)
         pi = dacos(-1.d0)

         z  = cj*SQRT(pi)*wpop(omega)

         END SUBROUTINE PDF
         
         
!=======================================================================
!    function det:calcalate the determinant of a matrix using LU decompostion
!
!=======================================================================
      complex(wp) function det(a,n)
        complex(wp),intent(in)::a(:,:)
        integer,intent(in)::n
        complex(wp),allocatable::b(:,:)
        integer,allocatable::c(:)
        integer::info,i
        allocate(b(n,n))
        allocate(c(n))
        b=a
        call getrf(b,c,info)
        det=1
        do i=1,n
            det=det*b(i,i)
        end do
        do i=1,n
            if(c(i)/=i) then
                det=-det
            end if
        end do
        deallocate(b)
        deallocate(c)
      end function det
      


END MODULE MATHFUNCTION


