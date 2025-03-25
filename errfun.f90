!$Id: errfun.f90 10 2016-08-15 19:43:56Z mexas $
! Copyright (c) 2016 Anton Shterenlikht, The University of Bristol, UK
!
! Module with error function and related functions.
! Two algorithms are implemented: Poppe and Wijers (wpop, erf_pop)
! and Zaghloul and Ali (wzag, erf_zag). The second algorithm is
! supposed to be more accurate for some values. However, the first
! algorithm can be over an order of magnitude faster.
!
! This code was written before the author became aware of
! M. R. Zaghloul, Remark on "Algorithm 916: Computing the
! Faddeyeva and Voight functions": efficiency improvements and
! Fortran translation, ACM Trans. Math. Software 42, Article 26, 2016.

module ERRFUN
use iso_fortran_env, only: real64
implicit none

integer, parameter :: rkind = real64
real( kind=rkind ), parameter :: zero = 0.0_rkind, one = 1.0_rkind,    &
   two = 2.0_rkind, half = 0.5_rkind, rmin = tiny( one ),              &
   eps0 = epsilon( one ), sqrt_log_rmin = sqrt( -log( rmin ) ),        &
   pi = 3.14159265358979323846264338327950288419_rkind,                &
   pi2 = pi * pi, one_sqrt_pi = one / sqrt( pi )
complex( kind=rkind ), parameter ::                                    &
   cmplxj = cmplx( zero, one, kind=rkind ),                            &
   cmplx0 = cmplx( zero, zero, kind=rkind )

private
public :: rkind, eps0, wpop, wzag, erf_pop, erf_zag

contains

!*********************************************************************72
elemental complex( kind=rkind ) function wzag( z, err )

! This is the Faddeyeva or the plasma dispersion function,
! w(z) = exp(-z**2) * erfc(-i*z). erfc(z) is the complex complementary
! error function of z. z is a complex number.
!
! Adapted from the Matlab code implementing TOMS
! Algorithm 916: http://www.netlib.org/toms/
!
! file: 916.zip
! ref: TOMS 38,2 (Dec 2011) Article: 15
! for: Computing the Faddeyeva and Voigt Functions
! by: Mofreh R. Zaghloul and Ahmed N. Ali
!
! Most of the code is calculation of equations (13)-(19) of the
! above paper.
!
! Inputs:
!    z - the argument of the function 
!  err - The desired accuracy (positive). err must be (err .le. 1.0e-4)
!        and  (err .ge. 0.06447*epsilon). For efficiency,
!        no checks are made!
!        The lowest accuracy and the fastest calculation are obtained
!        with err .eq. 1.0e-4. For higher accuracy use smaller err.

 complex( kind=rkind ), intent(in) :: z
 real( kind=rkind ), intent(in) :: err 

 real( kind=rkind ) :: a, half_a, a_sqr, two_a, two_a_sqr, &
   four_a_sqr, myerr, a_pi, two_a_pi, x, y, erfcsy, &
   xsign, ysign, x_sqr, y_sqr, two_yx, two_a_x, exp_x_sqr, &
   cos_2yx, sin_2yx, v_old, l_old, sigma1, sigma2, sigma3, &
   sigma4, sigma5, sigma4_5, delta3, delta5, exp1, exp2, exp3, &
   del2_tmp, del3_tmp, del3_3_tmp, den1, exp_del1, &
   exp3_den, exp3_3_den, del5, del3, two_exp_x_sqr_ysqr, aux13
 integer :: n, n3, n3_3

 x = real( z )
 y = aimag( z )

! For purely imaginary z, use intrinsic scaled complement of
! the error function, erfc_scaled (F2008 and beyond).
! Return immediately.
 if ( abs( x ) .eq. zero ) then
   wzag = erfc_scaled( y )
   return
 end if

      myerr = max( err, eps0 )
          a = sqrt( - pi2 / log( err / 2.0_rkind ) )
     half_a = half * a
      a_sqr = a**2
      two_a = 2 * a
  two_a_sqr = 2 * a_sqr
 four_a_sqr = 4 * a_sqr
       a_pi = a / pi
   two_a_pi = 2 * a_pi
     erfcsy = erfc_scaled( abs( y ) )
      xsign = sign( one, x )
      ysign = sign( one, y )
          x = abs( x )
          y = max( rmin, abs( y ) )
      x_sqr = x**2
      y_sqr = y**2
     two_yx = 2 * y * x
    two_a_x = two_a * x
  exp_x_sqr = exp( -x_sqr )
    cos_2yx = cos( two_yx )
    sin_2yx = sin( two_yx )
      v_old = exp_x_sqr * &
              ( erfcsy * cos_2yx + two_a_pi * sin( two_yx / 2)**2 / y )
      l_old = - erfcsy + a_pi / y
     sigma3 = rmin
     sigma5 = rmin
   sigma4_5 = zero
     delta3 = one
     delta5 = one
          n = 0
         n3 = ceiling( x / a )
       n3_3 = n3 - 1

outer: if ( (sqrt_log_rmin - x) .gt. 0 ) then
       sigma1 = rmin
       sigma2 = rmin
       sigma4 = rmin
         exp1 = exp( - two_a_x )
         exp2 = exp( four_a_sqr * n3 - 2 * two_a_x - two_a_sqr )
         exp3 = exp( - ( two_a_sqr * n3 - two_a_x - two_a_sqr ) )
     del2_tmp = one
     del3_tmp = exp( - ( a_sqr * n3**2 - two_a_x * n3                  &
                         - two_a_sqr * n3 + x_sqr + two_a_x + a_sqr ) )
   del3_3_tmp = exp( a_sqr - ( two_a_sqr * n3 - two_a_x ) )

   loop1: do
     if ( delta3 .lt. myerr .and. delta5 .lt. myerr .and.              &
                                        n .gt. 50 ) exit loop1
            n = n + 1
         den1 = a_sqr * n**2 + y_sqr
     exp_del1 = exp( -( a_sqr * n**2) ) / den1
     del2_tmp = del2_tmp * exp1

     minor: if ( n3_3 .ge. 1 ) then
         del3_tmp = del3_tmp * exp3
         exp3_den = del3_tmp * exp_del1 *                              &
                       ( den1 / ( a_sqr * n3**2 + y_sqr) )
       del3_3_tmp = del3_3_tmp * exp2
       exp3_3_den = exp3_den * del3_3_tmp *                            &
                     ( (a_sqr * n3**2 + y_sqr ) /                      &
                       (a_sqr * n3_3**2 + y_sqr ) )
             del5 = n3_3 * exp3_3_den + n3 * exp3_den
             del3 = exp3_3_den + exp3_den
     else
         del3_tmp = del3_tmp * exp3
             del3 = del3_tmp * exp_del1 *                              &
                           ( den1 / ( a_sqr * n3**2 + y_sqr ) )
             del5 = n3 * del3
     end if minor

     delta3 = del3 / sigma3
     delta5 = del5 / sigma5
     sigma1 = sigma1 + exp_del1
     sigma2 = sigma2 + del2_tmp * exp_x_sqr * exp_del1
     sigma3 = sigma3 + del3
     sigma4 = sigma4 + n * del2_tmp * exp_x_sqr * exp_del1
     sigma5 = sigma5 + del5

     if ( x .ge. 5.0e-4_rkind ) then
       sigma4_5 = - sigma4 + sigma5
     else
       sigma4_5 = sigma4_5 + 2 * n**2 * two_a_x * exp_x_sqr            &
                   * exp_del1 * ( one + 1.666666666666667e-1_rkind     &
                   * ( two_a_x * n )**2 + 8.333333333333333e-3_rkind   &
                   * ( two_a_x * n )**4 )
     end if

       n3 = n3 + 1
     n3_3 = n3_3 - 1

   end do loop1

   ! Second line of Eqn (13)
   aux13 = y * two_a_pi *                                              &
     ( - cos_2yx*exp_x_sqr*sigma1 + half*(sigma2 + sigma3) )
   mumu: if ( y .le. 5.0_rkind .and. two_yx .gt. rmin ) then
     wzag = v_old + aux13 + cmplxj * xsign             &
           * ( sin_2yx * exp_x_sqr * ( l_old + two_a_pi * y * sigma1 ) &
           + two_a_pi * half_a * sigma4_5 )
   else if ( y .le. 5.0_rkind .and. two_yx .le. rmin ) then
     wzag = v_old + aux13 + cmplxj * xsign             &
           * ( two * y * exp_x_sqr * ( x * l_old + x * two_a_pi * y    &
           * sigma1 ) + two_a_pi * half_a * sigma4_5 )
   else
     wzag = v_old + aux13 + cmplxj * xsign             &
           * ( sin_2yx * exp_x_sqr * min( zero, abs( l_old             &
           + ( two_a_pi * y * sigma1) ) ) + two_a_pi *half_a           &
           * sigma4_5 )
   end if mumu

 else if ( x .ge. sqrt_log_rmin .and. x .lt. 1.0e15_rkind ) then

         exp2 = exp( four_a_sqr * n3 - 2 * two_a_x - two_a_sqr )
   del3_3_tmp = exp( a_sqr + two_a_x - two_a_sqr * n3 )

   loop2: do
     if ( delta3 .lt. myerr .and. delta5 .lt. myerr .and.              &
                                        n .gt. 50 ) exit loop2
     n = n + 1
     if ( n3_3 .ge. 1 ) then
         exp3_den = exp( - ( a * n3 - x ) * ( a * n3 - x ) )           &
                     / ( a_sqr * n3**2 + y_sqr )
       del3_3_tmp = del3_3_tmp * exp2
       exp3_3_den = exp3_den * del3_3_tmp * ( ( a_sqr * n3**2 + y_sqr) &
                     / ( a_sqr * n3_3**2 + y_sqr ) )
             del5 = n3_3 * exp3_3_den + n3 * exp3_den
             del3 = exp3_3_den + exp3_den
     else
             del3 = exp( - (a * n3 - x)**2) / ( a_sqr * n3**2 + y_sqr )
             del5 = n3 * del3
     end if

     delta3 = del3 / sigma3
     delta5 = del5 / sigma5
     sigma3 = sigma3 + del3
     sigma5 = sigma5 + del5
         n3 = n3 + 1
       n3_3 = n3_3 - 1

   end do loop2

   wzag = v_old + y * a_pi * sigma3 + cmplxj * xsign * ( sin_2yx        &
         * exp_x_sqr * l_old + two_a_pi * half_a * sigma5 )

 else
   wzag = one_sqrt_pi * ( ( y + cmplxj * xsign * x ) / ( x_sqr + y_sqr ))
 end if outer

 if ( ysign .lt. zero ) then
   two_exp_x_sqr_ysqr = two * exp( - x_sqr + y_sqr )
   wzag = two_exp_x_sqr_ysqr * cos_2yx - real( wzag ) - cmplxj           &
         * ( - xsign * two_exp_x_sqr_ysqr * sin_2yx - aimag( wzag ) ) 
 end if

end function wzag

!*********************************************************************72
elemental complex( kind=rkind ) function erf_zag( z, err )

! This is an error function of a complex argument, which uses wzag(z).

 complex( kind=rkind ), intent(in) :: z
 real( kind=rkind ), intent(in) :: err 

 erf_zag = one - wzag( cmplxj * z, err ) * exp( - z**2 )

end function erf_zag

!*********************************************************************72
elemental complex( kind=rkind ) function wpop( z )

! A modified version of algorithm 680, rewritten in Fortran 2008.
! G.P.M. Poppe, C.M.J. Wijers, More efficient computation of
! the complex error-function, ACM Trans. Math. Software 16:38-46, 1990.
!  and
! G.P.M. Poppe, C.M.J. Wijers, Algorithm 680, Evaluation of the
! complex error function, ACM Trans. Math. Software 16:47, 1990.
!
! Given a complex number z, this function computes
! the value of the Faddeeva-function w(z) = exp(-z**2)*erfc(-i*z),
! where erfc is the complex complementary error-function and i
! means sqrt(-1).  The accuracy of the algorithm for z in the 1st
! and 2nd quadrant is 14 significant digits; in the 3rd and 4th
! it is 13 significant digits outside a circular region with radius
! 0.126 around a zero of the function.

 complex( kind=rkind ), intent(in) :: z

!  factor is 2/sqrt(pi)
 real( kind=rkind ), parameter ::                                      &
   factor = 1.12837916709551257388_rkind

 logical :: a, b

 real( kind=rkind ) :: xabs, yabs, x, y, qrho, xabsq, yquad, c, daux,  &
   h, h2, qlambda, rx, ry, sx, sy, tx, ty, u1, u2, v1, v2, w1, xaux,   &
   xquad, xsum, ysum, xi, yi, u, v

 integer :: i, j, kapn, n, np1, nu

! To avoid the complier uninitialised varning
         h2 = zero

         xi = real( z )
         yi = aimag( z )
       xabs = abs( xi )
       yabs = abs( yi )
          x = xabs / 6.3_rkind
          y = yabs / 4.4_rkind
       qrho = x**2 + y**2
      xabsq = xabs**2
      xquad = xabsq - yabs**2
      yquad = 2*xabs*yabs

          a = qrho .lt. 0.085264_rkind

  branch1: if ( a ) then

    ! If ( qrho .lt. 0.085264 ) then the Faddeeva-function is evaluated
    !  using a power-series (abramowitz/stegun, equation (7.1.5), p.297)
    !  n is the minimum number of terms needed to obtain the required
    !  accuracy

        qrho  = ( one - 0.85_rkind * y ) * sqrt( qrho )
            n = nint( 6.0_rkind + 72.0_rkind * qrho )
            j = 2 * n + 1
         xsum = one / real( j, kind=rkind )
         ysum = zero

        do i = n, 1, -1
          j = j - 2
          xaux = ( xsum * xquad - ysum * yquad) / real( i, kind=rkind )
          ysum = ( xsum * yquad + ysum * xquad) / real( i, kind=rkind )
          xsum = xaux + one / real( j, kind=rkind )
        end do

          u1 = - factor * ( xsum * yabs + ysum * xabs ) + one
          v1 = factor * ( xsum * xabs - ysum * yabs )
        daux = exp( -xquad )
          u2 = daux * cos( yquad )
          v2 = -daux * sin( yquad )
           u = u1 * u2 - v1 * v2
           v = u1 * v2 + v1 * u2
  else

    bran2: if ( qrho .gt. one ) then

      ! If ( qrho .gt. 1) then w(z) is evaluated using the laplace
      ! continued fraction. nu is the minimum number of terms needed
      ! to obtain the required accuracy.

           h = zero
        kapn = 0
        qrho = sqrt( qrho )
          nu = int( 3.0_rkind + (1442.0_rkind / ( 26.0_rkind * qrho &
                + 77.0_rkind )))

    else

      ! If ( qrho .ge. 0.085264 .and. qrho .le. one ) then
      ! w(z) is evaluated by a truncated Taylor expansion,
      ! where the Laplace continued fraction is used to calculate
      ! the derivatives of w(z). KAPN is the minimum number of terms
      ! in the Taylor expansion needed to obtain the required accuracy.
      ! NU is the minimum number of terms of the continued fraction
      ! needed to calculate the derivatives with the required accuracy.

          qrho = ( one - y ) * sqrt( one - qrho )
             h = 1.88_rkind * qrho
            h2 = two * h
          kapn = nint( 7.0_rkind + 34.0_rkind * qrho )
          nu   = nint( 16.0_rkind + 26.0_rkind * qrho )

    end if bran2

    b = h .gt. zero

! To avoid uninitialise compiler warning. qlambda is used
! only if (b), so can define to any value otherwise.
    qlambda = zero 
    if ( b ) qlambda = h2**kapn

    rx = zero
    ry = zero
    sx = zero
    sy = zero

    do n = nu, 0, -1
         np1 = n + 1
          tx = yabs + h + np1 * rx
          ty = xabs - np1 * ry
           c = half / (tx**2 + ty**2)
          rx = c * tx
          ry = c * ty
          if ( b .and. n .le. kapn ) then
            tx = qlambda + sx
            sx = rx*tx - ry*sy
            sy = ry*tx + rx*sy
            qlambda = qlambda / h2
          end if
    end do

        if ( h .eq. zero ) then
          u = factor * rx
          v = factor * ry
        else
          u = factor * sx
          v = factor * sy
        end if

        if ( yabs .eq. zero ) u = exp( -xabs**2 )

  end if branch1

! Evaluation of w(z) in the other quadrants

      if ( yi .lt. zero ) then

        if ( a ) then
          u2 = two * u2
          v2 = two * v2
        else
          xquad = -xquad
          w1 = two * exp( xquad )
          u2 = w1 * cos( yquad )
          v2 = -w1 * sin( yquad )
        end if

        u = u2 - u
        v = v2 - v
        if ( xi .gt. zero ) v = -v
      else
        if ( xi .lt. zero ) v = -v
      end if

 wpop = cmplx( u, v, kind=rkind )

end function wpop

!*********************************************************************72
elemental complex( kind=rkind ) function erf_pop( z )

! This is an error function of a complex argument, which uses wpop(z).

 complex( kind=rkind ), intent(in) :: z

 erf_pop = one - wpop( cmplxj * z ) * exp( - z**2 )

end function erf_pop

!*********************************************************************72
end module ERRFUN
