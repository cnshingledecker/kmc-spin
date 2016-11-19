MODULE mc_toolbox

  !  Here one can specify parameters to be used  by the functions
  INTEGER(KIND=4), PARAMETER :: GSEED=12345
  REAL(KIND=8), PARAMETER :: MEANVAL=0.00
  REAL(KIND=8), PARAMETER :: SDEV=2.00

  CONTAINS
    FUNCTION rgamma(a)
      ! This function returns a random number with a Gamma PDF using
      ! the method described by Marsaglia and Tsang 2000.
      IMPLICIT NONE

      ! Data dictionary: Calling parameters
      REAL(KIND=8), INTENT(IN) :: a
      REAL(KIND=8)             :: rgamma

      ! Data dictionary: Local variables
      REAL(KIND=8)             :: d,c,v
      REAL(KIND=8)             :: x,u

      d = a-1./3.
      c = 1./SQRT(9.*d)
      rgamma = 0.
      DO
        v = 0.
        DO WHILE ( v .LE. 0. )
          x = R8_NORMAL_AB(MEANVAL,SDEV)
          v = 1. + c*x
        END DO
!        PRINT *, 'x=',x,'and v=',v
        v = v*v*v
!        u = RANDOM()
        CALL RANDOM_NUMBER(u)
        IF ( u .LT. 1.-0.331*(x*x)*(x*x) ) THEN
!          PRINT *, 'Success-1!'
          rgamma = (d*v)
          GOTO 10
        ELSE IF ( DLOG(u) .LT. 0.5*x*x+d*(1.-v+DLOG(v)) ) THEN
!          PRINT *, 'Success-2!'
          rgamma = (d*v)
          GOTO 10
        END IF
      END DO
10    RETURN
    END FUNCTION rgamma

FUNCTION r8_normal_01 ()
!*****************************************************************************80
!
!! R8_NORMAL_01 returns a unit pseudonormal R8.
!
!  Discussion:
!
!    The standard normal probability distribution function (PDF) has
!    mean 0 and standard deviation 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 August 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, real ( kind = 8 ) R8_NORMAL_01, a normally distributed
!    random value.
!
  implicit none

  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  real ( kind = 8 ) r8_normal_01
  real ( kind = 8 ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = 8 ) x

!  r1 = RANDOM()
!  r2 = RAND()
  ! Here we'll use a couple of random numbers generated using the built in
  ! RNG of the compiler and previously defined seed
  CALL RANDOM_NUMBER(r1)
  CALL RANDOM_NUMBER(r2)
  x = sqrt( - 2.0D+00 * DLOG( r1 ) ) * cos( 2.0D+00 * r8_pi * r2 )

  r8_normal_01 = x

  return
END FUNCTION r8_normal_01

FUNCTION r8_uniform_01 ( seedy )

!*****************************************************************************80
!
!! R8_UNIFORM_01 returns a unit pseudorandom R8.
!
!  Discussion:
!
!    This routine implements the recursion
!
!      seed = 16807 * seed mod ( 2^31 - 1 )
!      r8_uniform_01 = seed / ( 2^31 - 1 )
!
!    The integer arithmetic never requires more than 32 bits,
!    including a sign bit.
!
!    If the initial seed is 12345, then the first three computations are
!
!      Input     Output      R8_UNIFORM_01
!      SEED      SEED
!
!         12345   207482415  0.096616
!     207482415  1790989824  0.833995
!    1790989824  2035175616  0.947702
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Second Edition,
!    Springer, 1987,
!    ISBN: 0387964673,
!    LC: QA76.9.C65.B73.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, December 1986, pages 362-376.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley, 1998,
!    ISBN: 0471134031,
!    LC: T57.62.H37.
!
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, 1969, pages 136-143.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.
!    On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R8_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
!
  implicit none

  integer ( kind = 4 ) k
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seedy

  PRINT *, '*******/Now in R8_UNIFORM_01/*******'
  PRINT *, 'Seed is:',seedy
  k = seedy / 127773
  PRINT *, 'k=',k
  seedy = 16807 * ( seedy - k * 127773 ) - k * 2836
  PRINT *, 'seedy now =',seedy


  if ( seedy < 0 ) then
    seedy = seedy + 2147483647
  end if
!
!  Although SEED can be represented exactly as a 32 bit integer,
!  it generally cannot be represented exactly as a 32 bit real number!
!
  r8_uniform_01 = real ( seedy, kind = 8 ) * 4.656612875D-10
  PRINT *, 'The uniform RN is:',r8_uniform_01
  return
END FUNCTION r8_uniform_01

FUNCTION r8_normal_ab ( a, b)
!*****************************************************************************80
!
!! R8_NORMAL_AB returns a scaled pseudonormal R8.
!
!  Discussion:
!
!    The normal probability distribution function (PDF) is sampled,
!    with mean A and standard deviation B.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 August 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, the mean of the PDF.
!
!    Input, real ( kind = 8 ) B, the standard deviation of the PDF.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, real ( kind = 8 ) R8_NORMAL_AB, a sample of the normal PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  real ( kind = 8 ) r8_normal_ab
  real ( kind = 8 ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = 8 ) x

!  r1 = RANDOM()
  CALL RANDOM_NUMBER(r1)
  CALL RANDOM_NUMBER(r2)
!  r2 = RANDOM()
  x = sqrt ( - 2.0D+00 * log ( r1 ) ) * cos ( 2.0D+00 * r8_pi * r2 )

  r8_normal_ab = a + b * x

  return
END FUNCTION r8_normal_ab

END MODULE mc_toolbox
