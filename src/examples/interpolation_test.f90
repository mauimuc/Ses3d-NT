!> @example interpolation_test.f90
!!
!! Example program which tests procedures for Lagrange interpolation of module
!! interpolation_mod.
!! * Interpolation of constant function
!! * Interpolation of linear function
!! * Interpolation of sin(x) function
!!
!! Compile simply using the Makefile `make interpolation_test` and run the
!! program by executing `./interpolation_test`. All tests have to pass through.
PROGRAM interpolation_test
    USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY : OUTPUT_UNIT, ERROR_UNIT
    USE parameters_mod, ONLY : rk => real_kind, pi
    USE interp_mod

    IMPLICIT NONE

    INTEGER :: i, n
    LOGICAL :: l
    REAL(rk), ALLOCATABLE :: knots(:), weights(:), f(:), x(:), f_(:)
    REAL(rk) :: e


    ! Interpolate constant function
    WRITE(UNIT=OUTPUT_UNIT, FMT='(/, "Interpolation of constant function ", /)')

    knots = [-1.0,0.0,1.0]
    f = const(knots)
    weights = calc_weights(knots)

    n = 50
    x = [(REAL(i, rk), i=-n, n)]
    n = SIZE(x)
    ! Add random perturbations to x
    x(2:n-1) = x(2:n-1) + random(n-2)
    x(:) = x/MAXVAL(ABS(x))

    f_ = [(SUM(calc_lagrange_polynomial(x(i),knots,weights)*f), i=1, SIZE(x))]

    l = ALL( f_(:)/const(x(:))-1.0 < 0.01 )

    ! Average interpolation error
    e = SUM(ABS(f_(:)/const(x(:))-1.0))/n

    CALL check( l, "Constant; error < 0.1%")



    ! Interpolate linear function
    WRITE(UNIT=OUTPUT_UNIT, FMT='(/, "Interpolation of linear function ", /)')

    knots = [-1.0, 1.0]
    f = linear(knots)
    weights = calc_weights(knots)

    n = 3
    x = [ (REAL(i,rk), i=-n, n) ]
    n = SIZE(x)
    ! Add random perturbations to x
    x(2:n-1) = x(2:n-1) + random(n-2)
    x(:) = x/MAXVAL(ABS(x))

    f_ = [(SUM(calc_lagrange_polynomial(x(i),knots,weights)*f), i=1, SIZE(x))]

    l = ALL( f_(:)/const(x(:))-1.0 < 0.01 )

    ! Average interpolation error
    e = SUM(ABS(f_(:)/const(x(:))-1.0))/n

    CALL check( l, "linear; error < 0.1%")



    ! Interpolate sin
    WRITE(UNIT=OUTPUT_UNIT, FMT='(/, "Interpolation of sin(x)", /)')

    n = 10
    knots = REAL([(i, i=-n/2, n/2)], rk ) / (n/2) * pi
    f = [SIN(knots)]
    weights = calc_weights(knots)

    n = 24
    x = [(REAL(i, rk), i=-n, n)]
    n = SIZE(x)
    ! Add random perturbations to x
    x(2:n-1) = x(2:n-1) + random(n-2)
    x(:) = x/MAXVAL(ABS(x)) * pi

    f_ = [(SUM(calc_lagrange_polynomial(x(i),knots,weights)*f), i=1, SIZE(x))]

    l = ALL( f_(:)/const(x(:))-1.0 < 0.01 )

    ! Average interpolation error
    e = SUM(ABS(f_(:)/const(x(:))-1.0))/n

    CALL check( l, "sin(x); error < 0.1%")



CONTAINS

    REAL(rk) ELEMENTAL FUNCTION const(x)
        REAL(rk), INTENT(IN) :: x
        const = 1.0_rk
    END FUNCTION

    REAL(rk) ELEMENTAL FUNCTION linear(x)
        REAL(rk), INTENT(IN) :: x
        linear = x
    END FUNCTION

    SUBROUTINE check(l, s)
        LOGICAL :: l
        CHARACTER(LEN=*) :: s

        IF ( l ) THEN
            WRITE(UNIT=OUTPUT_UNIT, FMT='(A,A)') '[PASS] ', TRIM(s)
        ELSE
            WRITE(UNIT=ERROR_UNIT, FMT='(A,A)') '[FAIL] ', TRIM(s)
            STOP 1
        END IF

    END SUBROUTINE

    !> Returns an array of s random values in the range of [-0.5:0.5]
    FUNCTION random(s) RESULT(r)
        ! Passed dummy arguments
        implicit none
        INTEGER, INTENT(IN) :: s
        ! Return value
        REAL(rk) :: r(s)
        ! Local variables
        INTEGER :: lun, n
        INTEGER, ALLOCATABLE :: seed(:)

        ! Obtain size of random seed
        CALL RANDOM_SEED(SIZE=n)
        ALLOCATE(seed(n))
        ! Invoke the OS random number generator to get a random seed
        OPEN(NEWUNIT=lun, FILE="/dev/urandom", ACCESS="stream", &
             FORM="unformatted")!, ACTION="read", STATUS="old" )
        READ(UNIT=lun) seed
        CLOSE(UNIT=lun)

        ! Put a random seed
        CALL RANDOM_SEED(PUT=seed)

        ! Obtain random values for each element
        CALL RANDOM_NUMBER(r(:))

        ! Scale to range -0.5 to 0.5
        r = r - 0.5_rk

    END FUNCTION

END PROGRAM interpolation_test


!> @file
!! $Date: 2013-11-04 18:34:49 +0100 (Mon, 04 Nov 2013) $
!! $Author: mauerberger $
!! $Revision: 766 $
!! @copyright GNU General Public License version 3 or later
