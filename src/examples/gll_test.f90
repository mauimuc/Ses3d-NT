!> @example gll_test.f90
!!
!! Example program which tests procedures of module gll_mod.
!!
!! @todo  Describe what it does ...
!!
!! Compile simply using the Makefile:
!!
!!     make gll_test
!!
!! Executing <tt>./gll_test</tt> runs all the tests.
PROGRAM gll_test 
    USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY : OUTPUT_UNIT, ERROR_UNIT
    USE gll_mod, ONLY : get_knots, get_weights

    IMPLICIT NONE 

    INTEGER :: i

    ! Check retuned array sizes
    if ( ANY( [ (SIZE(get_knots(i)), i=-2, 10 ) ] /= &
              [ [0,0], [ (i+1, i=0, 10 ) ] ] ) ) THEN 
        STOP '[FAIL] ...' 
    ELSE
        print *, '[PASS] ...'
    END IF
    
    if ( ANY( [ (SIZE(get_weights(i)), i=-2, 10 ) ] /= &
              [ [0,0], [ (i+1, i=0, 10 ) ] ] ) ) THEN 
        STOP '[FAIL] ...'
    ELSE
        print *, '[PASS] ...'
    END IF


    ! Check if first array entity is -1

    ! Check if last array entity is 1

    ! For even N check if the midpoint is 0

    ! Check if values match splib

CONTAINS

    
END PROGRAM 


!> @file
!! $Date: 2013-10-30 16:25:25 +0100 (Wed, 30 Oct 2013) $
!! $Author: mauerberger $
!! $Revision: 743 $
!! @copyright GNU General Public License version 3 or later
