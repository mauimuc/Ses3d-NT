!> @example coordinates_test.f90
!!
!! Example program which tests procedures of the module coordinates_mod. 
!! The following checks are carried out:
!!  * Consistency of Upper and lower coordinate bounds 
!!  * Simple coordinates check
!!  * Check if neighbouring elements are connecting seamlessly 
!!
!! @todo Vary c_min, c_max for each iteration
!! @todo Check midpoint for even lpd values 
!! @todo Check symmetry
!!
!! Compile by simply useing the Makefile `make coordinates_test` 
!! and execute via `./coordinates_test`. All tests should pass. 
PROGRAM coordinates_test
    USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY : OUTPUT_UNIT, ERROR_UNIT
    USE parameters_mod, ONLY : rk => real_kind
    USE coordinates_mod, ONLY : theta_coordinates, phi_coordinates, &
                                r_coordinates, pack_r2tor1

    IMPLICIT NONE

    INTEGER :: i, j, elms, lpd, scl
    LOGICAL :: l
    REAL(rk) :: c_min, c_max
    REAL(rk), ALLOCATABLE :: coords(:,:), coords_1d(:,:)
    CHARACTER(LEN=128) :: msg 

    ! Check if coordinates at lower and upper bounds are matching 
    c_min = -1.234136
    c_max = 14.237239
    DO i=1, 13, 3
        DO j=1, 7
            coords_1d = theta_coordinates(elms=i, lpd=j, min=c_min, max=c_max)
            WRITE(UNIT=msg, FMT='("Checking lower and upper coordinate bounds&
                & for elms=", I0, " and lpd=", I0)') i, j
            l = (coords_1d(1,1)==c_min .AND. coords_1d(i,j+1)==c_max)
            CALL check(l, msg)
        END DO
    END DO

    ! A very simple consistency check of element boundaries 
    DO elms=1, 50, 7
        c_min = 0.0
        c_max = REAL(elms,rk)
        DO j=1, 7
            coords = r_coordinates(elms=elms, lpd=j, min=c_min, max=c_max)
            coords = coords(:,[1,j+1])
            l = ALL( pack_r2tor1(coords) == REAL([(i, i=0, elms)],rk) )
            WRITE(UNIT=msg, FMT='("Simple coordinates check for elms=", I0, &
                &" lpd=",I0)') elms, j
            CALL check(l, msg)
        END DO
    END DO
    
    ! Check if element boundaries are connecting seamlessly
    DO elms=1, 235, 43
        DO lpd=1, 7
            coords = r_coordinates(elms=elms, lpd=lpd, min=-1.1_rk, max=1.3_rk)
            WRITE(UNIT=msg, FMT='("Element boundaries are connection &
                &seamlessly for elms=", I0, " lpd=",I0)') elms, lpd
            l = ALL( [(coords(i,lpd+1) == coords(i+1,1), i=1, elms-1) ] )
            CALL check(l, msg)
        END DO
    END DO

    !! Check symmetry

CONTAINS
    
    SUBROUTINE check(l,s)
        LOGICAL :: l
        CHARACTER(LEN=*) :: s

        IF ( l ) THEN
            WRITE(UNIT=OUTPUT_UNIT, FMT='(A,A)') '[PASS] ', TRIM(s)
        ELSE
            WRITE(UNIT=ERROR_UNIT, FMT='(A,A)') '[FAIL] ', TRIM(s)
            STOP 1
        END IF

    END SUBROUTINE
    
END PROGRAM coordinates_test


!> @file
!! $Date: 2013-11-01 21:20:02 +0100 (Fri, 01 Nov 2013) $
!! $Author: mauerberger $
!! $Revision: 756 $
!! @copyright GNU General Public License version 3 or later
