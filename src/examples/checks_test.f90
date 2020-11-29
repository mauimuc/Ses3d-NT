!> @example checks_test.f90
!!
!! Example program which tests procedures of the module checks_mod
!! The following checks are carried out:
!!  * Check if the source file exists
!!  * Check existence of multiple files
!!  * Check latitude
!!  * Check longitude
!!  * Check depth
!!
!! Compile with simply using the Makefile:
!!
!!      make checks_test
!!
!! Run the program by executing:
!!
!!      ./checks_test
PROGRAM checks_test
    USE parameters_mod, ONLY : OUTPUT_UNIT, ERROR_UNIT, &
                               real_kind, fnl ! File name length
    USE checks_mod, ONLY : file_exists, valid_lat, valid_lon, valid_depth

    IMPLICIT NONE

    INTEGER :: i
    REAL(real_kind) :: lat, lon, depth
    CHARACTER(LEN=*), PARAMETER :: file_name = 'checks_test.f90'
    CHARACTER(LEN=fnl), PARAMETER :: file_names(2) = [ 'checks_test.f90     ', &
                                                       'dummy_parameters.f90' ]

    ! Check if the source file exists
    WRITE(UNIT=OUTPUT_UNIT, FMT='(/,A)') "Testing logical function 'file_exists(file_name)' for single file:"
    IF ( file_exists(file_name) ) THEN
        WRITE(UNIT=OUTPUT_UNIT, FMT='(1x, 3A)') "File '", file_name, "' exists!"
    ELSE
        WRITE(UNIT=ERROR_UNIT, FMT='(1x, 3A)') "File ", file_name, " does NOT exist!"
        STOP 1
    END IF

    ! Check existence of multiple files
    WRITE(UNIT=OUTPUT_UNIT, FMT='(/,A)') "Testing logical function 'file_exists(file_name)' for multiple files:"
    IF ( ALL( file_exists(file_names) ) ) THEN
        WRITE(UNIT=OUTPUT_UNIT, FMT=*) "Files ", &
           ("'"//TRIM(file_names(i))//"', ", i=1, SIZE(file_names)), &
            "exist!"
    ELSE
        WRITE(UNIT=ERROR_UNIT, FMT=*) "At least one File out of ", &
           ("'"//TRIM(file_names(i))//"', ", i=1, SIZE(file_names)), &
           "does NOT exist!"
        STOP 1
    END IF


    ! Check whether latitude (in degrees) is in the range [-90,90]
    WRITE(UNIT=OUTPUT_UNIT, FMT='(/,A)') "Testing logical function 'valid_lat(lat)':"
    lat = 90.0_real_kind
    IF ( valid_lat(lat)) THEN
        WRITE(UNIT=OUTPUT_UNIT, FMT=*) "Latitude '",lat,"' is in the range!"
    ELSE
        WRITE(UNIT=ERROR_UNIT, FMT=*) "Latitude '",lat,"' is out of range!"
        STOP 1
    END IF

    lat = -90.0_real_kind
    IF ( valid_lat(lat)) THEN
        WRITE(UNIT=OUTPUT_UNIT, FMT=*) "Latitude '",lat,"' is in the range!"
    ELSE
        WRITE(UNIT=ERROR_UNIT, FMT=*) "Latitude '",lat,"' is out of range!"
        STOP 1
    END IF  


    ! Check whether longitude (in degrees) is in the range [-180,180]
    WRITE(UNIT=OUTPUT_UNIT, FMT='(/,A)') "Testing logical function 'valid_lon(lon)':"
    lon = 180.0_real_kind
    IF ( valid_lon(lon)) THEN
        WRITE(UNIT=OUTPUT_UNIT, FMT=*) "Longitude '",lon,"' is in the range!"
    ELSE
        WRITE(UNIT=ERROR_UNIT, FMT=*) "Longitude '",lon,"' is out of range!"
        STOP 1
    END IF

    lon = -180.0_real_kind
    IF ( valid_lon(lon)) THEN
        WRITE(UNIT=OUTPUT_UNIT, FMT=*) "Longitude '",lon,"' is in the range!"
    ELSE
        WRITE(UNIT=ERROR_UNIT, FMT=*) "Longitude '",lon,"' is out of range!"
        STOP 1
    END IF
  

    ! Check whether depth (in meters) is in the range [0,6371000]
    WRITE(UNIT=OUTPUT_UNIT, FMT='(/,A)') "Testing logical function 'valid_depth(depth)':"
    
    ! 1. Above the surface 
    depth = -TINY(0.0_real_kind) 
    IF ( valid_depth(depth) ) THEN
        STOP                           "[FAIL] Checking depth above the Surface!"
    ELSE
        WRITE(UNIT=OUTPUT_UNIT, FMT=*) "[PASS] Checking depth above the Surface!"
    END IF
    ! 2. At the surface
    depth = 0.0_real_kind
    IF ( valid_depth(depth) ) THEN
        WRITE(UNIT=OUTPUT_UNIT, FMT=*) "[PASS] Checking depth at the surface!"
    ELSE
        STOP                           "[FAIL] Checking depth at the surface!"
    END IF
    ! 2. Somewhere in between 
    depth = 1230.0_real_kind
    IF ( valid_depth(depth) ) THEN
        WRITE(UNIT=OUTPUT_UNIT, FMT=*) "[PASS] Checking depth inside the Earth!"
    ELSE
        STOP                           "[FAIL] Checking depth inside the Earth!"
    END IF
    ! 4. At the core
    depth = 6371000.0_real_kind
    IF ( valid_depth(depth) ) THEN
        STOP                           "[FAIL] Checking depth at the center!"
    ELSE
        WRITE(UNIT=OUTPUT_UNIT, FMT=*) "[PASS] Checking depth at the center!"
    END IF
    
END PROGRAM checks_test


!> @file
!! $Date: 2013-10-30 16:25:25 +0100 (Wed, 30 Oct 2013) $
!! $Author: mauerberger $
!! $Revision: 743 $
!! @copyright GNU General Public License version 3 or later
