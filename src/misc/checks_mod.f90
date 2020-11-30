!> @file
!! Ses3d-NT - simulation of elastic wave propagation in spherical sections
!!
!! (c) by Stefan Mauerberger
!!    and Maksym Melnyk
!!
!! This program is free software: you can redistribute it and/or modify
!! under the terms of the GNU General Public License as published by
!! the Free Software Foundation, either version 3 of the License, or
!! (at your option) any later version.
!!
!! This program is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!!
!! $Date: 2013-11-01 18:19:11 +0100 (Fri, 01 Nov 2013) $
!! $Author: mmelnyk $
!! $Revision: 754 $
!> @copyright GNU General Public License version 3 or later



!> This module checks whether the file(s) exist(s). Also checks whether
!! latitude, longitude and depth satisfy the corresponding ranges.

MODULE checks_mod

    IMPLICIT NONE

    PRIVATE

    !> Checks if files exist
    INTERFACE file_exists
        MODULE PROCEDURE file_exists_scalar
        MODULE PROCEDURE file_exists_array
    END INTERFACE

    PUBLIC :: file_exists, valid_lat, valid_lon, valid_depth

CONTAINS

    !> Checks multiple files existence
    FUNCTION file_exists_array(file_name)
        CHARACTER(LEN=*), INTENT(IN) :: file_name(:)
        LOGICAL :: file_exists_array(SIZE(file_name))
        INTEGER :: i
        DO i=1, SIZE(file_name)
            file_exists_array(i) = file_exists_scalar(file_name(i))
        END DO
    END FUNCTION file_exists_array

    !> Checks single files existence
    FUNCTION file_exists_scalar(file_name)
        CHARACTER(LEN=*), INTENT(IN) :: file_name
        LOGICAL :: file_exists_scalar
        INQUIRE( FILE=file_name, EXIST=file_exists_scalar )
    END FUNCTION file_exists_scalar

    !> Checks whether latitude (in degrees) is in the range [-90,90]
    LOGICAL ELEMENTAL FUNCTION valid_lat( lat )
        USE parameters_mod, ONLY : real_kind
        REAL(real_kind), INTENT(IN) :: lat

        IF ( ABS(lat) <= 90.0 ) THEN
            valid_lat = .TRUE.
        ELSE
            valid_lat = .FALSE.
        END IF

    END FUNCTION valid_lat


    !> Checks whether longitude (in degrees) is in the range [-180,180]
    LOGICAL ELEMENTAL FUNCTION valid_lon( lon )
        USE parameters_mod, ONLY : real_kind
        REAL(real_kind), INTENT(IN) :: lon

        IF ( ABS(lon) <= 180.0 ) THEN
            valid_lon = .TRUE.
        ELSE
            valid_lon = .FALSE.
        END IF

    END FUNCTION valid_lon


    !> Checks whether depth (in meters) is in the range (0,6371000]
    LOGICAL ELEMENTAL FUNCTION valid_depth( depth )
        USE parameters_mod, ONLY : real_kind, earth_radius
        REAL(real_kind), INTENT(IN) :: depth

        IF ( 0.0 <= depth .AND. depth < earth_radius ) THEN
            valid_depth = .TRUE.
        ELSE
            valid_depth = .FALSE.
        END IF

    END FUNCTION valid_depth

END MODULE checks_mod
