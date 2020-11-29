! Ses3d-NT - simulation of elastic wave propagation in spherical sections
!
! (c) by Stefan Mauerberger <mauerberger@geophysik.uni-muenchen.de>
!
! This program is free software: you can redistribute it and/or modify
! under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.


!> @brief Dummy global parameters without MPI and NetCDF dependencies
!> @details Module providing parameters and procedures for initializing
!! parameters at program start.
!! $Date: 2013-09-20 14:18:50 +0200 (Fri, 20 Sep 2013) $
!! $Author: mauerberger $
!! $Revision: 699 $
!> @copyright GNU General Public License version 3 or later
MODULE parameters_mod
    USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY : REAL32, REAL64, OUTPUT_UNIT, &
                                              ERROR_UNIT, IOSTAT_END

    IMPLICIT NONE

    PRIVATE

    ! Kind type parameters, preconnected output unit and end-of-file value
    PUBLIC :: REAL32, REAL64, OUTPUT_UNIT, ERROR_UNIT, IOSTAT_END

    !> Max number of observables
    !! (see receivers_mod and volume_snapshots_mod)
    INTEGER, PARAMETER, PUBLIC :: max_attributes = 30

    !> Specifying the real kind
    INTEGER, PARAMETER, PUBLIC :: real_kind = REAL32

    !> Specifying the length of file-name stings
    INTEGER, PARAMETER, PUBLIC :: file_name_length = 512
    INTEGER, PARAMETER, PUBLIC :: fnl=file_name_length !< Alias

    !> Specifying the length of formatter stings
    INTEGER, PARAMETER, PUBLIC :: formatter_string_length = 80
    INTEGER, PARAMETER, PUBLIC :: fsl=formatter_string_length !< Alias

    !> Specifying the length of I/O messages
    INTEGER, PARAMETER, PUBLIC :: iomsg_string_length = 255
    INTEGER, PARAMETER, PUBLIC :: isl=iomsg_string_length !< Alias

    !> Specifying the earth radius in meters
    REAL(real_kind), PARAMETER, PUBLIC :: earth_radius = 6371000.0

    !> Retrieve pi at machines precision at compile time
    REAL(real_kind), PARAMETER, PUBLIC :: pi = 2*ACOS( 0.0_real_kind )

END MODULE parameters_mod

