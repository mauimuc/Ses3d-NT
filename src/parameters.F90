!> @file
!! Ses3d-NT - simulation of elastic wave propagation in spherical sections
!!
!! (c) by Stefan Mauerberger
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
!! $Date: 2013-11-20 14:01:31 +0100 (Wed, 20 Nov 2013) $
!! $Author: mauerberger $
!! $Revision: 799 $
!! @copyright GNU General Public License version 3 or later



!> Global parameters
!!
!! Module providing constant parameters and procedures for initializing
!! parameters at program start.
MODULE parameters_mod
    USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY : REAL32, REAL64, REAL128, &
        OUTPUT_UNIT, ERROR_UNIT, IOSTAT_END
#ifdef NetCDF
    USE netcdf, ONLY : NF90_REAL4, NF90_REAL8
#endif

    IMPLICIT NONE

    PRIVATE

    ! Kind type parameters, pre-connected output unit and end-of-file value
    PUBLIC :: REAL32, REAL64, REAL128, OUTPUT_UNIT, ERROR_UNIT, IOSTAT_END

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

    !> MPI-type corresponding to parameters_mod::real_kind.
    !! For initialization call init_my_mpi_real() first
    INTEGER, PUBLIC, PROTECTED :: my_mpi_real = -1
    PUBLIC :: init_my_mpi_real

    !> MPI-type matching Fortran's atomic integer kind.
    !! For initialization call init_my_mpi_real() first
    INTEGER, PUBLIC, PROTECTED :: my_mpi_integer = -1
    PUBLIC :: init_my_mpi_integer

    !> String printed on program start.
    CHARACTER(LEN=*), PARAMETER, PUBLIC :: gpl_phrase = &
        "Ses3d-NT " // NEW_LINE('a') // &
        "Copyright (C) 2013 Stefan Mauerberger " // NEW_LINE('a') // &
        "This program comes with ABSOLUTELY NO WARRANTY. It is free software, &
        &and you are welcome to modify and redistribute it under the terms of &
        &GNU GPL version 3 or later. "

    !> Specifies the MPI master rank
    INTEGER, PARAMETER, PUBLIC :: mpi_master_rank = 0
    PUBLIC :: is_mpi_master_rank

#ifdef NetCDF
    !> The NetCDF real_kind specifier
    INTEGER, PARAMETER, PUBLIC :: my_nf90_real = NF90_REAL4
    !> Fortran's real_kind corresponding to my_nf90_real
    INTEGER, PARAMETER, PUBLIC :: my_nf90_real_kind = REAL32
#endif


CONTAINS

    !> This subroutine is ought to initialize parameters_mod::my_mpi_real as
    !! protected integer.
    !!
    !! This is necessary to guarantee the same floating-point representations
    !! in both Fortran and MPI.
    !!
    !! @exception MPI is not initialized
    SUBROUTINE init_my_mpi_real()
        USE error_mod
! Line control through pre-processor directives
#ifndef MPI_INCLUDE
        !! RECOMMENDED !!
        ! in case MPI_INCLUDE is NOT defined the mpi module will be used
        USE mpi, ONLY : MPI_SIZEOF, MPI_TYPE_MATCH_SIZE, MPI_INITIALIZED, MPI_TYPECLASS_REAL
#else
        !! DEPRECATED !!
        ! when MPI_INCLUDE is defined, mpi header-files will be included
        INCLUDE 'mpif.h'
#endif

        LOGICAL :: initialized = .FALSE.
        INTEGER :: my_mpi_real_kind
        INTEGER :: mpi_err ! no need to handle; all mpi-errors are treated as fatal

        CALL MPI_INITIALIZED( initialized, mpi_err )

        ! MPI_INIT() must be called first
        IF ( .NOT. initialized ) &
            CALL abort( 'Can not determine mpi_real_kind. &
                        & MPI is not yet initialized!' )

        ! For real(REAL_KIND) find corresponding MPI_KIND
! Line control through pre-processor directives
#ifndef MPI_INCLUDE
        CALL MPI_SIZEOF( 0.0_real_kind, &
                         my_mpi_real_kind, &
                         mpi_err )
#else
        ! MPI_SIZEOF procedure is not available in C/C++ header file
        my_mpi_real_kind = STORAGE_SIZE( 0.0_real_kind )/8 ! F08
#endif

        ! Set my_mpi_real
        CALL MPI_TYPE_MATCH_SIZE( MPI_TYPECLASS_REAL, &
                                  my_mpi_real_kind, &
                                  my_mpi_real, &
                                  mpi_err )

    END SUBROUTINE init_my_mpi_real


    !> This subroutine is to initialize parameters_mod::my_mpi_integer as
    !! protected integer.
    !!
    !! This is necessary to guarantee the same integer representations in both
    !! Fortran and MPI.
    !!
    !! @exception MPI is not initialized
    SUBROUTINE init_my_mpi_integer()
        USE error_mod
! line control through pre-processor directives
#ifndef MPI_INCLUDE
        !! RECOMMENDED !!
        ! in case MPI_INCLUDE is NOT defined the mpi module will be used
        USE mpi, ONLY : MPI_SIZEOF, MPI_TYPE_MATCH_SIZE, MPI_INITIALIZED, MPI_TYPECLASS_INTEGER
#else
        !! DEPRECATED !!
        ! when MPI_INCLUDE is defined, mpi header-files will be included
        INCLUDE 'mpif.h'
#endif

        LOGICAL :: initialized = .FALSE.
        INTEGER :: mpi_err ! no need to handle; all mpi-errors are treated as fatal
        INTEGER :: my_mpi_integer_kind

        CALL MPI_INITIALIZED( initialized, mpi_err )

        ! MPI_INIT() must be called first
        IF ( .NOT. initialized ) &
            CALL abort( 'MPI not yet initialized!' )

        ! For default INTEGER find corresponding MPI_KIND
#ifndef MPI_INCLUDE
        CALL MPI_SIZEOF( 0, &
                         my_mpi_integer_kind, &
                         mpi_err )
#else
        ! MPI_SIZEOF procedure is not available in C/C++ header file
        my_mpi_integer_kind = STORAGE_SIZE( 1 )/8 ! F08
#endif

        ! Set my_mpi_integer
        CALL MPI_TYPE_MATCH_SIZE( MPI_TYPECLASS_INTEGER, &
                                  my_mpi_integer_kind, &
                                  my_mpi_integer, &
                                  mpi_err )

    END SUBROUTINE init_my_mpi_integer


    !> Checks if current rank is MPI master rank.
    !!
    !! @result Returns .TRUE. if the MPI rank equals
    !!         parameters_mod::mpi_master_rank  else .FALSE.
    LOGICAL FUNCTION is_mpi_master_rank()
! line control through pre-processor directives
#ifndef MPI_INCLUDE
        !! RECOMMENDED !!
        ! in case MPI_INCLUDE is NOT defined the mpi module will be used
        USE mpi, ONLY : MPI_COMM_RANK, MPI_COMM_WORLD
#else
        !! DEPRECATED !!
        ! when MPI_INCLUDE is defined, mpi header-files will be included
        INCLUDE 'mpif.h'
#endif
        INTEGER :: mpi_err, my_mpi_rank

        CALL MPI_COMM_RANK( MPI_COMM_WORLD, my_mpi_rank, mpi_err )

        IF ( my_mpi_rank == mpi_master_rank ) THEN
            is_mpi_master_rank = .TRUE.
        ELSE
            is_mpi_master_rank = .FALSE.
        END IF

    END FUNCTION is_mpi_master_rank


END MODULE parameters_mod

