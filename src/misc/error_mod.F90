!> @file
!! Ses3d-NT - simulation of elastic wave propagation in spherical sections
!!
!! (c) by Stefan Mauerberger <mauerberger@geophysik.uni-muenchen.de>
!!    and Maksym Melnyk <mmelnyk@geophysik.uni-muenchen.de>
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
!! $Date: 2013-12-03 19:32:32 +0100 (Tue, 03 Dec 2013) $
!! $Author: mauerberger $
!! $Revision: 837 $
!> @copyright GNU General Public License version 3 or later



!> This module handling errors and warnings.
MODULE error_mod
        USE ISO_FORTRAN_ENV, ONLY : ERROR_UNIT
! line control through pre-processor directives
#ifndef MPI_INCLUDE
       !! RECOMMENDED !!
       !! in case MPI_INCLUDE is NOT defined the mpi module will be used
        USE mpi, ONLY : MPI_INITIALIZED, MPI_ABORT, MPI_COMM_SIZE, &
                        MPI_COMM_RANK, MPI_COMM_WORLD
#endif

    IMPLICIT NONE

#ifdef MPI_INCLUDE
       !! DEPRECATED !!
       !! when MPI_INCLUDE is defined, mpi header-files will be included
        INCLUDE 'mpif.h'
#endif

    PRIVATE

    !> @brief Writes warning messages.
    !> @details An interface differentiating between passing a single string 
    !! or an array of strings either with or without passing a unit identifier.
    !! If no unit is specified ERROR_UNIT is assumed. 
    INTERFACE warn 
        MODULE PROCEDURE warn_single 
        MODULE PROCEDURE warn_single_w_unit 
        MODULE PROCEDURE warn_multi 
        MODULE PROCEDURE warn_multi_w_unit 
    END INTERFACE warn


    !> @brief Causes immediate termination of program execution. 
    !> @details An interface differentiating between having a unit identifier 
    !! specified or not. If not ERROR_UNIT is assumed.
    INTERFACE abort
        MODULE PROCEDURE abort_w_unit
        MODULE PROCEDURE abort_
    END INTERFACE abort

    PUBLIC :: abort, warn

CONTAINS

    
    !> @brief Causes immediate termination of program execution. 
    !> @details First, the passed error message is written into a unit 
    !! identifier than the program gets terminated either using the MPI_Abort 
    !! call or Fortran's intrinsic STOP statement.
    !> @param msg Error message (A singe string of arbitrary length)
    !> @param unit A unit specifier
    SUBROUTINE abort_w_unit(msg, unit)
        CHARACTER(LEN=*), INTENT(IN) :: msg
        INTEGER, INTENT(IN) :: unit
        CHARACTER(LEN=*), PARAMETER :: fmt = &
            '("'//ACHAR(13)//'ERROR in rank ", I0, " of ", I0, ": ", A )'
        LOGICAL :: mpi_init
        INTEGER :: my_rank, mpi_err, n_cpus

        CALL MPI_INITIALIZED( mpi_init, mpi_err )
        IF ( .NOT. mpi_init ) THEN
            WRITE( UNIT=unit, FMT='("ERROR: ", A)' ) msg
            STOP 1

        ELSE
            CALL MPI_COMM_SIZE( MPI_COMM_WORLD, n_cpus, mpi_err )
            CALL MPI_COMM_RANK( MPI_COMM_WORLD, my_rank, mpi_err )

            WRITE( UNIT=unit, FMT=fmt ) my_rank, n_cpus, TRIM( msg )

            CALL MPI_ABORT( MPI_COMM_WORLD, 1, mpi_err )

        END IF

    END SUBROUTINE abort_w_unit



    !> Causes immediate termination of the program using abort_w_unit() 
    !! invoking ERROR_UNIT
    !> @param msg Error message (A singe string of arbitrary length)
    SUBROUTINE abort_(msg)
        CHARACTER(LEN=*), INTENT(IN) :: msg
        CALL abort_w_unit(msg, unit=ERROR_UNIT)
    END SUBROUTINE abort_


    
    !> @brief Writes a waring message into a unit identifier. 
    !> @details In case MPI is initialized the calling MPI rank is printed, too.
    !> @param msg The message (a singe string of arbitrary length)
    !> @param unit A unit identifier
    SUBROUTINE warn_single_w_unit(msg, unit)
        CHARACTER(LEN=*), INTENT(IN) :: msg
        INTEGER, INTENT(IN) :: unit
        CHARACTER(LEN=*), PARAMETER :: fmt = &
            '("'//ACHAR(13)//'WARNING in rank ", I0, " of ", I0, ": ", A )'
        LOGICAL :: mpi_init
        INTEGER :: my_rank, mpi_err, n_cpus

        CALL MPI_INITIALIZED( mpi_init, mpi_err )

        IF ( .NOT. mpi_init ) THEN
            WRITE(UNIT=unit, FMT='("WARNING: ", A)') msg

        ELSE
            CALL MPI_COMM_SIZE( MPI_COMM_WORLD, n_cpus, mpi_err )
            CALL MPI_COMM_RANK( MPI_COMM_WORLD, my_rank, mpi_err )

            WRITE(UNIT=unit, FMT=fmt) my_rank, n_cpus, TRIM(msg)

        END IF

    END SUBROUTINE warn_single_w_unit
    


    !> Writes a warning message into unit ERROR_UNIT using warn_single_w_unit()
    !> @param msg The message (a singe string of arbitrary length)
    SUBROUTINE warn_single(msg)
        CHARACTER(LEN=*), INTENT(IN) :: msg
        CALL warn_single_w_unit(msg=msg, unit=ERROR_UNIT)
    END SUBROUTINE warn_single
    

    
    !> Writes multiple waring messages into a unit identifier using 
    !! warn_single_w_unit()
    !> @param msg Array of messages (multiple strings of same length)
    !> @param unit A unit identifier 
    SUBROUTINE warn_multi_w_unit(msg, unit)
        CHARACTER(LEN=*), INTENT(IN) :: msg(:)
        INTEGER, INTENT(IN) :: unit
        INTEGER :: i
        DO i=1, SIZE(msg)
            CALL warn_single_w_unit(msg(i),unit)
        END DO
    END SUBROUTINE warn_multi_w_unit



    !> Writes multiple warnings into unit ERROR_UNIT using warn_multi_w_unit()
    !> @param msg Array of messages (multiple strings of same length)
    SUBROUTINE warn_multi(msg)
        CHARACTER(LEN=*), INTENT(IN) :: msg(:)
        CALL warn_multi_w_unit(msg=msg, unit=ERROR_UNIT)
    END SUBROUTINE warn_multi


END MODULE error_mod
