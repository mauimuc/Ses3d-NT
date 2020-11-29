!> @file
!! Ses3d-NT - simulation of elastic wave propagation in spherical sections
!!
!! (c) by Stefan Mauerberger <mauerberger@geophysik.uni-muenchen.de>
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
!! $Date: 2013-11-15 21:40:27 +0100 (Fri, 15 Nov 2013) $
!! $Author: mauerberger $
!! $Revision: 784 $
!! @copyright GNU General Public License version 3 or later



!> General module
MODULE general_mod
    USE error_mod
    USE parameters_mod
    USE grid_mod

    IMPLICIT NONE

    PRIVATE

    TYPE, EXTENDS( grid_typ ) :: general_typ
    CONTAINS
        PROCEDURE, NON_OVERRIDABLE :: init_general
        PROCEDURE :: init => init_general
        PROCEDURE :: print => print_general
        PROCEDURE :: log_file_name
        PROCEDURE :: log_file_true
        PROCEDURE :: event_name
        PROCEDURE :: workflow
    END TYPE general_typ

    PUBLIC :: general_typ

CONTAINS

    PURE FUNCTION workflow( self )
        CLASS(general_typ), INTENT(IN) :: self
        CHARACTER(LEN=LEN_TRIM(self%workflow_) ) :: workflow
        workflow = self%workflow_
    END FUNCTION workflow

    PURE FUNCTION event_name( self )
        CLASS(general_typ), INTENT(IN) :: self
        CHARACTER(LEN=LEN_TRIM(self%event_name_) ) :: event_name
        event_name = self%event_name_
    END FUNCTION event_name

    PURE FUNCTION log_file_name( self )
        CLASS(general_typ), INTENT(IN) :: self
        CHARACTER(LEN=LEN_TRIM(self%log_file_name_) ) :: log_file_name
        log_file_name = self%log_file_name_
    END FUNCTION log_file_name

    LOGICAL FUNCTION log_file_true( self )
        CLASS(general_typ), INTENT(IN) :: self

        IF ( self%log_file_dir_ /= 'N/A' ) THEN
            log_file_true = .TRUE.
        ELSE
            log_file_true = .FALSE.
        END IF

    END FUNCTION log_file_true


    !> @brief Subroutine which writes formated configuration into a unit
    SUBROUTINE print_general( self, unit )
        CLASS(general_typ), INTENT(IN) :: self
        INTEGER, INTENT(IN) :: unit
        CHARACTER(LEN=fsl) :: fmt

        fmt = '(A)'
        WRITE( UNIT=unit, FMT=fmt ) 'General:'
        fmt = '("    ",A,A)'
        WRITE( UNIT=unit, FMT=fmt ) 'Workflow = ', self%workflow()
        WRITE( UNIT=unit, FMT=fmt ) 'Event name = ', self%event_name()

        IF ( self%log_file_true() ) THEN
            fmt = '(4x,"Log-file name = ",A)'
            WRITE( UNIT=unit, FMT=fmt ) TRIM( self%log_file_name_ )
        END IF

    END SUBROUTINE print_general


    !> @brief Subroutine which initializes log-file names and if not specified date_time
    SUBROUTINE init_general( self )
        USE string_utilities_mod, ONLY : replace_char
        CLASS(general_typ), INTENT(INOUT) :: self
        CHARACTER(fnl) :: log_file_name

        ! Call parent init procedure first
        CALL self%init_grid()

        ! init log file names
        IF ( self%log_file_true() ) THEN

            ! If set, use event_name for log-file name prefix
            IF ( ALL( self%event_name_ /= [ '   ', 'N/A' ] ) ) THEN
                log_file_name = replace_char( self%event_name_, '/', '-' )
                log_file_name = replace_char( TRIM(log_file_name), ' ', '_' )
            ELSE
                log_file_name = 'logfile'
            END IF

            ! Construct mpi-rank dependent log-file names
            WRITE( UNIT=self%log_file_name_, FMT='(A,A,"_",I0,".log")' ) &
                TRIM(self%log_file_dir_), TRIM(log_file_name), self%mpi_rank()

        END IF

    END SUBROUTINE init_general


END MODULE general_mod
