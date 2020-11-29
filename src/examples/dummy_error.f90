! Ses3d-NT - simulation of elastic wave propagation in spherical sections
!
! (c) by Stefan Mauerberger <mauerberger@geophysik.uni-muenchen.de>
!    and Maksym Melnyk <mmelnyk@geophysik.uni-muenchen.de>
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


!> Module handling errors and warnings.
!! $Date: 2013-08-21 08:54:40 +0200 (Wed, 21 Aug 2013) $
!! $Author: mauerberger $
!! $Revision: 618 $
!> @copyright GNU General Public License version 3 or later 
MODULE error_mod
        USE ISO_FORTRAN_ENV, ONLY : ERROR_UNIT

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
    !! identifier than the program gets terminated using 
    !! Fortran's intrinsic STOP statement.
    !> @param msg Error message (A singe string of arbitrary length)
    !> @param unit A unit specifier
    SUBROUTINE abort_w_unit(msg, unit)
        CHARACTER(LEN=*), INTENT(IN) :: msg
        INTEGER, INTENT(IN) :: unit

        WRITE( UNIT=unit, FMT='("ERROR: ", A)' ) msg
        STOP 1

    END SUBROUTINE abort_w_unit



    !> Causes immediate termination of the program using abort_w_unit() 
    !! invoking ERROR_UNIT
    !> @param msg Error message (A singe string of arbitrary length)
    SUBROUTINE abort_(msg)
        CHARACTER(LEN=*), INTENT(IN) :: msg
        CALL abort_w_unit(msg, unit=ERROR_UNIT)
    END SUBROUTINE abort_


    
    !> Writes a waring message into a unit identifier. 
    !> @param msg The message (a singe string of arbitrary length)
    !> @param unit A unit identifier
    SUBROUTINE warn_single_w_unit(msg, unit)
        CHARACTER(LEN=*), INTENT(IN) :: msg
        INTEGER, INTENT(IN) :: unit

        WRITE(UNIT=unit, FMT='("WARNING: ", A)') msg

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
