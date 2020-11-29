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
!! $Date: 2013-11-15 16:19:14 +0100 (Fri, 15 Nov 2013) $
!! $Author: mauerberger $
!! $Revision: 781 $
!! @copyright GNU General Public License version 3 or later



!> MODULE DESCRIPTION
MODULE time_mod
    USE error_mod
    USE parameters_mod
    USE model_mod, ONLY : model_typ

    IMPLICIT NONE

    PRIVATE

    TYPE, EXTENDS( model_typ ) :: time_typ
    CONTAINS
        PROCEDURE, NON_OVERRIDABLE :: nt
        PROCEDURE, NON_OVERRIDABLE :: dt
        PROCEDURE, NON_OVERRIDABLE :: T
        PROCEDURE :: print => print_time
        PROCEDURE, NON_OVERRIDABLE :: set_fc_max
        PROCEDURE, NON_OVERRIDABLE :: fc_loc_max
        PROCEDURE, NON_OVERRIDABLE :: fc_glb_max
        PROCEDURE, NON_OVERRIDABLE :: set_cfl
        PROCEDURE, NON_OVERRIDABLE :: cfl_loc
        PROCEDURE, NON_OVERRIDABLE :: cfl_glb
        PROCEDURE :: date_time
    END TYPE time_typ

    PUBLIC :: time_typ

CONTAINS

    !--------------------------------------------------------------------------
    ! Returns number of time-steps
    !--------------------------------------------------------------------------
    ELEMENTAL INTEGER FUNCTION nt( self )
        CLASS(time_typ), INTENT(IN) :: self
        nt = self%nt_
    END FUNCTION nt


    !--------------------------------------------------------------------------
    ! Returns time-increment [s]
    !--------------------------------------------------------------------------
    ELEMENTAL REAL(real_kind) FUNCTION dt( self )
        CLASS(time_typ), INTENT(IN) :: self
        dt = self%dt_
    END FUNCTION dt


    !--------------------------------------------------------------------------
    ! Returns returns total simulation duration [s]
    !--------------------------------------------------------------------------
    ELEMENTAL REAL(real_kind) FUNCTION T( self )
        CLASS(time_typ), INTENT(IN) :: self
        T = self%dt_ * self%nt_
    END FUNCTION T


    !--------------------------------------------------------------------------
    ! Function returning date and time array
    !--------------------------------------------------------------------------
    PURE FUNCTION date_time( self )
        CLASS(time_typ), INTENT(IN) :: self
        INTEGER :: date_time(8)
        date_time = self%date_time_
    END FUNCTION


    !--------------------------------------------------------------------------
    !
    !--------------------------------------------------------------------------
    REAL(real_kind) FUNCTION fc_glb_max( self )
        CLASS(time_typ), INTENT(IN) :: self
        fc_glb_max = self%fc_glb_max_
    END FUNCTION fc_glb_max


    !--------------------------------------------------------------------------
    !
    !--------------------------------------------------------------------------
    REAL(real_kind) FUNCTION fc_loc_max( self )
        CLASS(time_typ), INTENT(IN) :: self
        fc_loc_max = self%fc_loc_max_
    END FUNCTION fc_loc_max


    !--------------------------------------------------------------------------
    !
    !--------------------------------------------------------------------------
    SUBROUTINE set_fc_max( self, fc_loc_max )
#ifndef MPI_INCLUDE
        !! RECOMMENDED !!
        ! in case MPI_INCLUDE is NOT defined the mpi module will be used
        USE mpi
#else
        !! DEPRECATED !!
        ! when MPI_INCLUDE is defined, mpi header-files will be included
        INCLUDE 'mpif.h'
#endif
        CLASS(time_typ), INTENT(INOUT) :: self
        REAL(real_kind) :: fc_loc_max
        INTEGER :: mpi_err

        ! Assign local values
        self%fc_loc_max_ = fc_loc_max

        ! Take the global minimum
        CALL MPI_ALLREDUCE( self%fc_loc_max_, self%fc_glb_max_, 1, &
                            my_mpi_real, MPI_MIN, self%mpi_comm_cart(), mpi_err )

    END SUBROUTINE set_fc_max


    !--------------------------------------------------------------------------
    !
    !--------------------------------------------------------------------------
    REAL(real_kind) FUNCTION cfl_glb( self )
        CLASS(time_typ), INTENT(IN) :: self
        cfl_glb = self%cfl_glb_
    END FUNCTION cfl_glb


    !--------------------------------------------------------------------------
    !
    !--------------------------------------------------------------------------
    REAL(real_kind) FUNCTION cfl_loc( self )
        CLASS(time_typ), INTENT(IN) :: self
        cfl_loc = self%cfl_loc_
    END FUNCTION cfl_loc


    !--------------------------------------------------------------------------
    !
    !--------------------------------------------------------------------------
    SUBROUTINE set_cfl( self, cfl_loc)
#ifndef MPI_INCLUDE
        !! RECOMMENDED !!
        ! in case MPI_INCLUDE is NOT defined the mpi module will be used
        USE mpi
#else
        !! DEPRECATED !!
        ! when MPI_INCLUDE is defined, mpi header-files will be included
        INCLUDE 'mpif.h'
#endif
        CLASS(time_typ), INTENT(INOUT) :: self
        REAL(real_kind) :: cfl_loc
        INTEGER :: mpi_err

        ! Assign local values
        self%cfl_loc_= cfl_loc

        ! Take the global minimum
        CALL MPI_ALLREDUCE( self%cfl_loc_, self%cfl_glb_, 1, &
                            my_mpi_real, MPI_MAX, self%mpi_comm_cart(), mpi_err )

    END SUBROUTINE set_cfl


    !--------------------------------------------------------------------------
    ! Performs a formatted write of most time-parameters into a unit-number
    !--------------------------------------------------------------------------
    SUBROUTINE print_time( self, unit )
        CLASS(time_typ), INTENT(IN) :: self
        INTEGER, INTENT(IN) :: unit
        CHARACTER(LEN=fsl) :: fmt

        fmt = '(A)'
        WRITE( UNIT=unit, FMT=fmt ) 'Time:'
        fmt = '("    dt=",G0,"s, nt=",I0," => T=",G0,"s")'
        WRITE( UNIT=unit, FMT=fmt ) self%dt_, self%nt_, self%dt_*self%nt_
        fmt = '("    CFL-number (local) = ",G0)'
        WRITE( UNIT=unit, FMT=fmt ) self%cfl_loc_
        fmt = '("    local, maximum exciting frequency = ",G0,"Hz")'
        WRITE( UNIT=unit, FMT=fmt ) self%fc_loc_max_

        fmt = '(4x, "Date: ", I0,"-",I0,"-",I0, &
              & "; Time: " , I0":",I0,":",I0,".",I0 )'
        WRITE( UNIT=unit, FMT=fmt ) self%date_time_(1:3), self%date_time_(5:8)

    END SUBROUTINE print_time


END MODULE time_mod
