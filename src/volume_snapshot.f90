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
!! $Date: 2013-12-03 19:32:32 +0100 (Tue, 03 Dec 2013) $
!! $Author: mauerberger $
!! $Revision: 837 $
!! @copyright GNU General Public License version 3 or later



!> Module handing Ses3d-NT's netCDF output
MODULE volume_snapshot_mod
    USE parameters_mod, ONLY : fnl, fsl, real_kind, OUTPUT_UNIT
    USE configuration_mod, ONLY : config => configuration

    IMPLICIT NONE

    PRIVATE

    !> A class which handles general volume-snapshot related things
    TYPE :: volume_snapshot_cls ! Components
    !PRIVATE
        INTEGER :: ts_start_ !< Timestep start
        INTEGER :: ts_end_ !< Timestep stop
        INTEGER :: ts_increment_ !< Timestep increment
        LOGICAL :: override_ !< Override existing files
        CHARACTER(fnl) :: prefix_ !< Filename prefix
    CONTAINS ! Type bound procedures
        PROCEDURE :: process
        PROCEDURE :: print => output_header
        PROCEDURE :: n_file_names
        PROCEDURE :: t_file_names
    END TYPE

    !> @see constructor
    INTERFACE volume_snapshot_cls
        MODULE PROCEDURE constructor
    END INTERFACE

    !> Object which does simply nothing
    TYPE(volume_snapshot_cls) :: null_volume_snapshot_obj = &
        volume_snapshot_cls( -1, -1, -1, .FALSE., '' )

    PUBLIC :: volume_snapshot_cls, get_field

CONTAINS

    !========================!
    ! Structure constructors !
    !========================!

    !> Structure constructor to instantiate a volume_snampshot_cls object
    !! with consistency checks and default values
    !> @todo Write errors and warnings to log-file and flush
    !> @todo Write errors and warnings only to master_rank
    FUNCTION constructor(ts_start, ts_increment, ts_end, &
                         override, prefix) RESULT(new)
        USE error_mod, ONLY : abort
        USE string_utilities_mod, ONLY : i2c
        USE parser_mod, ONLY : default_integer, default_character
        INTEGER, INTENT(INOUT) :: ts_start, ts_end, ts_increment
        LOGICAL, INTENT(IN) :: override
        CHARACTER(LEN=fnl), INTENT(INOUT) :: prefix
        TYPE(volume_snapshot_cls) :: new

        ! In case assign default values
        IF ( ts_end == default_integer ) &
            ts_end = config%nt()
        IF ( ts_start == default_integer ) &
            ts_start = 1
        IF ( prefix == default_character ) &
            prefix = './'

        ! Consistency checks
        IF ( ts_start < 1 ) THEN
            CALL abort( 'timestep_start out of range' )
            new = null_volume_snapshot_obj
            RETURN
        END IF
        IF ( ts_increment < 1 .OR. ts_increment > config%nt() ) THEN
            CALL abort( 'timestep_increment out of range' )
            new = null_volume_snapshot_obj
            RETURN
        END IF
        IF ( ts_end < ts_start ) THEN
            CALL abort( 'timestep_end out of range' )
            new = null_volume_snapshot_obj
            RETURN
        END IF

        ! Invoke the default structure constructor
        new = volume_snapshot_cls(ts_start_=ts_start, &
                                  ts_increment_=ts_increment, &
                                  ts_end_=ts_end, &
                                  override_=override, &
                                  prefix_=prefix )

    END FUNCTION



    !=======================!
    ! Type-Bound-Procedures !
    !=======================!

    !> Returns true if ...
    !> @todo function description and annotation
    LOGICAL ELEMENTAL FUNCTION process(self, it)
        CLASS(volume_snapshot_cls), INTENT(IN) :: self
        INTEGER, INTENT(IN) :: it

        IF ( ( it >= self%ts_start_ ) .AND. &
             ( it <= self%ts_end_ ) .AND. &
             ( MOD( it - self%ts_start_, self%ts_increment_ ) == 0 ) ) THEN
            process = .TRUE.
        ELSE
            process = .FALSE.
        END IF

    END FUNCTION


    !> Writes some general volume-snapshot related information into a unit.
    !> @param unit Optional. Default is OUTPUT_UNIT.
    SUBROUTINE output_header(self, unit)
        CLASS(volume_snapshot_cls), INTENT(IN) :: self
        INTEGER, INTENT(IN), OPTIONAL :: unit ! Passed unit
        INTEGER :: lun ! Local unit number
        CHARACTER(LEN=fsl) :: fmt ! Formatter string

        ! Determine output-unit
        lun = OUTPUT_UNIT
        IF ( PRESENT( unit ) ) &
            lun = unit

        fmt = '( 4x, A, " = ", I0 )'
        WRITE( UNIT=lun, FMT=fmt ) 'timestep_start', self%ts_start_
        WRITE( UNIT=lun, FMT=fmt ) 'timestep_end', self%ts_end_
        WRITE( UNIT=lun, FMT=fmt ) 'timestep_increment', self%ts_increment_
        fmt = '( 4x, A, " = ", L2 )'
        WRITE( UNIT=lun, FMT=fmt ) 'override', self%override_
        fmt = '( 4x, A, " = ", A )'

    END SUBROUTINE

    !> Returns the number of filename written during a simulation
    ELEMENTAL INTEGER FUNCTION n_file_names(self)
        CLASS(volume_snapshot_cls), INTENT(IN) :: self
        INTEGER :: i, time(config%nt())

        time(:) = [(i, i=1, config%nt())]

        n_file_names = COUNT(process(self, time(:)))

    END FUNCTION

    !> Returns the timesteps when files are written
    PURE FUNCTION t_file_names(self)
        CLASS(volume_snapshot_cls), INTENT(IN) :: self
        INTEGER :: t_file_names(n_file_names(self))
        INTEGER :: i, time(config%nt())

        time(:) = [(i, i=1, config%nt())]

        t_file_names(:) = PACK(time(:), process(self, time(:)))

    END FUNCTION



    !===================!
    ! Module procedures !
    !===================!

    !> Copies the demanded field 'name' to 'field'
    !> @deprecated Should be replaced by channel_cls class
    PURE FUNCTION get_field( name ) RESULT( field )
        USE string_utilities_mod, ONLY : lower_case
        USE elastic_vars_mod, ONLY : vx, vy, vz, exy, exz, eyx, eyz, ezy, ezx,&
                                     sxx_pml, syy_pml, szz_pml, dxux, dxuy, &
                                     dxuz, dyux, dyuy, dyuz, dzux, dzuy, dzuz
        USE model_paras_mod, ONLY : rho, lambda, mu, a, b, c, qq
        USE grad_mod, ONLY : grad_rho, grad_cs, grad_cp, grad_csv, grad_csh
        USE geometric_paras_mod, ONLY : r, sin_theta, cot_theta, r_sin_theta, &
                                        taper
        CHARACTER(*), INTENT(IN) :: name
        REAL(real_kind) :: field( 0:config%nx_loc(), &
                                  0:config%ny_loc(), &
                                  0:config%nz_loc(), &
                                  0:config%lpd(), &
                                  0:config%lpd(), &
                                  0:config%lpd() )

        SELECT CASE( lower_case( TRIM(name) ) )

            CASE( 'vx' )
                field(:,:,:,:,:,:) = vx

            CASE( 'vy' )
                field(:,:,:,:,:,:) = vy

            CASE( 'vz' )
                field(:,:,:,:,:,:) = vz

            CASE( 'exx' )
                field(:,:,:,:,:,:) = dxux

            CASE( 'exy' )
                field(:,:,:,:,:,:) = ( dxuy + dyux ) / 2.0

            CASE( 'exz' )
                field(:,:,:,:,:,:) = ( dxuz + dzux ) / 2.0

            CASE( 'eyy' )
                field(:,:,:,:,:,:) = dyuy

            CASE( 'eyz' )
                field(:,:,:,:,:,:) = ( dyuz + dzuy ) / 2.0

            CASE( 'ezz' )
                field(:,:,:,:,:,:) = dzuz

            CASE( 'rot_vx' )
                field(:,:,:,:,:,:) = ( -vy / r + eyz / r_sin_theta + ezy )

            CASE( 'rot_vy' )
                field(:,:,:,:,:,:) = ( -ezx + ( vx - exz ) / r )

            CASE( 'rot_vz' )
                field(:,:,:,:,:,:) = ( exy + vy*cot_theta - eyx/sin_theta ) / r

            CASE( 'pressure' )
                ! XXX is this pressure or pressure-rate?
                field(:,:,:,:,:,:) = -1.0*( sxx_pml + syy_pml + szz_pml )/3.0

            CASE( 'taper' )
                field(:,:,:,:,:,:) = taper

            CASE( 'rho' )
                field(:,:,:,:,:,:) = rho

            CASE( 'lambda' )
                field(:,:,:,:,:,:) = lambda

            CASE( 'vs', 'cs' )
                field(:,:,:,:,:,:) = SQRT( mu / rho )

            CASE( 'vp', 'cp' )
                field(:,:,:,:,:,:) = SQRT( (lambda + 2*mu) / rho )

            CASE( 'mu' )
                field(:,:,:,:,:,:) = mu

            CASE( 'a' )
                field(:,:,:,:,:,:) = a

            CASE( 'b' )
                field(:,:,:,:,:,:) = b

            CASE( 'c' )
                field(:,:,:,:,:,:) = c

            CASE( 'q' )
                field(:,:,:,:,:,:) = qq

            CASE( 'grad_rho' )
                field(:,:,:,:,:,:) = grad_rho

            CASE( 'grad_cs' )
                field(:,:,:,:,:,:) = grad_cs

            CASE( 'grad_cp' )
                field(:,:,:,:,:,:) = grad_cp

            CASE( 'grad_csv')
                field(:,:,:,:,:,:) = grad_csv

            CASE( 'grad_csh')
                field(:,:,:,:,:,:) = grad_csh

        END SELECT

    END FUNCTION get_field

END MODULE volume_snapshot_mod
