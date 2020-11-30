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
!! $Date: 2013-11-25 20:20:46 +0100 (Mon, 25 Nov 2013) $
!! $Author: mauerberger $
!! $Revision: 814 $
!! @copyright GNU General Public License version 3 or later



!> Module handing Ses3d-NT's netCDF output
MODULE output_netcdf_mod
    USE error_mod, ONLY : warn
    USE parameters_mod, ONLY : real_kind, fnl, fsl, OUTPUT_UNIT, &
                               mpi_master_rank, is_mpi_master_rank
    USE configuration_mod, ONLY : config => configuration
    USE volume_snapshot_mod, ONLY : volume_snapshot_cls

    IMPLICIT NONE

    PRIVATE

    !> A class which manges writing netCDF output
    TYPE, EXTENDS(volume_snapshot_cls) :: output_netcdf_cls
        INTEGER, ALLOCATABLE :: varids_(:) !< netCDF variable IDs
        CHARACTER(8), ALLOCATABLE :: attributes_(:) !< Channels to be written
    !! @deprecated Should be replaced by an array of channel_cls class
    CONTAINS
#ifdef NetCDF
        PROCEDURE, PRIVATE :: def_var
        PROCEDURE, PRIVATE :: put_var
#endif
        PROCEDURE :: write_netcdf
        PROCEDURE :: print => output_header
        PROCEDURE :: file_name
    END TYPE output_netcdf_cls

    !> @see from_namelist
    INTERFACE output_netcdf_cls
        MODULE PROCEDURE from_namelist
    END INTERFACE output_netcdf_cls


    PUBLIC :: output_netcdf_cls, parse_output_netcdf_file

CONTAINS

    !========================!
    ! Structure constructors !
    !========================!

    !> Structure constructor to instantiate a output_netcdf_cls object
    !! from a parser_mod::output_netcdf_typ
    !> @todo rename output_netcdf_typ
    !> @param netcdf_in ...
    !> @result output_netcdf_cls ...
    FUNCTION from_namelist(netcdf_in) RESULT(new)
        USE parser_mod, ONLY : output_netcdf_typ
        TYPE(output_netcdf_typ), INTENT(INOUT) :: netcdf_in
        TYPE(output_netcdf_cls) :: new

        ! For consistency checks invoke the volume_snapshot_cls class
        ! structure constructor fist
        new%volume_snapshot_cls = volume_snapshot_cls( &
            ts_start=netcdf_in%timestep_start, &
            ts_end=netcdf_in%timestep_end, &
            ts_increment=netcdf_in%timestep_increment, &
            override=netcdf_in%override, &
            prefix=netcdf_in%prefix )

        ! Allocate and assign attributes
        new%attributes_ = PACK( netcdf_in%attributes, &
            netcdf_in%attributes(:) /= '' )

    END FUNCTION from_namelist



    !=======================!
    ! Type bound procedures !
    !=======================!

    !> Adds netCDF-output related information to the print procedure
    SUBROUTINE output_header( self, unit )
        CLASS(output_netcdf_cls), INTENT(IN) :: self
        INTEGER, INTENT(IN), OPTIONAL :: unit ! Passed unit

        INTEGER :: i, n, lun ! Local unit number
        CHARACTER(LEN=fsl) :: fmt ! Formatter string

        ! Determine output-unit
        lun = OUTPUT_UNIT
        IF ( PRESENT( unit ) ) &
            lun = unit

        WRITE(UNIT=lun, FMT='("netCDF-snapshots: ")')

        CALL self%volume_snapshot_cls%print(unit)

        n = SIZE(self%attributes_)

        ! Make formatter string
        WRITE(UNIT=fmt, FMT='(A,I0,A)') &
            '( 4x, "attributes = ", ',n ,'(A, :, ", ") )'
        WRITE(UNIT=lun, FMT=fmt ) (TRIM(self%attributes_(i)), i=1, n)

    END SUBROUTINE output_header

    !> Returns the netCDF-output filename for timestep it
    ELEMENTAL CHARACTER(LEN=fnl) FUNCTION file_name(self, it)
        CLASS(output_netcdf_cls), INTENT(IN) :: self
        INTEGER, INTENT(IN) :: it

        CHARACTER(fsl) :: fmt
        INTEGER :: i

        ! Make formatter string
        WRITE(UNIT=fmt,FMT='( A,I0,A )') &
            '(A, ', SIZE(self%attributes_), '(A,"_"), I0, ".nc" )'

        WRITE(UNIT=file_name,FMT=fmt) TRIM(self%prefix_), &
            ( TRIM(self%attributes_(i)), i=1, SIZE(self%attributes_) ), it

    END FUNCTION

#ifdef NetCDF
    !> Subroutine defining the netCDF variables
    SUBROUTINE def_var( self, &
           ncid, dimidx , dimidy, dimidz, dimidt )
        USE netcdf_parallel_io_mod, ONLY : my_nf90_def_var
        CLASS(output_netcdf_cls), INTENT(INOUT) :: self
        INTEGER, INTENT(IN) :: ncid, dimidx, dimidy, dimidz, dimidt

        IF ( .NOT. ALLOCATED( self%varids_ ) ) &
            ALLOCATE( self%varids_( SIZE(self%attributes_) ) )

        CALL my_nf90_def_var( ncid, dimidx, dimidy, dimidz, dimidt, &
                              self%attributes_, self%varids_ )

    END SUBROUTINE def_var

    ! Subroutine for putting the netCDF variable values
    SUBROUTINE put_var( self, ncid )
        USE netcdf_parallel_io_mod, ONLY : my_nf90_put_var
        USE volume_snapshot_mod, ONLY : get_field
        CLASS(output_netcdf_cls), INTENT(IN) :: self
        INTEGER, INTENT(IN) :: ncid
        INTEGER :: i

        DO i=1, SIZE( self%varids_ )
            CALL my_nf90_put_var( ncid, self%varids_(i), &
                                  get_field( self%attributes_(i) ) )
        END DO

    END SUBROUTINE put_var
#endif

    !> Writes a netCDF file
    SUBROUTINE write_netcdf( self, it )
#ifdef NetCDF
        USE netcdf_parallel_io_mod, ONLY : open_netcdf_file, my_nf90_end_def,&
                                           my_nf90_close, add_metadata, nf90_noerr
#endif
        CLASS(output_netcdf_cls), INTENT(INOUT) :: self
        INTEGER, INTENT(IN) :: it
        INTEGER :: ncid, dimid(4), nf90_err
        REAL :: time

        IF ( .NOT. self%process(it) ) &
            RETURN

#ifdef NetCDF
        time = REAL(it*config%dt(), KIND(time))
        CALL open_netcdf_file( self%file_name(it), &
                               ncid, dimid, time, &
                               override=self%override_, &
                               nf90_err=nf90_err )

        IF ( nf90_err /= nf90_noerr ) THEN
            IF ( is_mpi_master_rank() ) &
                CALL warn( 'Error while opening netCDF file ' // &
                            self%file_name(it) )
            RETURN
        END IF

        CALL self%def_var( ncid, dimid(1), dimid(2), dimid(3), dimid(4) )

        CALL add_metadata( ncid, it )

        CALL my_nf90_end_def( ncid )

        CALL self%put_var( ncid )

        CALL my_nf90_close( ncid )
#else
        IF ( is_mpi_master_rank() ) &
            CALL warn( "Can't write netCDF file! &
                &Ses3d-NT is compiled without netCDF support." )
        IF ( config%log_file_true() ) &
            CALL warn( "Can't write netCDF file! &
                &Ses3d-NT is compiled without netCDF support.", &
                unit=config%log_file_unit )
#endif

    END SUBROUTINE write_netcdf



    !===================!
    ! Module procedures !
    !===================!

    !> Parses config.-file for &output_raw NAMELIST-groups
    SUBROUTINE parse_output_netcdf_file(filename, output_netcdf)
        USE parser_mod, ONLY : count_namelists, parse_namelist, &
                               output_netcdf_typ
        USE bcast_mod, ONLY : bcast_namelist, bcast_integer
        USE string_utilities_mod, ONLY : lower_case
        CHARACTER(LEN=*), INTENT(IN) :: filename
        TYPE(output_netcdf_cls), INTENT(OUT), ALLOCATABLE :: output_netcdf(:)
        TYPE(output_netcdf_typ) :: netcdf
        INTEGER :: lun, n, i

#ifndef NetCDF
        ! FIXME: Those messages are written even so no netcdf output is definded
        IF ( is_mpi_master_rank() ) &
            CALL warn( 'No netCDF files will be written. &
                &Ses3d-NT is compiled without netCDF support.' )
        IF ( config%log_file_true() ) &
            CALL warn( "No netCDF files will be written. &
                &Ses3d-NT is compiled without netCDF support.", &
                unit=config%log_file_unit )
#endif

        ! Open config.-file
        IF ( is_mpi_master_rank() ) &
            OPEN(NEWUNIT=lun, FILE=filename, ACTION='READ', STATUS='OLD')

        !----------------
        ! Parse receivers
        !----------------
        IF ( is_mpi_master_rank() ) &
            ! Count number of NAMELIST groups
            n = count_namelists(unit=lun, nml_cls=netcdf)

        ! Broadcast number of netcdf outputs
        CALL bcast_integer(i=n, mr=mpi_master_rank)

        ! Allocate array of netcdf outputs
        ALLOCATE( output_netcdf(n) )

        IF ( is_mpi_master_rank() ) &
            REWIND(UNIT=lun)

        DO i=1, n
            IF ( is_mpi_master_rank() ) THEN
                ! Parse NAMELIST
                CALL parse_namelist(unit=lun, nml_cls=netcdf)
                ! Consistency checks
                ! TODO
                IF ( netcdf%timestep_end == HUGE(1) ) &
                    netcdf%timestep_end = config%nt()
                netcdf%attributes(:) = lower_case(netcdf%attributes(:))
                netcdf%attributes(:) = ADJUSTL(netcdf%attributes(:))
            END IF
            ! Broadcast values
            CALL bcast_namelist(netcdf, mpi_master_rank)
            ! Invoke structure constructor
            output_netcdf(i) = output_netcdf_cls(netcdf)
        END DO

        ! Disconnect unit
        IF ( is_mpi_master_rank() ) &
            CLOSE(UNIT=lun)

    END SUBROUTINE parse_output_netcdf_file

END MODULE output_netcdf_mod
