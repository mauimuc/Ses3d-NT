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
!! $Date: 2013-12-03 19:32:32 +0100 (Tue, 03 Dec 2013) $
!! $Author: mauerberger $
!! $Revision: 837 $
!! @copyright GNU General Public License version 3 or later



!> Module handing Ses3d-NT's raw output
MODULE output_raw_mod
    USE parameters_mod, ONLY : real_kind, fnl, fsl, OUTPUT_UNIT, &
                               mpi_master_rank, is_mpi_master_rank
    USE configuration_mod, ONLY : config => configuration
    USE volume_snapshot_mod, ONLY : volume_snapshot_cls

    IMPLICIT NONE

    PRIVATE

    !> A class which manges writing raw binary output
    TYPE, EXTENDS(volume_snapshot_cls) :: output_raw_cls
        INTEGER, ALLOCATABLE :: fh_ !< MPI file handle
        REAL(real_kind), ALLOCATABLE :: field_(:,:,:,:,:,:) !< Temporary field
        CHARACTER(LEN=8) :: attribute_ !< Channel to be written
    !! @deprecated Should be replaced by an instance of channel_cls class
    CONTAINS
        PROCEDURE :: mpi_write => write_volume_snapshot
        PROCEDURE :: print => output_header
        PROCEDURE :: file_name
    END TYPE output_raw_cls

    !> @see from_namelist
    INTERFACE output_raw_cls
        MODULE PROCEDURE from_namelist
    END INTERFACE


    PUBLIC :: parse_output_raw_file, output_raw_cls

CONTAINS

    !========================!
    ! Structure constructors !
    !========================!

    !> Structure constructor to instantiate a output_raw_cls object
    !! from a output_raw_typ
    !> @todo rename output_raw_typ
    !> @todo Having the additional argument attribute is not nice
    !> @param raw_in ...
    !> @param attribute ...
    !> @result output_raw_cls ...
    FUNCTION from_namelist(raw_in, attribute) RESULT(new)
        USE parser_mod, ONLY : output_raw_typ
        TYPE(output_raw_typ), INTENT(INOUT) :: raw_in
        CHARACTER(LEN=8), INTENT(IN) :: attribute
        TYPE(output_raw_cls) :: new

        ! For consistency checks invoke the volume_snapshot_cls class
        ! structure constructor fist
        new%volume_snapshot_cls = volume_snapshot_cls( &
            ts_start=raw_in%timestep_start, &
            ts_end=raw_in%timestep_end, &
            ts_increment=raw_in%timestep_increment, &
            override=raw_in%override, &
            prefix=raw_in%prefix)

        ! Assign attribute
        new%attribute_   = attribute

    END FUNCTION from_namelist



    !=======================!
    ! Type bound procedures !
    !=======================!

    !> Adds raw-output related information to the print procedure
    SUBROUTINE output_header( self, unit )
        CLASS(output_raw_cls), INTENT(IN) :: self
        INTEGER, INTENT(IN), OPTIONAL :: unit ! Passed unit

        INTEGER :: lun ! Local unit number

        ! Determine output-unit
        lun = OUTPUT_UNIT
        IF ( PRESENT( unit ) ) &
            lun = unit

        WRITE( UNIT=lun, FMT='("raw-snapshots: ")' )

        CALL self%volume_snapshot_cls%print(unit)

        WRITE( UNIT=lun, FMT='( 4x, "attribute = ", A )') TRIM(self%attribute_)

    END SUBROUTINE output_header

    !> Returns the raw-output filename for timestep it
    ELEMENTAL CHARACTER(LEN=fnl) FUNCTION file_name(self, it)
        CLASS(output_raw_cls), INTENT(IN) :: self
        INTEGER, INTENT(IN) :: it

        WRITE(UNIT=file_name, FMT='(A, A, "_", I0, ".vol" )') &
            TRIM(self%prefix_), TRIM(self%attribute_), it

    END FUNCTION

    !> Stores the current value for attribute as raw binary file
    !! @todo does not work reliable if a write has not finished before the
    !!       simulation ends
    SUBROUTINE write_volume_snapshot(self, it)
        USE mpi_parallel_io_mod, ONLY : mpi_write_parallel_begin, &
                                        mpi_write_parallel_end
        USE volume_snapshot_mod, ONLY : get_field
        CLASS(output_raw_cls), INTENT(INOUT) :: self
        INTEGER, INTENT(IN) :: it

        ! In case nothing needs to be done -> return to main program fg
        IF ( .NOT. self%process( it ) ) &
            RETURN

        ! Wait till previous write has been finished
        IF( ALLOCATED(self%fh_) ) THEN
            CALL mpi_write_parallel_end( self%fh_, self%field_ )
            DEALLOCATE( self%fh_, &
                        self%field_ )
        END IF

        ! Allocate memory to buffer field
        IF ( .NOT. ALLOCATED(self%field_) ) &
            ALLOCATE( self%field_( 0:config%nx_loc(),&
                                   0:config%ny_loc(),&
                                   0:config%nz_loc(),&
                                   0:config%lpd(),&
                                   0:config%lpd(),&
                                   0:config%lpd() ) )
        ! Buffer field
        self%field_ = get_field( self%attribute_ )

        ! Start asynchronously writing 'field' to file 'file_name'
        self%fh_ = mpi_write_parallel_begin( &
                    file_name=self%file_name(it), &
                    field=self%field_,&
                    override=self%override_ )

        ! In case the last time-step has been reached wait till file is written to disk
        IF ( it == config%nt() ) &
            CALL mpi_write_parallel_end( self%fh_, self%field_ )

    END SUBROUTINE write_volume_snapshot



    !===================!
    ! Module procedures !
    !===================!

    !> Parses config.-file for &output_raw NAMELIST-groups
    SUBROUTINE parse_output_raw_file(output_raw, filename)
        USE parser_mod, ONLY : count_namelists, parse_namelist, &
                               output_raw_typ
        USE bcast_mod, ONLY : bcast_namelist, bcast_integer
        USE string_utilities_mod, ONLY : lower_case
        CHARACTER(LEN=*), INTENT(IN) :: filename
        TYPE(output_raw_cls), INTENT(OUT), ALLOCATABLE :: output_raw(:)
        TYPE(output_raw_typ) :: raw
        INTEGER :: lun, i, n, j, m
        INTEGER, ALLOCATABLE :: n_attr(:)
        CHARACTER(8) :: attribute

        ! Open receiver file
        IF ( is_mpi_master_rank() ) &
            OPEN(NEWUNIT=lun, FILE=filename, ACTION='READ', STATUS='OLD')

        !----------------
        ! Parse receivers
        !----------------
        IF ( is_mpi_master_rank() ) THEN
            ! Count number of NAMELIST groups
            n = count_namelists(unit=lun, nml_cls=raw)

            ALLOCATE(n_attr(n))

            REWIND(UNIT=lun)
            DO i=1, n
                ! Parse NAMELIST
                CALL parse_namelist(unit=lun, nml_cls=raw)
                raw%attributes(:) = lower_case(raw%attributes(:))
                raw%attributes(:) = ADJUSTL(raw%attributes(:))
                n_attr(i) = SIZE( &
                    PACK(raw%attributes(:), raw%attributes(:) /= ''))
            END DO
        END IF

        CALL bcast_integer(i=n, mr=mpi_master_rank)
        IF ( .NOT. ALLOCATED(n_attr) ) &
            ALLOCATE(n_attr(n))

        DO i=1, SIZE(n_attr)
            CALL bcast_integer(i=n_attr(i), mr=mpi_master_rank)
        END DO

        ! Allocate array of receivers
        ALLOCATE( output_raw(SUM(n_attr)) )

        IF ( is_mpi_master_rank() ) &
            REWIND(UNIT=lun)
        m=0
        DO i=1, n
            IF ( is_mpi_master_rank() ) THEN
                ! Parse NAMELIST
                CALL parse_namelist(unit=lun, nml_cls=raw)
                ! Consistency checks
                ! TODO
                IF ( raw%timestep_end == -1 ) &
                    raw%timestep_end = config%nt()
                raw%attributes(:) = lower_case(raw%attributes(:))
                raw%attributes(:) = ADJUSTL(raw%attributes(:))
            END IF
            ! Broadcast values
            CALL bcast_namelist(raw, mpi_master_rank)
            DO j=1, n_attr(i)
                m = m+1
                attribute = raw%attributes(j)
                ! Invoke structure constructor
                output_raw(m) = output_raw_cls(raw,attribute)
            END DO
        END DO

        ! Disconnect unit
        IF ( is_mpi_master_rank() ) &
            CLOSE(UNIT=lun)

    END SUBROUTINE parse_output_raw_file

END MODULE output_raw_mod
