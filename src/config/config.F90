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
!! $Date: 2013-12-02 18:27:28 +0100 (Mon, 02 Dec 2013) $
!! $Author: mauerberger $
!! $Revision: 833 $
!! @copyright GNU General Public License version 3 or later



!> Configuration module
MODULE configuration_mod
    USE error_mod, ONLY : warn, abort
    USE parameters_mod
    USE time_mod

    IMPLICIT NONE

    PRIVATE

    TYPE, EXTENDS( time_typ ) :: configuration_typ
    CONTAINS
        PROCEDURE :: print => print_configuration
    END TYPE

    ! object to store configuration in
    TYPE( configuration_typ ) :: configuration

    INTERFACE configuration_typ
        MODULE PROCEDURE parse_conf_file
    END INTERFACE


    PUBLIC :: configuration, configuration_typ

CONTAINS

    SUBROUTINE print_configuration( self, unit )
        CLASS(configuration_typ), INTENT(IN) :: self
        INTEGER, INTENT(IN) :: unit

        WRITE( UNIT=unit, FMT='(A)') REPEAT( '-', 80 )
        WRITE( UNIT=unit, FMT='(A)') 'Configuration:'
        WRITE( UNIT=unit, FMT='(A)') REPEAT( '-', 80 )
        CALL self%grid_typ%print(unit)
        CALL self%model_typ%print(unit)
        CALL self%general_typ%print(unit)
        CALL self%time_typ%print(unit)
        WRITE( UNIT=unit, FMT='(A /)') REPEAT( '-', 80 )

    END SUBROUTINE print_configuration


    !> @brief This function opens, parses and closes the configuration file.
    !! In return it provides an object of type configuration_typ.
    !> @param par_file_name Expects a Ses3d config.-file incl. full path
    !> @return Object of type configuration_typ
    FUNCTION parse_conf_file( conf_file_name ) RESULT( config )
        USE parser_mod, ONLY : parse_namelist, count_namelists, &
                               time_typ, grid_typ, model_typ, general_typ
        USE bcast_mod, ONLY : bcast_namelist
        USE string_utilities_mod, ONLY : i2c
        USE checks_mod, ONLY : file_exists, valid_lat, valid_lon
        TYPE(configuration_typ) :: config
        CHARACTER(LEN=*), INTENT(IN) :: conf_file_name
        INTEGER :: conf_file_unit, n
        TYPE(grid_typ) :: grid
        TYPE(general_typ) :: general
        TYPE(time_typ) :: time
        TYPE(model_typ) :: model

        ! Open config file
        IF ( is_mpi_master_rank() ) THEN
            OPEN(NEWUNIT=conf_file_unit, FILE=conf_file_name, &
                 STATUS='old', ACTION='read' )
        END IF


        !------------
        ! Parse grid
        !------------
        IF ( is_mpi_master_rank() ) THEN
            ! Count number of NAMELIST groups
            n = count_namelists( unit=conf_file_unit, nml_cls=grid)
            IF ( n /= 1 ) &
                CALL abort( 'Found '//i2c(n)//' &grid-groups in config-file!' )
            ! Parse NAMELIST
            REWIND(unit=conf_file_unit)
            CALL parse_namelist(nml_cls=grid, unit=conf_file_unit)
            ! Consistency checks
            IF ( grid%nx < 4 ) &
                CALL abort('Number of x-elements must be greater than 3!')
            IF ( grid%ny < 4 ) &
                CALL abort('Number of y-elements must be greater than 3!')
            IF ( grid%nz < 4 ) &
                CALL abort('Number of z-elements must be greater than 3!')
            IF ( grid%lpd < 2 .OR. grid%lpd > 8 ) &
                CALL abort('Lagrange polynomial degree must be between 3 and 7!')
            IF ( ANY ( 2*grid%taper_width+1 > [ grid%nx, grid%ny, grid%nz ] ) ) &
                CALL abort( 'Relaxing boundaries are too wide!' )
        END IF
        ! Broadcast values
        CALL bcast_namelist(grid, mpi_master_rank)
        ! Assign values
        config%nx_glb_ = grid%nx
        config%ny_glb_ = grid%ny
        config%nz_glb_ = grid%nz
        config%lpd_    = grid%lpd
        config%taper_width = grid%taper_width
        config%taper_slope = grid%taper_slope



        !--------------
        ! Parse general
        !--------------
        IF ( is_mpi_master_rank() ) THEN
            ! Count number of NAMELIST groups
            n = count_namelists( unit=conf_file_unit, nml_cls=general)
            IF ( n/=1 ) &
                CALL abort( 'Found '//i2c(n)//' &general-groups in config-file!' )
            ! Parse NAMELIST
            REWIND(unit=conf_file_unit)
            CALL parse_namelist(nml_cls=general, unit=conf_file_unit)
            ! Consistency checks
            IF ( .NOT. file_exists(general%log_file_dir) ) &
                CALL abort( "Log-file directory '"//TRIM(general%log_file_dir)//"' not found!" )
            IF ( general%event_name == ' ' ) &
                general%event_name = 'ses3d-nt event'
            ! Workflow will implicitly checked in main program
        END IF
        ! Broadcast values
        CALL bcast_namelist(general, mpi_master_rank)
        ! Assign values
        config%workflow_     = general%workflow
        config%log_file_dir_ = TRIM(general%log_file_dir)
        config%event_name_   = general%event_name



        !--------------
        ! Parse model
        !--------------
        IF ( is_mpi_master_rank() ) THEN
            ! Count number of NAMELIST groups
            n = count_namelists( unit=conf_file_unit, nml_cls=model)
            IF ( n/=1 ) &
                CALL abort( 'Found '//i2c(n)//' &model-groups in config-file!' )
            ! Parse NAMELIST
            REWIND(unit=conf_file_unit)
            CALL parse_namelist(nml_cls=model, unit=conf_file_unit)
            ! Consistency checks
            IF ( model%lat_min > model%lat_max ) &
                CALL abort( 'lat_max < lat_min!' )
            IF ( model%lon_min > model%lon_max ) &
                CALL abort( 'lon_max < lon_min!' )
            IF ( model%rad_min > model%rad_max ) &
                CALL abort( 'rad_max < rad_min!' )
            IF ( .NOT. ALL( valid_lat( [model%lat_min, model%lat_max] ) ) ) &
                CALL abort( 'Invalid latitude' )
            IF ( .NOT. ALL( valid_lon( [model%lon_min, model%lon_max] ) ) ) &
                CALL abort( 'Invalid longitude' )
            IF ( ANY( [model%rad_min, model%rad_max] < 1.0 ) ) &
                CALL abort( 'Invalid Radius' )
        END IF
        ! Broadcast values
        CALL bcast_namelist(model, mpi_master_rank)
        ! Assign values
        config%lat_max_    = model%lat_max
        config%lat_min_    = model%lat_min
        config%lon_max_    = model%lon_max
        config%lon_min_    = model%lon_min
        config%z_glb_max_  = model%rad_max
        config%z_glb_min_  = model%rad_min
        config%model_type_ = model%model_type
        config%rhoinv_     = model%rhoinv
        config%lambda_     = model%lambda
        config%mu_         = model%mu
        config%a_          = model%a
        config%b_          = model%b
        config%c_          = model%c
        config%q_          = model%q
        config%override_   = model%override



        !--------------
        ! Parse time
        !--------------
        IF ( is_mpi_master_rank() ) THEN
            ! Count number of NAMELIST groups
            n = count_namelists( unit=conf_file_unit, nml_cls=time)
            IF ( n/=1 ) &
                CALL abort( 'Found '//i2c(n)//' &time-groups in config-file!' )
            ! Parse NAMELIST
            REWIND(unit=conf_file_unit)
            CALL parse_namelist(nml_cls=time, unit=conf_file_unit)
            ! Consistency checks
            IF ( ( 2 > time%nt ) .OR. ( time%nt > 50000 ) ) &
                CALL abort( 'Non-conforming nt value in &time-group!' )
            IF ( time%dt <= 0.0 ) &
                CALL abort( 'Non-conforming dt value in &time-group!' )
            IF ( ANY( time%date_time /= -1 ) ) THEN
                !         date_and_time     mm, dd, D_utc,  h,  m,  s, ms
                IF ( ANY( time%date_time(2:) < [ 1, 1,-1439, 0, 0, 0, 0] ).OR.&
                     ANY( time%date_time(2:) > [12,31, 1439,23,59,59,99] ) ) &
                    CALL abort( 'Inappropriate value in date_time array! ' )
            ELSE
                CALL DATE_AND_TIME( values=time%date_time )
            END IF
        END IF
        ! Broadcast values
        CALL bcast_namelist(time, mpi_master_rank)
        ! Assign values
        config%dt_           = time%dt
        config%nt_           = time%nt
        config%date_time_(:) = time%date_time(:)


        ! Close config.-file
        IF ( is_mpi_master_rank() ) THEN
            CLOSE( UNIT=conf_file_unit )
        END IF

    END FUNCTION parse_conf_file


END MODULE configuration_mod
