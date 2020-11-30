!> @file
!! Ses3d-NT - simulation of elastic wave propagation in spherical sections
!!
!! (c) by Andreas Fichtner
!!        Tarje Nissen-Meyer
!!        Heiner Igel
!!        Stefan Mauerberger
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



!> Ses3d-NT's main program
PROGRAM ses3d_nt
    USE parameters_mod, ONLY : real_kind, fsl, fnl, OUTPUT_UNIT, gpl_phrase, &
                               is_mpi_master_rank, init_my_mpi_real, &
                               init_my_mpi_integer
    USE error_mod, ONLY : abort
    USE checks_mod, ONLY : file_exists
! Line control through pre-processor directive
#ifndef MPI_INCLUDE
    ! in case MPI_INCLUDE is NOT set the Fortran mpi module will be used
    USE mpi
#endif
    USE configuration_mod, ONLY : config => configuration, configuration_typ
    USE geometric_paras_mod, ONLY : allocate_geometric_parameters, &
                                    init_geometric_parameters
    USE model_paras_mod, ONLY : allocate_model_parameters, &
                                make_model_parameters, &
                                init_auxiliary_model_parameters, &
                                write_model_parameters, &
                                print_model_parameters, &
                                cfl_loc, f_c_loc
    USE elastic_vars_mod, ONLY : kinetic_energy, v_ptp_glb, v_ptp_loc, &
                                 allocate_elastic_variables
    USE receiver_mod, ONLY : receiver_cls, parse_receiver_file
    USE source_mod, ONLY : source_cls, parse_source_file
    USE grad_mod, ONLY : ses3d_grad, allocate_gradient_variables, gradient_cls, gradient
    USE output_raw_mod, ONLY : output_raw_cls, parse_output_raw_file
    USE output_netcdf_mod, ONLY : output_netcdf_cls, parse_output_netcdf_file
    USE evolution_mod, ONLY : ses3d_evolution, init_temp

    IMPLICIT NONE

#ifdef MPI_INCLUDE
    ! in case MPI_INCLUDE is set use the mpi header file !! DEPRECATED !!
    INCLUDE 'mpif.h'
#endif

    CHARACTER(LEN=fnl) :: conf_file_name
    INTEGER :: i, ierr, it
    CHARACTER(LEN=fsl) :: fmt
    REAL(real_kind) :: ptp(2)
    TYPE(source_cls), ALLOCATABLE :: sources(:) ! Array which holds all source-typ objects
    TYPE(receiver_cls), ALLOCATABLE :: receivers(:) ! Array which holds all receiver-objects
    TYPE(output_netcdf_cls), ALLOCATABLE :: output_netcdf(:)
    TYPE(output_raw_cls), ALLOCATABLE :: output_raw(:)


    ! Initializes the MPI execution environment
    CALL MPI_INIT( ierr )
    ! Initializes my_mpi_real which corresponds to REAL(real_kind)
    CALL init_my_mpi_real()
    ! Initializes my_mpi_integer which corresponds to atomic INTEGER
    CALL init_my_mpi_integer()


    !======================================================================
    ! Give a short notice to stdout at program start
    !======================================================================
    IF ( is_mpi_master_rank() ) &
        WRITE(OUTPUT_UNIT, '(/, A, /)') gpl_phrase


    !======================================================================
    ! Check parameter file name
    !======================================================================

    ! Get number of command line arguments
    IF ( COMMAND_ARGUMENT_COUNT() < 1 ) &
        CALL abort( 'Config-file missing!' )

    ! Retrieve 1st argument that was passed on the command line
    CALL GET_COMMAND_ARGUMENT( 1, conf_file_name )

    ! Check if config-file exists
    IF ( .NOT. file_exists( conf_file_name ) ) &
        CALL abort( 'Config-File not found: ' // TRIM( conf_file_name ) )


    !======================================================================
    ! Parse configuration file
    !======================================================================

    ! Read configuration
    config = configuration_typ( conf_file_name )
    ! Initialize configuration
    CALL config%init()

    ! Open log-files
    IF ( config%log_file_true() ) &
        OPEN( NEWUNIT=config%log_file_unit, &
              FILE=config%log_file_name(), &
              ACTION='WRITE')

    ! XXX this is just a temporary hack
    CALL allocate_geometric_parameters()
    CALL init_geometric_parameters()

    ! XXX this is just a temporary hack
    CALL allocate_model_parameters()
    CALL make_model_parameters()
    CALL init_auxiliary_model_parameters()

    ! XXX this is just a temporary hack
    CALL config%set_fc_max( f_c_loc() )
    CALL config%set_cfl( cfl_loc() )
    IF ( is_mpi_master_rank() ) THEN
        write(*,'("Maximum frequency:", EN12.2, "Hz")') config%fc_glb_max()
        write(*,'("Courant number:", F5.2)') config%cfl_glb()
    END IF

    ! XXX this is just a temporary hack
    CALL allocate_elastic_variables()

    ! Print log file
    IF ( config%log_file_true() ) THEN
        CALL config%print( config%log_file_unit )
        CALL print_model_parameters( config%log_file_unit )
    END IF

    ! XXX this is just a temporary hack
    CALL init_temp()

    ! Parse receivers
    CALL parse_receiver_file( conf_file_name, receivers )
    ! Write receiver-data to log-file
    IF ( config%log_file_true() ) THEN
        DO i=1, SIZE(receivers)
            CALL receivers(i)%print( config%log_file_unit )
        END DO
    END IF
    ! Screen output
    fmt = '( "In rank ", I0, ": Found ", I0, " receivers." )'
    IF ( SIZE( receivers(:) ) > 0 ) &
        WRITE( UNIT=OUTPUT_UNIT, FMT=fmt ) config%mpi_rank(), SIZE( receivers(:) )


    ! Parse sources
    CALL parse_source_file(conf_file_name, sources)
    ! Write source-data to log-file
    IF ( config%log_file_true() ) THEN
        DO i=1, SIZE(sources)
            CALL sources(i)%print( config%log_file_unit )
        END DO
    END IF
    ! Screen output
    fmt = '( "In rank ", I0, ": Found ", I0, " sources." )'
    IF ( SIZE( sources(:) ) > 0 ) &
        WRITE( UNIT=OUTPUT_UNIT, FMT=fmt ) config%mpi_rank(), SIZE( sources(:) )


    ! Parse netcdf output
    CALL parse_output_netcdf_file(conf_file_name, output_netcdf)
    ! Write netcdf output to log-file
    IF ( config%log_file_true() ) THEN
        DO i=1, SIZE(output_netcdf)
            CALL output_netcdf(i)%print( config%log_file_unit )
        END DO
    END IF
    ! Screen output
    fmt = '( "netCDF-output: ", I0, " detected" )'
    IF ( is_mpi_master_rank() .AND. SIZE(output_netcdf(:)) /= 0 ) &
        WRITE( UNIT=OUTPUT_UNIT, FMT=fmt ) SIZE( output_netcdf(:) )


    ! Parse raw output
    CALL parse_output_raw_file(output_raw, filename=conf_file_name)
    ! Write netcdf output to log-file
    IF ( config%log_file_true() ) THEN
        DO i=1, SIZE(output_raw)
            CALL output_raw(i)%print( config%log_file_unit )
        END DO
    END IF
    ! Screen output
    fmt = '( "raw-output: ", I0, " detected" )'
    IF ( is_mpi_master_rank() .AND. SIZE(output_raw(:)) /= 0 ) &
        WRITE( UNIT=OUTPUT_UNIT, FMT=fmt ) SIZE( output_raw(:) )


    ! Flushes log-file-unit to output
    IF ( config%log_file_true() ) &
        FLUSH( UNIT=config%log_file_unit )


    ! Synchronize all ranks
    CALL MPI_BARRIER( MPI_COMM_WORLD, ierr)


    !======================================================================
    ! Start time evolution
    !======================================================================

    fmt = '( A, /, "Processing workflow ''", A, "'' of event ''", A, "''", /, A )'
    ! Write to log-file
    IF ( config%log_file_true() ) &
        WRITE( UNIT=config%log_file_unit, FMT=fmt ) REPEAT('=',50), &
            config%workflow(), config%event_name(), REPEAT('=',50)
    ! Screen output
    IF (config%mpi_rank()==0) &
        WRITE( UNIT=OUTPUT_UNIT, FMT=fmt ) REPEAT('=',50), &
            config%workflow(), config%event_name(), REPEAT('=',50)


    SELECT CASE( config%workflow() )

        !=======================
        ! Forward time stepping
        !=======================
        CASE( 'forward' )

            fmt = '(A,"Iteration ",I0,", ",I0,"%: ",G0," < v < "G0 )'

            DO it=1, config%nt()

                ! Calculate time increment
                CALL ses3d_evolution()
                ! Add sources
                DO i=1, SIZE(sources)
                    CALL sources(i)%add( it )
                END DO

                ! Record synthetics
                CALL receivers(:)%record( it )
                ! Write raw-output
                DO i=1, SIZE(output_raw)
                    CALL output_raw(i)%mpi_write( it )
                END DO
                ! Write netcdf-output
                DO i=1, SIZE(output_netcdf)
                    CALL output_netcdf(i)%write_netcdf( it )
                END DO

                ! Write to log-file
                IF ( config%log_file_true() ) &
                    WRITE(UNIT=config%log_file_unit, FMT=fmt) &
                        '', it, 100*it/config%nt(), v_ptp_loc()
                ! Screen output
                ptp(:) = v_ptp_glb()
                IF (config%mpi_rank()==0) &
                    WRITE( UNIT=OUTPUT_UNIT ,FMT=fmt, ADVANCE='no' ) &
                        ACHAR(13), it, 100*it/config%nt(), ptp

            END DO

            ! Write synthetics
            DO i=1, SIZE(receivers)
                CALL receivers(i)%write()
            END DO


        !=======================
        ! Adjoint time stepping
        !=======================
        CASE( 'adjoint' )

            gradient = gradient_cls(conf_file_name)

            CALL allocate_gradient_variables()

            fmt = '("Iteration ",I0,": ",G0," < v < "G0 )'

            DO it=1, config%nt()

                ! Calculate time increment
                CALL ses3d_evolution()
                ! Add sources
                DO i=1, SIZE(sources)
                    CALL sources(i)%add( it )
                END DO

                ! Record synthetics
                CALL receivers%record( it )
                ! Write raw-output
                DO i=1, SIZE(output_raw)
                    CALL output_raw(i)%mpi_write( it )
                END DO
                ! Write netcdf-output
                DO i=1, SIZE(output_netcdf)
                    CALL output_netcdf(i)%write_netcdf( it )
                END DO

                ! Write to log-file
                IF ( config%log_file_true() ) &
                    WRITE(UNIT=config%log_file_unit, FMT=fmt) it, v_ptp_loc()
                ! Screen output
                ptp(:) = v_ptp_glb()
                IF (config%mpi_rank()==0) &
                    WRITE( UNIT=OUTPUT_UNIT ,FMT=fmt ) it, ptp

                ! Compute Frechet derivatives
                CALL ses3d_grad( it )

            END DO


        !=======================
        ! Model generation
        !=======================
        CASE ( 'model' )
            ! This approach comes with tones of problems. In addition it
            ! exposes some design flaws.
            CALL write_model_parameters()


        CASE DEFAULT
            CALL abort( 'Unknown workflow ' // config%workflow() )

    END SELECT


    ! Synchronize all ranks
    CALL MPI_BARRIER( MPI_COMM_WORLD, ierr)


    fmt = '( /, A, /, "Workflow '''// config%workflow() // ''' successfully completed! ", /, A )'
    ! Write to log-file
    IF ( config%log_file_true() ) &
        WRITE( UNIT=config%log_file_unit, FMT=fmt ) &
            REPEAT('=',50), REPEAT('=',50)
    ! Screen output
    IF (config%mpi_rank()==0) &
        WRITE( UNIT=OUTPUT_UNIT, FMT=fmt ) REPEAT('=',50), REPEAT('=',50)


    IF ( config%log_file_true() ) &
        CLOSE( UNIT=config%log_file_unit ) ! Log-file

    ! Terminates MPI execution environment
    CALL MPI_FINALIZE( ierr )


END PROGRAM ses3d_nt
