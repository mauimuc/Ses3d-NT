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



!> MODEL PARAMETERS
!!
!! All parameters are expected to be in SI unit system:
!!  * rho [kg/m^3]
!!  * cp, ... [m/s]
!!  * mu, lambda [Nm]
!! ...
MODULE model_paras_mod
   USE parameters_mod

    IMPLICIT NONE

    ! Model parameters
    REAL(real_kind), PROTECTED, ALLOCATABLE :: rhoinv(:,:,:,:,:,:), &
                                               mu(:,:,:,:,:,:), &
                                               lambda(:,:,:,:,:,:), &
                                               a(:,:,:,:,:,:), &
                                               b(:,:,:,:,:,:), &
                                               c(:,:,:,:,:,:), &
                                               qq(:,:,:,:,:,:)
    ! Auxiliary model parameters
    REAL(real_kind), PROTECTED, ALLOCATABLE :: kappa(:,:,:,:,:,:), &
                                               two_mu(:,:,:,:,:,:), &
                                               cp(:,:,:,:,:,:), &
                                               cs(:,:,:,:,:,:), &
                                               rho(:,:,:,:,:,:), &
                                               tau(:,:,:,:,:,:), &
                                               mu_tau(:,:,:,:,:,:)
    ! Mass matrix
    REAL(real_kind), PROTECTED, ALLOCATABLE :: mm(:,:,:,:,:,:)

CONTAINS

    !--------------------------------------------------------------------------
    ! allocate memory for model parameters
    !--------------------------------------------------------------------------
    SUBROUTINE allocate_model_parameters
        USE configuration_mod, ONLY : config => configuration

        ASSOCIATE ( nx  => config%nx_loc(), &
                    ny  => config%ny_loc(), &
                    nz  => config%nz_loc(), &
                    lpd => config%lpd()   )

        ! Allocate memory for model parameters
        ALLOCATE( rhoinv(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd), SOURCE=0.0_real_kind )
        ALLOCATE(     mu(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd), SOURCE=0.0_real_kind )
        ALLOCATE( lambda(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd), SOURCE=0.0_real_kind )
        ALLOCATE(      a(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd), SOURCE=0.0_real_kind )
        ALLOCATE(      b(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd), SOURCE=0.0_real_kind )
        ALLOCATE(      c(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd), SOURCE=0.0_real_kind )
        ALLOCATE( qq(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd), SOURCE=HUGE(0.0_real_kind) )
        ! Allocate memory for auxiliary model parameters
        ALLOCATE(  kappa(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd), SOURCE=0.0_real_kind )
        ALLOCATE( two_mu(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd), SOURCE=0.0_real_kind )
        ALLOCATE(     cp(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd), SOURCE=0.0_real_kind )
        ALLOCATE(     cs(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd), SOURCE=0.0_real_kind )
        ALLOCATE(    rho(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd), SOURCE=0.0_real_kind )
        ALLOCATE(    tau(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd), SOURCE=0.0_real_kind )
        ALLOCATE( mu_tau(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd), SOURCE=0.0_real_kind )
        ! Allocate memory for mass-matrix
        ALLOCATE(     mm(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd), SOURCE=0.0_real_kind )

        END ASSOCIATE

    END SUBROUTINE allocate_model_parameters


    !--------------------------------------------------------------------------
    ! Subroutine for reading model parameters
    !--------------------------------------------------------------------------
    SUBROUTINE read_model_parameters( )
        USE mpi_parallel_io_mod, ONLY : mpi_read_parallel
        USE configuration_mod, ONLY : config => configuration

        ! TODO testing for N/A isn't a nice solution
        ! TODO no WRITE - except for exceptions - must occur here

        ! Elastic parameters
        IF ( config%mpi_rank() == 0 ) &
            WRITE(*,*) 'reading: ', config%file_name_rhoinv()
        rhoinv(:,:,:,:,:,:) = mpi_read_parallel( config%file_name_rhoinv() )
        IF ( config%mpi_rank() == 0 ) &
            WRITE(*,*) 'reading: ', config%file_name_mu()
        mu(:,:,:,:,:,:)     = mpi_read_parallel( config%file_name_mu() )
        IF ( config%mpi_rank() == 0 ) &
            WRITE(*,*) 'reading: ', config%file_name_lambda()
        lambda(:,:,:,:,:,:) = mpi_read_parallel( config%file_name_lambda() )
        ! Parameters for an-isotropy
        IF ( config%file_name_a() /= 'N/A' ) THEN
            IF ( config%mpi_rank() == 0 ) &
                WRITE(*,*) 'reading: ', config%file_name_a()
            a(:,:,:,:,:,:) = mpi_read_parallel( config%file_name_a() )
        END IF
        IF ( config%file_name_b() /= 'N/A' ) THEN
            IF ( config%mpi_rank() == 0 ) &
                WRITE(*,*) 'reading: ', config%file_name_b()
            b(:,:,:,:,:,:) = mpi_read_parallel( config%file_name_b() )
        END IF
        IF ( config%file_name_c() /= 'N/A' ) THEN
            IF ( config%mpi_rank() == 0 ) &
                WRITE(*,*) 'reading: ', config%file_name_c()
            c(:,:,:,:,:,:) = mpi_read_parallel( config%file_name_c() )
        END IF
        ! Attenuation
        IF ( config%file_name_q() /= 'N/A' ) THEN
            IF ( config%mpi_rank() == 0 ) &
                WRITE(*,*) 'reading: ', config%file_name_q()
            qq(:,:,:,:,:,:) = mpi_read_parallel( config%file_name_q() )
        END IF

    END SUBROUTINE read_model_parameters


    !--------------------------------------------------------------------------
    ! Subroutine for writing model parameters
    !--------------------------------------------------------------------------
    SUBROUTINE write_model_parameters
        USE mpi_parallel_io_mod, ONLY : mpi_write_parallel
        USE configuration_mod, ONLY : config => configuration

        ! TODO testing for N/A isn't a nice solution
        ! TODO no WRITE - except for exceptions - must occur here

        ! Elastic parameters
        IF ( config%file_name_rhoinv() /= 'N/A' ) THEN
            CALL mpi_write_parallel( file_name=config%file_name_rhoinv(), &
                                     field=rhoinv, &
                                     override=config%override() )
            IF ( config%mpi_rank() == 0 ) &
                WRITE(*,*) 'writing: ', config%file_name_rhoinv()
        END IF
        IF ( config%file_name_mu() /= 'N/A' ) THEN
            CALL mpi_write_parallel( file_name=config%file_name_mu(), &
                                     field=mu, &
                                     override=config%override() )
            IF ( config%mpi_rank() == 0 ) &
                WRITE(*,*) 'writing: ', config%file_name_mu()
        END IF
        IF ( config%file_name_lambda() /= 'N/A' ) THEN
            CALL mpi_write_parallel( file_name=config%file_name_lambda(), &
                                     field=lambda, &
                                     override=config%override() )
            IF ( config%mpi_rank() == 0 ) &
                WRITE(*,*) 'writing: ', config%file_name_lambda()
        END IF
        ! Parameters for an-isotropy
        IF ( config%file_name_a() /= 'N/A' ) THEN
            CALL mpi_write_parallel( config%file_name_a(), &
                                     field=a, &
                                     override=config%override() )
            IF ( config%mpi_rank() == 0 ) &
                WRITE(*,*) 'writing: ', config%file_name_a()
        END IF
        IF ( config%file_name_b() /= 'N/A' ) THEN
            CALL mpi_write_parallel( config%file_name_b(), &
                                     field=b, &
                                     override=config%override() )
            IF ( config%mpi_rank() == 0 ) &
                WRITE(*,*) 'writing: ', config%file_name_b()
        END IF
        IF ( config%file_name_c() /= 'N/A' ) THEN
            CALL mpi_write_parallel( config%file_name_c(), &
                                     field=c, &
                                     override=config%override() )
            IF ( config%mpi_rank() == 0 ) &
                WRITE(*,*) 'writing: ', config%file_name_c()
        END IF
        ! Attenuation
        IF ( config%file_name_q() /= 'N/A' ) THEN
            CALL mpi_write_parallel( config%file_name_q(), &
                                     field=qq, &
                                     override=config%override() )
            IF ( config%mpi_rank() == 0 ) &
                WRITE(*,*) 'writing: ', config%file_name_q()
        END IF

    END SUBROUTINE write_model_parameters


    SUBROUTINE make_model_parameters( )
        USE configuration_mod, ONLY : config => configuration
        USE geometric_paras_mod, ONLY : r
        USE homogeneous_mod, ONLY : homogeneous
        USE ak135_f_mod, ONLY : ak135_f
        USE prem_mod, ONLY : prem, prem_iso
        REAL(real_kind), ALLOCATABLE :: cph(:,:,:,:,:,:), &
                                        cpv(:,:,:,:,:,:), &
                                        csh(:,:,:,:,:,:), &
                                        csv(:,:,:,:,:,:), &
                                        eta(:,:,:,:,:,:)


        SELECT CASE ( config%model_type() )

            ! Homogeneous model
            CASE ( 1 )
                IF ( config%mpi_rank() == 0 ) &
                    WRITE(*,*) 'Homogeneous'
                CALL homogeneous( RHO=rho, LAMBDA=lambda, MU=mu, &
                                  A=a, B=b, C=c )
                                  !Qkappa=, qmu= )
                rhoinv(:,:,:,:,:,:) = 1.0/rho


            ! PREM
            CASE( 11 )
                IF ( config%mpi_rank() == 0 ) &
                    WRITE(*,*) 'PREM anisotropic with crust and with ocean'
                cph = rho; cpv = rho; csh = rho; csv = rho; eta = rho
                CALL prem( R=r, RHO=rho, VPV=cpv, VPH=cph, &
                           VSV=csv, VSH=csh, ETA=eta )
                ! Transform into Ses3d's parametrization
                rho(:,:,:,:,:,:) = 1000.0 * rho
                cph(:,:,:,:,:,:) = 1000.0 * cph
                cpv(:,:,:,:,:,:) = 1000.0 * cpv
                csh(:,:,:,:,:,:) = 1000.0 * csh
                csv(:,:,:,:,:,:) = 1000.0 * csv
                rhoinv(:,:,:,:,:,:) = 1.0 / rho
                mu(:,:,:,:,:,:) = rho * csh**2
                lambda(:,:,:,:,:,:) = rho * cph**2 - 2.0 * mu
                a(:,:,:,:,:,:) = rho * cph**2 - lambda - 2.0 * mu
                b(:,:,:,:,:,:) = rho * csv**2 - mu
                c(:,:,:,:,:,:) = eta * ( lambda + a ) - lambda

            ! PREM isotropic with crust and WITH ocean
            CASE (12)
                IF ( config%mpi_rank() == 0 ) &
                    WRITE(*,*) 'PREM isotropic with crust and with ocean'
                CALL prem_iso( R=r, RHO=rho, VS=cs, VP=cp, &
                               nocrust=.FALSE., noocean=.FALSE.)
                ! Transform into Ses3d's parametrization
                rho(:,:,:,:,:,:) = 1000.0 * rho
                cp(:,:,:,:,:,:) = 1000.0 * cp
                cs(:,:,:,:,:,:) = 1000.0 * cs
                rhoinv(:,:,:,:,:,:) = 1.0 / rho
                mu(:,:,:,:,:,:) = rho * cs**2
                lambda(:,:,:,:,:,:) = rho * cp**2 - 2.0 * mu
                a(:,:,:,:,:,:) = 0.0
                b(:,:,:,:,:,:) = 0.0
                c(:,:,:,:,:,:) = 0.0

            ! PREM with crust and WITHOUT ocean
            CASE (3,13)
                IF ( config%mpi_rank() == 0 ) &
                    WRITE(*,*) 'PREM isotropic with crust and withOUT ocean'
                CALL prem_iso( R=r, RHO=rho, VS=cs, VP=cp, &
                               nocrust=.FALSE., noocean=.TRUE.)
                ! Transform into Ses3d's parametrization
                rho(:,:,:,:,:,:) = 1000.0 * rho
                cp(:,:,:,:,:,:) = 1000.0 * cp
                cs(:,:,:,:,:,:) = 1000.0 * cs
                rhoinv(:,:,:,:,:,:) = 1.0 / rho
                mu(:,:,:,:,:,:) = rho * cs**2
                lambda(:,:,:,:,:,:) = rho * cp**2 - 2.0 * mu
                a(:,:,:,:,:,:) = 0.0
                b(:,:,:,:,:,:) = 0.0
                c(:,:,:,:,:,:) = 0.0

            ! PREM simplified crust
            CASE (14)
                IF ( config%mpi_rank() == 0 ) &
                    WRITE(*,*) 'PREM isotropic simplyfied crust without ocean'
                CALL prem_iso( R=r, RHO=rho, VS=cs, VP=cp, &
                               nocrust=.TRUE., noocean=.TRUE.)
                ! Transform into Ses3d's parametrization
                rho(:,:,:,:,:,:) = 1000.0 * rho
                cp(:,:,:,:,:,:) = 1000.0 * cp
                cs(:,:,:,:,:,:) = 1000.0 * cs
                rhoinv(:,:,:,:,:,:) = 1.0 / rho
                mu(:,:,:,:,:,:) = rho * cs**2
                lambda(:,:,:,:,:,:) = rho * cp**2 - 2.0 * mu
                a(:,:,:,:,:,:) = 0.0
                b(:,:,:,:,:,:) = 0.0
                c(:,:,:,:,:,:) = 0.0

            ! PREM simplified crust (by A. Fichtner)
            CASE (2,15)
                IF ( config%mpi_rank() == 0 ) &
                    WRITE(*,*) 'PREM isotropic simplyfied crust (by A. Fichtner)'
                CALL prem_iso( R=r, RHO=rho, VS=cs, VP=cp, &
                               nocrust=.TRUE., noocean=.TRUE.)
                ! Replace values in first element by constant values
                ! this causes a discontinuity between first an second element
                rho(:,:,0,:,:,:) = 3.00
                cp(:,:,0,:,:,:) = 7.20
                cs(:,:,0,:,:,:) = 4.10
                ! Transform into Ses3d's parametrization
                rho(:,:,:,:,:,:) = 1000.0 * rho
                cp(:,:,:,:,:,:) = 1000.0 * cp
                cs(:,:,:,:,:,:) = 1000.0 * cs
                rhoinv(:,:,:,:,:,:) = 1.0 / rho
                mu(:,:,:,:,:,:) = rho * cs**2
                lambda(:,:,:,:,:,:) = rho * cp**2 - 2.0 * mu
                a(:,:,:,:,:,:) = 0.0
                b(:,:,:,:,:,:) = 0.0
                c(:,:,:,:,:,:) = 0.0

            CASE ( 21 )
                IF ( config%mpi_rank() == 0 ) &
                    WRITE(*,*) 'AK135 full '
                ! Earth radius' in ak135 is in kilometers not meters
                CALL ak135_f( r/1000.0, RHO=rho, VP=cp, VS=cs )!, &
                           !QMU=qmu, QKAPPA=qkappa )
                ! Transform into Ses3d's parametrization
                rho(:,:,:,:,:,:) = 1000.0 * rho
                cp(:,:,:,:,:,:) = 1000.0 * cp
                cs(:,:,:,:,:,:) = 1000.0 * cs
                rhoinv(:,:,:,:,:,:) = 1.0 / rho
                mu(:,:,:,:,:,:) = rho * cs**2
                lambda(:,:,:,:,:,:) = rho * cp**2 - 2.0 * mu

            CASE DEFAULT
                CALL read_model_parameters( )

        END SELECT

    END SUBROUTINE make_model_parameters


    !--------------------------------------------------------------------------
    ! Subroutine printing values of local model-parameters
    ! Optional a unit number may be specified
    !--------------------------------------------------------------------------
    SUBROUTINE print_model_parameters( unit )
        USE, INTRINSIC :: ISO_FORTRAN_ENV , ONLY : OUTPUT_UNIT
        INTEGER, INTENT(IN), OPTIONAL :: unit
        INTEGER :: lun ! Local unit number
        CHARACTER(LEN=fsl) :: fmt ! Formatter string

        ! Determine output-unit
        lun = OUTPUT_UNIT
        IF ( PRESENT( unit ) ) &
            lun = unit

        WRITE( UNIT=lun, FMT='(A)' ) 'Local model parameters: '

        ! Make formatter string
        fmt = '( 4x, A, ":  min = ", G0, ", max = ", G0 )'
        WRITE( UNIT=lun, FMT=fmt ) 'rho',    minval( rho    ), &
                                             maxval( rho    )
        WRITE( UNIT=lun, FMT=fmt ) 'lambda', minval( lambda ), &
                                             maxval( lambda )
        WRITE( UNIT=lun, FMT=fmt ) 'mu',     minval( mu     ), &
                                             maxval( mu     )

        WRITE( UNIT=lun, FMT=fmt ) 'a', minval( a ), &
                                        maxval( a )
        WRITE( UNIT=lun, FMT=fmt ) 'b', minval( b ), &
                                        maxval( b )
        WRITE( UNIT=lun, FMT=fmt ) 'c', minval( c ), &
                                        maxval( c )

        WRITE( UNIT=lun, FMT=fmt ) 'q', minval( qq ), &
                                        maxval( qq )

    END SUBROUTINE print_model_parameters


    !--------------------------------------------------------------------------
    ! Initialization of auxiliary model parameters
    ! Initialization of mass matrix
    !--------------------------------------------------------------------------
    SUBROUTINE init_auxiliary_model_parameters
        USE geometric_paras_mod, ONLY : r, sin_theta, wx_wy_wz, jac
        USE communicate_fields_mod, ONLY : communicate_local_boundaries, &
                                           communicate_global_boundaries

        ! Just to save some floating point operations
        two_mu(:,:,:,:,:,:) = 2.0 * mu
        ! generates kappa and seismic velocities
        kappa(:,:,:,:,:,:) = lambda + 2.0 * mu / 3.0
        WHERE ( rhoinv > 0.0 )
            rho(:,:,:,:,:,:) = 1.0 / rhoinv
        ELSE WHERE
            rho(:,:,:,:,:,:) = 0.0
        END WHERE
        cp(:,:,:,:,:,:) = SQRT( ( lambda + 2.0 * mu ) * rhoinv )
        cs(:,:,:,:,:,:) = SQRT( mu * rhoinv )

        ! XXX where do those values come from
        tau(:,:,:,:,:,:) = 6.0814 * ( qq**(-1.015) )
        mu_tau(:,:,:,:,:,:) = ( 1.0 + tau ) * mu

        ! Generates local mass matrices
        mm(:,:,:,:,:,:) = -jac * r**2 * sin_theta * rho * wx_wy_wz
        ! Makes mass matrix global
        CALL communicate_local_boundaries( MM )
        CALL communicate_global_boundaries( MM )

    END SUBROUTINE init_auxiliary_model_parameters


    !--------------------------------------------------------------------------
    ! returns local minimal and maximal model velocities
    !--------------------------------------------------------------------------
    FUNCTION model_velocities_ptp_loc() RESULT( ptp )
        REAL(real_kind) :: ptp(2)

        ptp(1) = MINVAL( MERGE( cp, cs, cs > cp ) )
        ptp(2) = MAXVAL( MERGE( cp, cs, cs < cp ) )

    END FUNCTION model_velocities_ptp_loc


    !--------------------------------------------------------------------------
    ! Returns global minimum and maximum model velocities
    !--------------------------------------------------------------------------
    ! TODO check whether my_mpi_real and config%mpi_comm_cart() are initialized
    FUNCTION model_velocities_ptp_glb() result( ptp_glb )
        USE configuration_mod, ONLY : config => configuration
#ifndef MPI_INCLUDE
        !! RECOMMENDED !!
        ! in case MPI_INCLUDE is NOT defined the mpi module will be used
        USE mpi
#else
        !! DEPRECATED !!
        ! when MPI_INCLUDE is defined, mpi header-files will be included
        INCLUDE 'mpif.h'
#endif
        REAL(real_kind) :: ptp_loc(2), ptp_glb(2)
        INTEGER :: mpi_err

        ptp_loc(:) = model_velocities_ptp_loc()

        CALL MPI_ALLREDUCE( ptp_loc(2), ptp_glb(2), 1, &
                            my_mpi_real, MPI_MAX, config%mpi_comm_cart(), mpi_err )

        CALL MPI_ALLREDUCE( ptp_loc(1), ptp_glb(1), 1, &
                            my_mpi_real, MPI_MIN, config%mpi_comm_cart(), mpi_err )

    END FUNCTION model_velocities_ptp_glb


    !--------------------------------------------------------------------------
    ! Returns the local maximum frequency corresponding to a sampling rate of
    ! two elements per minimal wavelength
    !--------------------------------------------------------------------------
    FUNCTION f_c_loc()
        USE configuration_mod, ONLY : config => configuration
        USE geometric_paras_mod, ONLY : r, r_sin_theta
        REAL(real_kind) :: s, f_c_loc, f_c(3)
        REAL(real_kind), ALLOCATABLE :: v_mean(:,:,:,:,:), arc_lgt(:,:,:,:,:)

        ! Sampling: two elements per minimal wavelength
        s = 2.0

        ! Calculates the elements arc-length in x-direction [m]
        arc_lgt = r(:,:,:,0,:,:) * config%dx_elm()
        ! Elements' mean velocity in x-direction [m/s]
        ! Velocity here meas the minimum out of vp and vs
        ! This is obtained by merging cs with cp
        v_mean = SUM( MERGE( cp, cs, cs > cp ), DIM=4 ) / SIZE( cs, DIM=4 )
        f_c(1) = MINVAL( v_mean/s/arc_lgt, v_mean > 0.0 )

        ! Calculate arc-length in phi-direction [m]
        arc_lgt = r_sin_theta(:,:,:,:,0,:) * config%dy_elm()
        ! Mean-velocities in y-direction
        v_mean = SUM( MERGE( cp, cs, cs > cp ), DIM=5 ) / SIZE( cs, DIM=5 )
        f_c(2) = MINVAL( v_mean/s/arc_lgt, v_mean > 0.0 )

        ! Mean-velocities in z-direction
        v_mean = SUM( MERGE( cp, cs, cs > cp ), DIM=6 ) / SIZE( cs, DIM=6 )
        f_c(3) = MINVAL( v_mean/s/config%dz_elm(), v_mean > 0.0 )

        ! Return minimum
        f_c_loc = MINVAL( f_c )

    END FUNCTION f_c_loc

    !--------------------------------------------------------------------------
    ! Returns the global maximum frequency
    ! Collective and blocking! Can only be called on all ranks simultaneously
    !--------------------------------------------------------------------------
    ! TODO check whether my_mpi_real and config%mpi_comm_cart() are initialized
    FUNCTION f_c_glb()
        USE configuration_mod, ONLY : config => configuration
#ifndef MPI_INCLUDE
        !! RECOMMENDED !!
        ! in case MPI_INCLUDE is NOT defined the mpi module will be used
        USE mpi
#else
        !! DEPRECATED !!
        ! when MPI_INCLUDE is defined, mpi header-files will be included
        INCLUDE 'mpif.h'
#endif
        REAL(real_kind) :: fc_loc, f_c_glb
        INTEGER :: mpi_err

        ! Fetch local values first
        fc_loc = f_c_loc()

        ! Take the global minimum
        CALL MPI_ALLREDUCE( fc_loc, f_c_glb, 1, &
                            my_mpi_real, MPI_MIN, config%mpi_comm_cart(), mpi_err )

    END FUNCTION f_c_glb



    !--------------------------------------------------------------------------
    ! Calculates the CFL-number locally
    !--------------------------------------------------------------------------
    FUNCTION cfl_loc()
        USE configuration_mod, ONLY : config => configuration
        USE geometric_paras_mod, ONLY : phi, theta, r, r_sin_theta
        REAL(real_kind) :: cfl_loc, cfl(3)
        REAL(real_kind), ALLOCATABLE :: arc_lgt(:,:,:,:,:,:), v_max(:,:,:,:,:,:)

        ! Obtain maximum velocity by merging vs and vp velocities
        v_max = cp  ! Workaround concerning a gfortran 4.7.2 bug
        v_max(:,:,:,:,:,:) = MERGE( cp, cs, cs < cp )

        ASSOCIATE( lpd => config%lpd(), &
                   dt  => config%dt() )

        ! arc-length in theta-direction [m]
        arc_lgt = r(:,:,:,1:,:,:) * &
            ( theta(:,:,:,1:,:,:) - theta(:,:,:,:lpd-1,:,:) )
        cfl(1) = MAXVAL( v_max(:,:,:,1:,:,:) / arc_lgt ) * dt

        ! arc-length in phi-direction
        arc_lgt = r_sin_theta(:,:,:,:,1:,:) * &
            ( phi(:,:,:,:,1:,:) - phi(:,:,:,:,:lpd-1,:) )
        cfl(2) = MAXVAL( v_max(:,:,:,:,1:,:) / arc_lgt ) * dt

        ! Spacing in radial-direction
        arc_lgt = r(:,:,:,:,:,:lpd-1) - r(:,:,:,:,:,1:)
        cfl(3) = MAXVAL( v_max(:,:,:,:,:,1:) / arc_lgt ) * dt

        END ASSOCIATE

        ! Return maximum
        cfl_loc = MAXVAL( cfl )

    END FUNCTION cfl_loc

    !--------------------------------------------------------------------------
    ! Returns the global CFL-number
    ! Collective and blocking! Can only be called on all ranks simultaneously
    !--------------------------------------------------------------------------
    ! TODO check whether my_mpi_real and config%mpi_comm_cart() are initialized
    FUNCTION cfl_glb()
        USE configuration_mod, ONLY : config => configuration
#ifndef MPI_INCLUDE
        !! RECOMMENDED !!
        ! in case MPI_INCLUDE is NOT defined the mpi module will be used
        USE mpi
#else
        !! DEPRECATED !!
        ! when MPI_INCLUDE is defined, mpi header-files will be included
        INCLUDE 'mpif.h'
#endif
        REAL(real_kind) :: cfl_glb, courant_loc
        INTEGER :: mpi_err

        ! Fetch local values first
        courant_loc = cfl_loc()

        ! Return the global maximum
        CALL MPI_ALLREDUCE( courant_loc, cfl_glb, 1, &
            my_mpi_real, MPI_MAX, config%mpi_comm_cart(), mpi_err )

    END FUNCTION cfl_glb


END MODULE model_paras_mod
