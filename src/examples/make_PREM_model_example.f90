!> @example make_PREM_model_example.f90
!!
!! Example program to generate homogeneous models
!!
!! @todo Describe what it does ...
!!
!! Compile simply using the Makefile:
!!
!!      make make make_PREM_model_example
!!
!! ...
PROGRAM make_PREM_model_example
    USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY : OUTPUT_UNIT, ERROR_UNIT
    USE parameters_mod, ONLY : real_kind, fnl
    USE parser_mod, ONLY : grid_typ, model_typ, parse_namelist, default_character
    USE coordinates_mod, ONLY : r_coordinates
    USE prem_mod, ONLY : prem_iso

    IMPLICIT NONE

    INTEGER :: lun
    TYPE(grid_typ) :: grid
    TYPE(model_typ) :: model
    REAL(real_kind), ALLOCATABLE :: rho(:,:,:,:,:,:), &
                                    vs(:,:,:,:,:,:), &
                                    vp(:,:,:,:,:,:), &
                                    mu(:,:,:,:,:,:), &
                                    lambda(:,:,:,:,:,:), &
                                    r(:,:,:,:,:,:)
    CHARACTER(LEN=fnl) :: fn

    ! Check if a command line argument was passed
    IF ( COMMAND_ARGUMENT_COUNT() < 1 ) THEN
        WRITE(UNIT=ERROR_UNIT, FMT='(/,A)') 'Missing config file!'
        WRITE(UNIT=ERROR_UNIT, FMT='(A)') 'Usuage: make_PREM_model_exampe &
            &[config-file]'
        STOP 1
    END IF

    ! Get filename from 1st passed command line argument
    CALL GET_COMMAND_ARGUMENT(1, fn)

    ! Parse NAMELIST groups &grid and &model
    OPEN(NEWUNIT=lun, FILE=fn)
        CALL parse_namelist(lun, grid)
        CALL parse_namelist(lun, model)
    CLOSE(UNIT=lun)

    ! Check if filename are specified in conf.-file
    IF ( ANY( [model%rhoinv,model%lambda,model%mu] == default_character ) ) THEN
        WRITE(UNIT=ERROR_UNIT, FMT='(/,A)') 'Filenames for "rhoinv", "lambda" &
            &and "mu" not specified in &model-group!'
        STOP 1
    END IF

    !
    WRITE(UNIT=ERROR_UNIT, FMT='(/,"### Generating model PREM ###",/)')

    ! Allocate memory
    ! TODO: To save 1/3 of memory consumption calc. mu and lambda on the fly
    ALLOCATE(   rho(grid%nx,grid%ny,grid%nz,grid%lpd+1,grid%lpd+1,grid%lpd+1),&
                 vs(grid%nx,grid%ny,grid%nz,grid%lpd+1,grid%lpd+1,grid%lpd+1),&
                 vp(grid%nx,grid%ny,grid%nz,grid%lpd+1,grid%lpd+1,grid%lpd+1),&
                 mu(grid%nx,grid%ny,grid%nz,grid%lpd+1,grid%lpd+1,grid%lpd+1),&
             lambda(grid%nx,grid%ny,grid%nz,grid%lpd+1,grid%lpd+1,grid%lpd+1),&
                  r(grid%nx,grid%ny,grid%nz,grid%lpd+1,grid%lpd+1,grid%lpd+1) )

    ! ...
    r(:,:,:,:,:,:) = r_coordinates(elms_theta=grid%nx, &
                                   elms_phi=grid%ny, &
                                   elms_r=grid%nz, &
                                   lpd=grid%lpd, &
                                   min=model%rad_min, &
                                   max=model%rad_max )

    ! ...
    CALL prem_iso( R=r, RHO=rho, VS=vs, VP=vp, &
                   nocrust=.FALSE., noocean=.FALSE.)
    ! Transform into Ses3d's parametrization
    rho(:,:,:,:,:,:) = 1000.0 * rho
    vp(:,:,:,:,:,:) = 1000.0 * vp
    vs(:,:,:,:,:,:) = 1000.0 * vs
    mu(:,:,:,:,:,:) = rho * vs**2
    lambda(:,:,:,:,:,:) = rho * vp**2 - 2.0 * mu

    ! Write rhoinv
    OPEN(NEWUNIT=lun, FILE=TRIM(model%rhoinv), ACCESS='stream')
    WRITE(UNIT=OUTPUT_UNIT, FMT='("Writing rhoinv to ''", A, "''")') &
        TRIM(model%rhoinv)
    WRITE(UNIT=lun) 1.0/rho
    CLOSE(UNIT=lun)

    ! Write mu
    OPEN(NEWUNIT=lun, FILE=TRIM(model%mu), ACCESS='stream')
    WRITE(UNIT=OUTPUT_UNIT, FMT='("Writing mu to ''", A, "''")') &
        TRIM(model%mu)
    WRITE(UNIT=lun) mu
    CLOSE(UNIT=lun)

    ! Write lambda
    OPEN(NEWUNIT=lun, FILE=TRIM(model%lambda), ACCESS='stream')
    WRITE(UNIT=OUTPUT_UNIT, FMT='("Writing lambda to ''", A, "''")') &
        TRIM(model%lambda)
    WRITE(UNIT=lun) lambda
    CLOSE(UNIT=lun)

    WRITE(UNIT=OUTPUT_UNIT, FMT='()')

END PROGRAM make_PREM_model_example


!> @file
!! $Date: 2013-10-30 16:25:25 +0100 (Wed, 30 Oct 2013) $
!! $Author: mauerberger $
!! $Revision: 743 $
!! @copyright GNU General Public License version 3 or later
