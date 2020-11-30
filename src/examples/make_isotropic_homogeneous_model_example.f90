!> @example make_homogeneous_model_example.f90
!!
!! Example program to generate homogeneous models
!!
!! @todo  Describe what it does ...
!!
!! Compile simply using the Makefile:
!!
!!      make make make_homogeneous_model_example
!!
!! ...
PROGRAM make_isotropic_homogeneous_model_example
    USE parameters_mod, ONLY : real_kind, fnl
    USE parser_mod, ONLY : grid_typ, parse_namelist

    IMPLICIT NONE

    INTEGER :: lun
    TYPE(grid_typ) :: grid
    REAL(real_kind), ALLOCATABLE :: param(:,:,:,:,:,:)
    REAL(real_kind) :: rho, vp, vs
    CHARACTER(LEN=fnl) :: fn

    IF ( COMMAND_ARGUMENT_COUNT() /= 4 ) THEN
        WRITE(UNIT=*, FMT='(A)') 'Missing parameters!'
        WRITE(UNIT=*, FMT='(A)') 'Usuage: make_isotropic_homogeneous_model_&
            &example [rho] [vs] [vp] [config.-file]'
        STOP
    END IF

    CALL GET_COMMAND_ARGUMENT(1, fn)
    READ(fn, FMT=*) rho
    CALL GET_COMMAND_ARGUMENT(2, fn)
    READ(fn, FMT=*) vs
    CALL GET_COMMAND_ARGUMENT(3, fn)
    READ(fn, FMT=*) vp

    CALL GET_COMMAND_ARGUMENT(4, fn)
    OPEN(NEWUNIT=lun, FILE=fn)
        CALL parse_namelist(lun, grid)
    CLOSE(UNIT=lun)

    ALLOCATE( param(grid%nx,grid%ny,grid%nz,grid%lpd+1,grid%lpd+1,grid%lpd+1) )

    ! rhoinv
    param(:,:,:,:,:,:) = 1.0 / rho

    OPEN(NEWUNIT=lun, FILE='./rhoinv', ACCESS='stream')
    WRITE(UNIT=lun) param(:,:,:,:,:,:)
    CLOSE(UNIT=lun)

    ! mu
    param(:,:,:,:,:,:) = rho * vs**2

    OPEN(NEWUNIT=lun, FILE='./mu', ACCESS='stream')
    WRITE(UNIT=lun) param
    CLOSE(UNIT=lun)

    ! lambda
    param(:,:,:,:,:,:) = rho * ( vp**2 - 2 * vs**2 )

    OPEN(NEWUNIT=lun, FILE='./lambda', ACCESS='stream')
    WRITE(UNIT=lun) param
    CLOSE(UNIT=lun)

END PROGRAM


!> @file
!! $Date: 2013-11-26 21:39:48 +0100 (Tue, 26 Nov 2013) $
!! $Author: mauerberger $
!! $Revision: 817 $
!! @copyright GNU General Public License version 3 or later
