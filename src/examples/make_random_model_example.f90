!> @example make_random_model_example.f90
!!
!! Example program adding random perturbations on top of a homogeneous model. 
!! Inside each element parameters remain constant.  
!! 
!! @todo Describe what it does ...
!!
!! To compile simply use the Makefile: 
!!
!!     make make_random_model_example
!!
!! ...
PROGRAM make_homogeneous_model_example
    USE parameters_mod, ONLY : real_kind, fnl
    USE parser_mod, ONLY : grid_typ, parse_namelist

    IMPLICIT NONE

    INTEGER :: lun
    TYPE(grid_typ) :: grid
    REAL(real_kind), ALLOCATABLE :: param(:,:,:,:,:,:)
    REAL(real_kind), PARAMETER :: rho = 2800.0, &
                                  vp  = 6000.0, &
                                  vs  = 3000.0, &
                                  rate = 0.1
    CHARACTER(LEN=fnl) :: fn

    IF ( COMMAND_ARGUMENT_COUNT() /= 1 ) THEN
        WRITE(UNIT=*, FMT='(A)') 'Missing parameter!'
        WRITE(UNIT=*, FMT='(A)') 'Usuage: make_random_model_example [config.-file]'
        STOP
    END IF

    CALL GET_COMMAND_ARGUMENT(1, fn) 

    OPEN(NEWUNIT=lun, FILE=fn)
        CALL parse_namelist(lun, grid)
    CLOSE(UNIT=lun)

    ! Allocate memory 
    ALLOCATE( param(grid%nx,grid%ny,grid%nz,grid%lpd+1,grid%lpd+1,grid%lpd+1) )

    ! rhoinv
    !!!!!!!!
    
    ! Get field with randeom values in the range 1.0+-0.1
    CALL random_element_values(param, rate)
    param(:,:,:,:,:,:) = param / rho  

    OPEN(NEWUNIT=lun, FILE='./rhoinv', ACCESS='stream')
    WRITE(UNIT=lun) param(:,:,:,:,:,:)
    CLOSE(UNIT=lun)


    ! mu
    !!!!

    ! Get field with randeom values in the range 1.0+-0.1
    CALL random_element_values(param, rate)
    param(:,:,:,:,:,:) = param * rho * vs**2

    OPEN(NEWUNIT=lun, FILE='./mu', ACCESS='stream')
    WRITE(UNIT=lun) param
    CLOSE(UNIT=lun)


    ! lambda
    !!!!!!!!

    ! Get field with randeom values in the range 1.0+-0.1
    CALL random_element_values(param, rate)
    param(:,:,:,:,:,:) = param * rho * ( vp**2 - 2 * vs**2 )

    OPEN(NEWUNIT=lun, FILE='./lambda', ACCESS='stream')
    WRITE(UNIT=lun) param
    CLOSE(UNIT=lun)

CONTAINS 

    !> All values of 'field' - a whole 3d field - are filled with 1.0 plus 's' 
    !! percent of random values     
    SUBROUTINE random_element_values(field, s)
        ! Passed dummy arguments  
        REAL(real_kind), INTENT(IN) :: s
        REAL(real_kind), INTENT(OUT) :: field(:,:,:,:,:,:)
        ! Local variables
        REAL(real_kind) :: random(SIZE(field, dim=1), &
                                  SIZE(field, dim=2), &
                                  SIZE(field, dim=3))
        INTEGER :: lun, i, j, k, n
        INTEGER, ALLOCATABLE :: seed(:)

        ! Obtain size of random seed
        CALL RANDOM_SEED(SIZE=n)
        ALLOCATE(seed(n))
        ! Invoke the OS random number generator to get a random seed
        OPEN(NEWUNIT=lun, FILE="/dev/urandom", ACCESS="stream", &
             FORM="unformatted", ACTION="read", STATUS="old" )
        READ(UNIT=lun) seed
        CLOSE(UNIT=lun)

        ! Put a random seed
        CALL RANDOM_SEED(PUT=seed)

        ! Obtain random values for each element 
        CALL RANDOM_NUMBER(random)

        ! Scale values to the range 1.0 +- s
        random(:,:,:) = 1.0 + 2.0*s*(random - 0.5)

        ! Translate random values to the whole 3d field 
        ! covering all collocation points
        ! Inside each element the value remains the same 
        DO i=1, grid%nx
            DO j=1, grid%ny
                DO k=1, grid%nz
                    field(i,j,k,:,:,:) = random(i,j,k) 
                END DO
            END DO
        END DO

    END SUBROUTINE

END PROGRAM make_homogeneous_model_example


!> @file
!! $Date: 2013-10-30 16:25:25 +0100 (Wed, 30 Oct 2013) $
!! $Author: mauerberger $
!! $Revision: 743 $
!! @copyright GNU General Public License version 3 or later 
