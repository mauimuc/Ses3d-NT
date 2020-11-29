! SES3D - simulation of elastic wave propagation in spherical sections

! (c) by Andreas Fichtner <A.Fichtner@uu.nl>
!        Tarje Nissen-Meyer <tarjen@ethz.ch>
!        Heiner Igel <igel@geophysik.uni-muenchen.de>
!        Stefan Mauerberger <mauerberger@geophysik.uni-muenchen.de>

! This program is free software: you can redistribute it and/or modify
! under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.

! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.

! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

! Last changed: $Date: 2013-09-16 15:37:13 +0200 (Mon, 16 Sep 2013) $
! By: $Author: mauerberger $
! Revision: $Revision: 686 $

!------------------------------------------------------------------------------
! Module to compute Frechet Kernels
!------------------------------------------------------------------------------
MODULE grad_mod
    USE error_mod, ONLY : abort
    USE parameters_mod, ONLY : real_kind, fnl, mpi_master_rank, is_mpi_master_rank
    USE configuration_mod, ONLY : config => configuration

    IMPLICIT NONE

    PRIVATE

    TYPE :: gradient_cls
        INTEGER :: timestep_start, timestep_increment, timestep_end
        CHARACTER(LEN=fnl) :: prefix
        CHARACTER(LEN=8), ALLOCATABLE :: attributes(:)
    CONTAINS
        !PROCEDURE :: file_name
        PROCEDURE :: file_names
    END TYPE

    REAL(real_kind), ALLOCATABLE :: grad_rho(:,:,:,:,:,:), &
                                    grad_cs(:,:,:,:,:,:), &
                                    grad_cp(:,:,:,:,:,:), &
                                    grad_csh(:,:,:,:,:,:), &
                                    grad_csv(:,:,:,:,:,:)

    INTERFACE gradient_cls
        MODULE PROCEDURE parse_gradient_file
    END INTERFACE

    TYPE(gradient_cls) :: gradient

    INTEGER, ALLOCATABLE :: saving_vector(:)

    PUBLIC :: ses3d_grad, allocate_gradient_variables, &
              grad_rho, grad_cs, grad_cp, grad_csh, grad_csv, gradient_cls, gradient
CONTAINS

    FUNCTION file_names(self, channel)
        USE string_utilities_mod, ONLY : i2c
        CLASS(gradient_cls), INTENT(IN) :: self
        CHARACTER(LEN=*), INTENT(IN) :: channel
        CHARACTER(LEN=fnl) :: file_names( n_files(self) )
        INTEGER :: sampels( n_files(self) )
        INTEGER :: i

        sampels(:) = samples(self)

        DO i=1, n_files(self)
            file_names(i) = TRIM(self%prefix) // TRIM(channel) // '_' // &
                i2c(sampels(i)) // '.vol'
        END DO

    END FUNCTION file_names

    PURE FUNCTION samples(self)
        CLASS(gradient_cls), INTENT(IN) :: self
        INTEGER :: samples( n_files(self) )
        INTEGER :: k
        samples(:) = [ (k, k=self%timestep_start, self%timestep_end, &
                           self%timestep_increment) ]
    END FUNCTION samples

    PURE INTEGER FUNCTION n_files(self)
        CLASS(gradient_cls), INTENT(IN) :: self
        INTEGER :: k
        n_files = SIZE([ (k, k=self%timestep_start, self%timestep_end, &
                             self%timestep_increment) ])
    END FUNCTION n_files


    SUBROUTINE allocate_gradient_variables

        ASSOCIATE( nx  => config%nx_loc(), &
                   ny  => config%ny_loc(), &
                   nz  => config%nz_loc(), &
                   lpd => config%lpd() )

        ALLOCATE( grad_rho(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd), SOURCE=0.0_real_kind )
        ALLOCATE( grad_cs (0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd), SOURCE=0.0_real_kind )
        ALLOCATE( grad_cp (0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd), SOURCE=0.0_real_kind )
        ALLOCATE( grad_csv(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd), SOURCE=0.0_real_kind )
        ALLOCATE( grad_csh(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd), SOURCE=0.0_real_kind )

        END ASSOCIATE

    END SUBROUTINE allocate_gradient_variables

SUBROUTINE ses3d_grad( it )
    USE elastic_vars_mod, ONLY : vx, vy, vz, dxux, dyux, dzux, dxuy, dyuy, &
                                 dzuy, dxuz, dyuz, dzuz
    USE model_paras_mod, ONLY : cp, cs, rho
    USE geometric_paras_mod, ONLY : jac
    USE mpi_parallel_io_mod, ONLY : mpi_read_parallel

    CHARACTER(LEN=fnl) :: fn
    INTEGER, INTENT(IN) :: it
    INTEGER :: i
    REAL(real_kind), ALLOCATABLE :: vx_fw (:,:,:,:,:,:), &
                                    vy_fw (:,:,:,:,:,:), &
                                    vz_fw (:,:,:,:,:,:), &
                                    exx_fw(:,:,:,:,:,:), &
                                    eyy_fw(:,:,:,:,:,:), &
                                    ezz_fw(:,:,:,:,:,:), &
                                    exy_fw(:,:,:,:,:,:), &
                                    exz_fw(:,:,:,:,:,:), &
                                    eyz_fw(:,:,:,:,:,:), &
                                    div_u_fw(:,:,:,:,:,:), &
                                    div_u_rw(:,:,:,:,:,:), &
                                    e_ddot_e(:,:,:,:,:,:)

    ASSOCIATE( nt  => config%nt(), &
               dt  => config%dt(), &
               nx  => config%nx_loc(), &
               ny  => config%ny_loc(), &
               nz  => config%nz_loc(), &
               lpd => config%lpd() )

    ALLOCATE( vx_fw (0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd), &
              vy_fw (0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd), &
              vz_fw (0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd), &
              exx_fw(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd), &
              eyy_fw(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd), &
              ezz_fw(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd), &
              exy_fw(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd), &
              exz_fw(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd), &
              eyz_fw(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd), &
              div_u_fw(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd), &
              div_u_rw(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd), &
              e_ddot_e(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd) )

    !======================================================================
    ! Initialisation of saving vector to read and restore fields, needed
    ! for adjoint simulations
    !======================================================================
    IF ( .NOT. ALLOCATED(saving_vector) ) THEN
        saving_vector = [ ( i, i=1, nt ) ] ! lhs-(re)allocation F2003
        WHERE ( MOD(saving_vector,gradient%timestep_increment) /= 0 )
            saving_vector(:) = 0
        END WHERE
    END IF


    !======================================================================
    ! The sensitivity densities assume that all measurements are made on
    ! velocity fields !!!
    !======================================================================

    IF ( saving_vector(nt+1-it)>0 ) THEN

        !- read velocity fields ---------------------------------------
        WRITE( UNIT=fn, FMT='(A,A,"_",I0,".vol")') TRIM(gradient%prefix), 'vx', saving_vector(nt+1-it)
        vx_fw = mpi_read_parallel( fn )
        WRITE( UNIT=fn, FMT='(A,A,"_",I0,".vol")') TRIM(gradient%prefix), 'vy', saving_vector(nt+1-it)
        vy_fw = mpi_read_parallel( fn )
        WRITE( UNIT=fn, FMT='(A,A,"_",I0,".vol")') TRIM(gradient%prefix), 'vz', saving_vector(nt+1-it)
        vz_fw = mpi_read_parallel( fn )

        !- read strain rates ------------------------------------------
        !exx_fw =
        !eyy_fw =
        !ezz_fw =

        !exy_fw =
        !exz_fw =
        !eyz_fw =

        !==============================================================
        ! compute Frechet kernels
        !==============================================================

        div_u_fw(:,:,:,:,:,:) = exx_fw + eyy_fw + ezz_fw
        div_u_rw(:,:,:,:,:,:) = dxux + dyuy + dzuz

        e_ddot_e(:,:,:,:,:,:) = dxux*exx_fw + dyuy*eyy_fw + dzuz*ezz_fw + &
                                (dxuy + dyux)*exy_fw + &
                                (dxuz + dzux)*exz_fw + &
                                (dyuz + dzuy)*eyz_fw

         grad_cp(:,:,:,:,:,:) = grad_cp + gradient%timestep_increment*dt* &
                                (2.0*rho*cp*div_u_fw*div_u_rw)/Jac

         grad_cs(:,:,:,:,:,:) = grad_cs + gradient%timestep_increment*dt* &
                                (4.0*rho*cs*(e_ddot_e-div_u_rw*div_u_fw))/Jac

        ! XXX: What is the difference?
        grad_rho(:,:,:,:,:,:) = grad_rho + gradient%timestep_increment*dt* &
                                (vx_fw*vx + vy_fw*vy + vz_fw*vz)/Jac
        !grad_rho(:,:,:,:,:,:) = grad_rho + gradient%timestep_increment*dt* &
        !                        ( (cp**2 - 2.0*cs**2)*div_u_fw*div_u_rw + &
        !                          2.0*cs**2*e_ddot_e )/Jac

        grad_csh(:,:,:,:,:,:) = grad_csh + gradient%timestep_increment*dt* &
                                ( 4.0*rho*cs*( e_ddot_e - div_u_rw*div_u_fw - &
                                               (dxuz + dzux)*exz_fw - &
                                               (dyuz + dzuy)*eyz_fw ) )/jac

        grad_csv(:,:,:,:,:,:) = grad_csv + gradient%timestep_increment*dt*&
                                ( 4.0*rho*cs*( (dxuz + dzux)*exz_fw + &
                                               (dyuz + dzuy)*eyz_fw ) )/jac

    ENDIF

    END ASSOCIATE

END SUBROUTINE ses3d_grad

    !-----------------------------------------------------------------
    ! this function parse all occurring &receiver-groups in a file
    ! and returns an array of receiver_cls objects with several meta-
    ! data already specified
    !-----------------------------------------------------------------
    FUNCTION parse_gradient_file(gradient_file_name) RESULT(gradient)
        USE parser_mod, ONLY : count_namelists, parse_namelist, &
                               grad_typ
        USE bcast_mod, ONLY : bcast_namelist
        USE string_utilities_mod, ONLY : i2c
#ifndef MPI_INCLUDE
       !! RECOMMENDED !!
       ! in case MPI_INCLUDE is NOT defined the mpi module will be used
       USE mpi
#else
       !! DEPRECATED !!
       ! when MPI_INCLUDE is defined, mpi header-files will be included
       INCLUDE 'mpif.h'
#endif
        CHARACTER(LEN=*), INTENT(IN) :: gradient_file_name
        TYPE(gradient_cls) :: gradient
        INTEGER :: gradient_file_unit
        TYPE(grad_typ) :: grad
        INTEGER :: n
        CHARACTER(LEN=8), ALLOCATABLE :: attributes(:)

        ! Open gradient file
        IF ( is_mpi_master_rank() ) &
            OPEN( NEWUNIT=gradient_file_unit, FILE=gradient_file_name, &
                  ACTION='READ', STATUS='OLD' )

        !------------
        ! Parse grid
        !------------

        ! Count number of NAMELIST groups
        IF ( is_mpi_master_rank() ) THEN
            n = count_namelists( unit=gradient_file_unit, nml_cls=grad)
            IF ( n /= 1 ) &
                CALL abort( 'Found '//i2c(n)//' &grad-groups in config-file!' )
            ! Parse NAMELIST
            REWIND( UNIT=gradient_file_unit )
            CALL parse_namelist(unit=gradient_file_unit, nml_cls=grad)
            ! Consistency checks
            IF ( 1 > SIZE(PACK(grad%attributes,grad%attributes /= ' ')) )&
                CALL abort('No gradients specified')

            ! all done -> close file
            CLOSE( gradient_file_unit )
        END IF
        ! Broadcast values
        CALL bcast_namelist(grad, mpi_master_rank)
        ! Assign values
        attributes = PACK(grad%attributes,grad%attributes /= ' ')
        gradient = gradient_cls(timestep_start=grad%timestep_start, &
                                timestep_increment=grad%timestep_increment, &
                                timestep_end=grad%timestep_end, &
                                prefix='./', &
                                attributes=attributes(:) )

    END FUNCTION parse_gradient_file

END MODULE grad_mod
