! SES3D - simulation of elastic wave propagation in spherical sections

! (c) by Stefan Mauerberger

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

! Last changed: $Date: 2013-02-20 18:02:56 +0100 (Wed, 20 Feb 2013) $
! By: $Author: mauerberger $
! Revision: $Revision: 353 $


!------------------------------------------------------------------------------
! ELASTIC VARIABLES
!
! ...
!------------------------------------------------------------------------------
MODULE elastic_vars_mod
    USE parameters_mod

    IMPLICIT NONE

    PRIVATE

    ! Displacement velocity fields
    REAL(real_kind), ALLOCATABLE :: vx(:,:,:,:,:,:), &
                                    vy(:,:,:,:,:,:), &
                                    vz(:,:,:,:,:,:)
    ! Weak form stress
    REAL(real_kind), ALLOCATABLE :: sx(:,:,:,:,:,:), &
                                    sy(:,:,:,:,:,:), &
                                    sz(:,:,:,:,:,:)
    ! Strain fields
    REAL(real_kind), ALLOCATABLE :: exx(:,:,:,:,:,:), &
                                    exy(:,:,:,:,:,:), &
                                    exz(:,:,:,:,:,:), &
                                    eyx(:,:,:,:,:,:), &
                                    eyy(:,:,:,:,:,:), &
                                    eyz(:,:,:,:,:,:), &
                                    ezx(:,:,:,:,:,:), &
                                    ezy(:,:,:,:,:,:), &
                                    ezz(:,:,:,:,:,:)
    ! Strong form stress fields
    REAL(real_kind), ALLOCATABLE :: sxx_pml(:,:,:,:,:,:), &
                                    sxy_pml(:,:,:,:,:,:), &
                                    sxz_pml(:,:,:,:,:,:), &
                                    syy_pml(:,:,:,:,:,:), &
                                    syz_pml(:,:,:,:,:,:), &
                                    szz_pml(:,:,:,:,:,:)
    ! ...
    REAL(real_kind), ALLOCATABLE :: dxux(:,:,:,:,:,:), &
                                    dxuy(:,:,:,:,:,:), &
                                    dxuz(:,:,:,:,:,:), &
                                    dyux(:,:,:,:,:,:), &
                                    dyuy(:,:,:,:,:,:), &
                                    dyuz(:,:,:,:,:,:), &
                                    dzux(:,:,:,:,:,:), &
                                    dzuy(:,:,:,:,:,:), &
                                    dzuz(:,:,:,:,:,:)
    ! ...
    REAL(real_kind), ALLOCATABLE :: Mxx(:,:,:,:,:,:,:), &
                                    Mxy(:,:,:,:,:,:,:), &
                                    Mxz(:,:,:,:,:,:,:), &
                                    Myy(:,:,:,:,:,:,:), &
                                    Myz(:,:,:,:,:,:,:), &
                                    Mzz(:,:,:,:,:,:,:)

    ! Source fields
    REAL(real_kind), ALLOCATABLE :: src_xx(:,:,:,:,:,:), &
                                    src_xy(:,:,:,:,:,:), &
                                    src_xz(:,:,:,:,:,:), &
                                    src_yy(:,:,:,:,:,:), &
                                    src_yz(:,:,:,:,:,:), &
                                    src_zz(:,:,:,:,:,:), &
                                    src_x(:,:,:,:,:,:), &
                                    src_y(:,:,:,:,:,:), &
                                    src_z(:,:,:,:,:,:)

    PUBLIC :: vx, vy, vz, sx, sy, sz, exx, exy, exz, eyx, eyy, eyz, ezx, ezy, &
              ezz, sxx_pml, sxy_pml, sxz_pml, syy_pml, syz_pml, szz_pml, dxux,&
              dxuy, dxuz, dyux, dyuy, dyuz, dzux, dzuy, dzuz, Mxx, Mxy, Mxz, &
              Myy, Myz, Mzz, src_xx, src_xy, src_xz, src_yy, src_yz, src_zz, &
              src_x, src_y, src_z, &
              allocate_elastic_variables, v_ptp_loc, v_ptp_glb, kinetic_energy

CONTAINS

    !--------------------------------------------------------------------------
    !
    !--------------------------------------------------------------------------
    SUBROUTINE allocate_elastic_variables
        USE configuration_mod, ONLY : config => configuration
        USE geometric_paras_mod, only : nrdiss

        ! Just aliases to keep expressions concise
        ASSOCIATE( nx  => config%nx_loc(), &
                   ny  => config%ny_loc(), &
                   nz  => config%nz_loc(), &
                   lpd => config%lpd() )

        ! Allocate fields for displacement-velocity
        ALLOCATE( vx(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd), SOURCE=0.0_real_kind )
        ALLOCATE( vy(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd), SOURCE=0.0_real_kind )
        ALLOCATE( vz(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd), SOURCE=0.0_real_kind )

        ! Allocate fields for displacement-velocity
        ALLOCATE( sx(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd), SOURCE=0.0_real_kind )
        ALLOCATE( sy(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd), SOURCE=0.0_real_kind )
        ALLOCATE( sz(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd), SOURCE=0.0_real_kind )

        ! Allocate strain fields
        ALLOCATE( exx(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd), SOURCE=0.0_real_kind )
        ALLOCATE( exy(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd), SOURCE=0.0_real_kind )
        ALLOCATE( exz(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd), SOURCE=0.0_real_kind )
        ALLOCATE( eyx(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd), SOURCE=0.0_real_kind )
        ALLOCATE( eyy(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd), SOURCE=0.0_real_kind )
        ALLOCATE( eyz(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd), SOURCE=0.0_real_kind )
        ALLOCATE( ezx(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd), SOURCE=0.0_real_kind )
        ALLOCATE( ezy(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd), SOURCE=0.0_real_kind )
        ALLOCATE( ezz(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd), SOURCE=0.0_real_kind )

        ! Allocate stress fields
        ALLOCATE( sxx_pml(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd), source=0.0_real_kind )
        ALLOCATE( sxy_pml(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd), source=0.0_real_kind )
        ALLOCATE( sxz_pml(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd), source=0.0_real_kind )
        ALLOCATE( syy_pml(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd), source=0.0_real_kind )
        ALLOCATE( syz_pml(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd), source=0.0_real_kind )
        ALLOCATE( szz_pml(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd), source=0.0_real_kind )

        ! Allocate strain-rates
        ALLOCATE( dxux(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd), SOURCE=0.0_real_kind )
        ALLOCATE( dxuy(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd), SOURCE=0.0_real_kind )
        ALLOCATE( dxuz(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd), SOURCE=0.0_real_kind )
        ALLOCATE( dyux(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd), SOURCE=0.0_real_kind )
        ALLOCATE( dyuy(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd), SOURCE=0.0_real_kind )
        ALLOCATE( dyuz(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd), SOURCE=0.0_real_kind )
        ALLOCATE( dzux(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd), SOURCE=0.0_real_kind )
        ALLOCATE( dzuy(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd), SOURCE=0.0_real_kind )
        ALLOCATE( dzuz(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd), SOURCE=0.0_real_kind )

        ! Allocate memory-variables for visco-elasticity
        ALLOCATE( Mxx(1:nrdiss,0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd), SOURCE=0.0_real_kind )
        ALLOCATE( Mxy(1:nrdiss,0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd), SOURCE=0.0_real_kind )
        ALLOCATE( Mxz(1:nrdiss,0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd), SOURCE=0.0_real_kind )
        ALLOCATE( Myy(1:nrdiss,0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd), SOURCE=0.0_real_kind )
        ALLOCATE( Myz(1:nrdiss,0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd), SOURCE=0.0_real_kind )
        ALLOCATE( Mzz(1:nrdiss,0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd), SOURCE=0.0_real_kind )

        ! Allocate moment-tensor source fields
        ALLOCATE( src_xx(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd), source=0.0_real_kind )
        ALLOCATE( src_yy(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd), source=0.0_real_kind )
        ALLOCATE( src_zz(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd), source=0.0_real_kind )
        ALLOCATE( src_xy(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd), source=0.0_real_kind )
        ALLOCATE( src_xz(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd), source=0.0_real_kind )
        ALLOCATE( src_yz(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd), source=0.0_real_kind )
        ! Allocate monopole-force fields
        ALLOCATE( src_x(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd), source=0.0_real_kind )
        ALLOCATE( src_y(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd), source=0.0_real_kind )
        ALLOCATE( src_z(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd), source=0.0_real_kind )

        END ASSOCIATE

    END SUBROUTINE allocate_elastic_variables


    !--------------------------------------------------------------------------
    ! Returns local peak-to-peak displacement-velocity values
    !--------------------------------------------------------------------------
    PURE FUNCTION v_ptp_loc()
        REAL(real_kind) :: v_ptp_loc(2)

        v_ptp_loc(1) = MIN( MINVAL(vx), MINVAL(vy), MINVAL(vz) )
        v_ptp_loc(2) = MAX( MAXVAL(vx), MAXVAL(vy), MAXVAL(vz) )

    END FUNCTION v_ptp_loc

    !--------------------------------------------------------------------------
    ! Returns global peak-to-peak displacement-velocity values
    !--------------------------------------------------------------------------
    FUNCTION v_ptp_glb()
        USE configuration_mod, ONLY : config => configuration
! line control through pre-processor directives
#ifndef MPI_INCLUDE
        !! RECOMMENDED !!
        ! in case MPI_INCLUDE is NOT defined the mpi module will be used
        USE mpi!, only : mpi_allreduce, mpi_max, mpi_min
#else
        !! DEPRECATED !!
        ! when MPI_INCLUDE is defined, mpi header-files will be included
        INCLUDE 'mpif.h'
#endif

        REAL(real_kind) :: v_ptp_glb(2), loc(2)
        INTEGER :: mpi_err

        loc(:) = v_ptp_loc()

        CALL MPI_ALLREDUCE( loc(1), v_ptp_glb(1), 1, my_mpi_real, &
                            MPI_MIN, config%mpi_comm_cart(), mpi_err)
        CALL MPI_ALLREDUCE( loc(2), v_ptp_glb(2), 1, my_mpi_real, &
                            MPI_MAX, config%mpi_comm_cart(), mpi_err)

    END FUNCTION v_ptp_glb


    !--------------------------------------------------------------------------
    ! Compute kinetic energy
    !--------------------------------------------------------------------------
    FUNCTION kinetic_energy() RESULT( E_kin_glb )
        USE configuration_mod, ONLY : config => configuration
        USE geometric_paras_mod, ONLY : rr_sin_theta, wx_wy_wz, jac
        USE model_paras_mod, ONLY : rho
! line control through pre-processor directives
#ifndef MPI_INCLUDE
        !! RECOMMENDED !!
        ! in case MPI_INCLUDE is NOT defined the mpi module will be used
        USE mpi!, ONLY : MPI_ALLREDUCE, MPI_SUM
#else
        !! DEPRECATED !!
        ! when MPI_INCLUDE is defined, mpi header-files will be included
        INCLUDE 'mpif.h'
#endif

        REAL(real_kind) :: E_kin_loc, E_kin_glb
        INTEGER :: mpi_ierr

        ! compute kinetic energy locally
        ! STILL NOT EFFICIENT! r**2 * sin_theta * rho * wx_wy_wz is CONSTANT
        E_kin_loc = 0.5 * jac * &
            SUM( ( vx**2 + vy**2 + vz**2 ) * rr_sin_theta * rho * wx_wy_wz )

        ! mpi sum up to get global kinetic energy
        CALL MPI_ALLREDUCE( E_kin_loc, E_kin_glb, 1, my_mpi_real, &
                            MPI_SUM, config%mpi_comm_cart(), mpi_ierr )

    END FUNCTION kinetic_energy


END MODULE elastic_vars_mod
