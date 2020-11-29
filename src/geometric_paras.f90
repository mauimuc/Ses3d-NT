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



!> Geometric parameters
!!
!! ...
MODULE geometric_paras_mod
    USE parameters_mod, ONLY : real_kind

    IMPLICIT NONE

    PRIVATE

    ! For elements Jacobian
    REAL(real_kind) :: jac

    ! Local memory variables
    ! TODO should be somewhere in parameters_mod maybe as hidden parameter
    INTEGER, PARAMETER :: nrdiss=3  ! Number of relaxation mechanisms
    REAL(real_kind) :: tau_p(1:nrdiss) = [ 1000.00, 70.00, 5.25 ]

    ! Coordinate fields
    REAL(real_kind), PROTECTED, ALLOCATABLE :: theta(:,:,:,:,:,:), &
                                               phi(:,:,:,:,:,:), &
                                               r(:,:,:,:,:,:)
    ! Auxiliary coordinate fields
    REAL(real_kind), PROTECTED, ALLOCATABLE :: sin_theta(:,:,:,:,:,:), &
                                               cot_theta(:,:,:,:,:,:), &
                                               r_sin_theta(:,:,:,:,:,:), &
                                               rr_sin_theta(:,:,:,:,:,:), &
                                               wx_wy_wz(:,:,:,:,:,:)
    ! Relaxing boundaries
    REAL(real_kind), PROTECTED, ALLOCATABLE :: prof_x(:,:,:,:,:,:), &
                                               prof_y(:,:,:,:,:,:), &
                                               prof_z(:,:,:,:,:,:), &
                                               prof(:,:,:,:,:,:), &
                                               taper(:,:,:,:,:,:)

    PUBLIC :: nrdiss, jac, tau_p, &
              theta, phi, r, &
              sin_theta, cot_theta, r_sin_theta, rr_sin_theta, &
              wx_wy_wz, &
              prof_x, prof_y, prof_z, prof, taper, &
              allocate_geometric_parameters, init_geometric_parameters


CONTAINS

    !--------------------------------------------------------------------------
    ! Allocates memory for geometric parameters
    !--------------------------------------------------------------------------
    SUBROUTINE allocate_geometric_parameters
        USE configuration_mod, ONLY : config => configuration

        ! Aliases just to keep expressions concise
        ASSOCIATE( nx  => config%nx_loc(), &
                   ny  => config%ny_loc(), &
                   nz  => config%nz_loc(), &
                   lpd => config%lpd() )

        ! Allocate 3D fields for auxiliary coordinate values
        ALLOCATE(       theta(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd),source=0.0_real_kind)
        ALLOCATE(         phi(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd),source=0.0_real_kind)
        ALLOCATE(           r(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd),source=0.0_real_kind)
        ALLOCATE(   sin_theta(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd),source=0.0_real_kind)
        ALLOCATE(   cot_theta(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd),source=0.0_real_kind)
        ALLOCATE( r_sin_theta(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd),source=0.0_real_kind)
        ALLOCATE(rr_sin_theta(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd),source=0.0_real_kind)
        ! Allocate auxiliary field for GLL integration
        ALLOCATE(wx_wy_wz(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd),source=0.0_real_kind)
        ! Allocate fields for relaxing boundaries
        ALLOCATE(prof_x(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd),source=0.0_real_kind)
        ALLOCATE(prof_y(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd),source=0.0_real_kind)
        ALLOCATE(prof_z(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd),source=0.0_real_kind)
        ALLOCATE(  prof(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd),source=0.0_real_kind)
        ALLOCATE( taper(0:nx,0:ny,0:nz,0:lpd,0:lpd,0:lpd),source=0.0_real_kind)

        END ASSOCIATE

    END SUBROUTINE allocate_geometric_parameters


    !--------------------------------------------------------------------------
    ! Initializes geometric parameters
    !--------------------------------------------------------------------------
    SUBROUTINE init_geometric_parameters
        USE configuration_mod, only : config => configuration
        USE gll_mod, ONLY : get_weights
        USE coordinates_mod, ONLY : theta_coordinates, phi_coordinates, &
                                    r_coordinates
        ! Just local indices
        INTEGER :: i, j, k, n
        ! Local variables (indicated by trailing underscores)
        REAL(real_kind), ALLOCATABLE :: delta_x_(:,:), delta_y_(:,:), &
                                        delta_z_(:,:), x_(:,:), y_(:,:), &
                                        z_(:,:), weights_(:)

        ! Just some aliases to keep expressions concise
        ASSOCIATE( alpha  => config%taper_slope, &
                   pml    => config%taper_width, &
                   nx_loc => config%nx_loc(), &
                   ny_loc => config%ny_loc(), &
                   nz_loc => config%nz_loc(), &
                   lpd    => config%lpd() )

        ! Make elements Jacobi determinant
        jac = config%dx_elm() * config%dy_elm() * config%dz_elm() / 8.0


        ! Temporarily allocate memory for 1d coordinate fields
        ALLOCATE( x_(0:nx_loc,0:lpd))
        ALLOCATE( y_(0:ny_loc,0:lpd))
        ALLOCATE( z_(0:nz_loc,0:lpd))

        ! Get x-coordinate axis
        x_(:,:) = theta_coordinates(elms=nx_loc+1, lpd=lpd, min=config%x_loc_min(), max=config%x_loc_max())
        ! Translate x-coordinates to 3d field
        DO i=0, nx_loc
            DO j=0, lpd
                theta(i,:,:,j,:,:) = x_(i,j)
            END DO
        END DO

        ! Get y-coordinate axis
        y_(:,:) = theta_coordinates(elms=ny_loc+1, lpd=lpd, min=config%y_loc_min(), max=config%y_loc_max())
        ! Translate y-coordinates to 3d field
        DO i=0, ny_loc
            DO j=0, lpd
                phi(:,i,:,:,j,:) = y_(i,j)
            END DO
        END DO

        ! Get z-coordinate axis
        ! In Ses3d-NT r coordinates are upside down. Thus, min and max have to be interchanged
        z_(:,:) = theta_coordinates(elms=nz_loc+1, lpd=lpd, min=config%z_loc_max(), max=config%z_loc_min())
        ! Translate z-coordinates to 3d field
        DO i=0, nz_loc
            DO j=0, lpd
                r(:,:,i,:,:,j) = z_(i,j)
            END DO
        END DO

        ! Due to optimization make auxiliary coordinates
        sin_theta(:,:,:,:,:,:) = SIN( theta )
        cot_theta(:,:,:,:,:,:) = COS( theta ) / SIN( theta )
        r_sin_theta(:,:,:,:,:,:) = r * sin_theta
        rr_sin_theta(:,:,:,:,:,:) = r * r_sin_theta


        ! Temporarily allocate memory for GLL integration weights
        ALLOCATE( weights_(0:lpd) )
        ! Fetch weights for GLL-integration
        weights_(:) = get_weights( lpd )

        ! Due to optimization make 3d auxiliary weights
        DO i=0, lpd ! x-direction
            DO j=0, lpd ! y-direction
                DO k=0, lpd ! z-direction
                    wx_wy_wz(:,:,:,i,j,k) = weights_(i)*weights_(j)*weights_(k)
                END DO
            END DO
        END DO

        ! Clean up
        DEALLOCATE( weights_ )


        ! x-direction damping profiles
        prof_x(:,:,:,:,:,:) = 0.0

        ! Left x-boundary damping profile
        IF ( config%mpi_north() < 0 ) THEN
            ! Temporarily allocate memory
            ALLOCATE( delta_x_(0:pml-1,0:lpd) )
            ! Calculate 1d damping profile
            delta_x_(:,:) = config%x_loc_min() + pml*config%dx_elm() - &
                            x_(0:pml-1,0:lpd)
            delta_x_(:,:) = alpha*( delta_x_(:,:)/(pml*config%dx_elm()) )**2
            ! Translate to 3d fields
            DO i=0,pml-1
                DO n=0,lpd
                    prof_x(i,:,:,n,:,:) = delta_x_(i,n)
                ENDDO
            ENDDO
            ! Clean up for latter use
            DEALLOCATE( delta_x_ )
        ENDIF

        ! Right x-boundary damping profile
        IF ( config%mpi_south() < 0 ) THEN
            ! Temporarily allocate memory
            ALLOCATE( delta_x_(nx_loc-pml+1:nx_loc,0:lpd) )
            ! Calculate 1d damping profile
            delta_x_(:,:) = config%x_loc_max() - pml*config%dx_elm() - &
                            x_(nx_loc-pml+1:,:)
            delta_x_(:,:) = alpha*( delta_x_(:,:)/(pml*config%dx_elm()) )**2
            ! Translate to 3d fields
            DO i=config%nx_loc()-pml+1, config%nx_loc()
                DO n=0, lpd
                    prof_x(i,:,:,n,:,:) =  delta_x_(i,n)
                ENDDO
            ENDDO
            ! Clean up
            DEALLOCATE( delta_x_ )
        ENDIF


        ! y-direction damping profiles
        prof_y(:,:,:,:,:,:) = 0.0

        ! Left y-boundary damping profile
        IF ( config%mpi_west() < 0 ) THEN
            ! Temporarily allocate memory
            ALLOCATE( delta_y_(0:pml-1,0:lpd) )
            ! Calculate 1d damping profile
            delta_y_(:,:) = config%y_loc_min() + pml*config%dy_elm() - &
                            y_(0:pml-1,0:lpd)
            delta_y_(:,:) = alpha*( delta_y_(:,:)/(pml*config%dy_elm()) )**2
            ! Translate to 3d fields
            DO i=0, pml-1
                DO n=0, lpd
                    prof_y(:,i,:,:,n,:) = delta_y_(i,n)
                ENDDO
            ENDDO
            ! Clean up for latter use
            DEALLOCATE( delta_y_ )
        ENDIF

        ! Right y-boundary damping profile
        IF ( config%mpi_east() < 0) THEN
            ! Temporarily allocate memory
            ALLOCATE( delta_y_(ny_loc-pml+1:ny_loc,0:lpd) )
            ! Calculate 1d damping profile
            delta_y_(:,:) = config%y_loc_max() - pml*config%dy_elm() - &
                            y_(ny_loc-pml+1:,:)
            delta_y_(:,:) = alpha*( delta_y_(:,:)/(pml*config%dy_elm()) )**2
            ! Translate to 3d fields
            DO i=config%ny_loc()-pml+1, config%ny_loc()
                DO n=0, lpd
                    prof_y(:,i,:,:,n,:) = delta_y_(i,n)
                ENDDO
            ENDDO
            ! Clean up
            DEALLOCATE( delta_y_ )
        ENDIF


        ! Upper z-boundary damping profile
        ! Why the hell are z-coordinates upside down?
        prof_z(:,:,:,:,:,:) = 0.0
        IF ( config%mpi_bottom() < 0 ) THEN
            ! Free surface; Nothing to do
        ENDIF

        ! Lower z-boundary damping profile
        ! Why the hell are z-coordinates upside down?
        IF ( config%mpi_top() < 0 ) THEN
            ! Temporarily allocate memory
            ALLOCATE( delta_z_(nz_loc-pml+1:nz_loc,0:lpd) )
            ! Calculate 1d damping profile
            delta_z_(:,:) = config%z_loc_min() + pml*config%dz_elm() - &
                            z_(nz_loc-pml+1:,:)
            delta_z_(:,:) = alpha*( delta_z_(:,:)/(pml*config%dz_elm()) )**2
            ! Translate to 3d fields
            DO k=nz_loc-pml+1, nz_loc
                DO n=0, lpd
                    prof_z(:,:,k,:,:,n) = delta_z_(k,n)
                ENDDO
            ENDDO
            ! Clean up
            DEALLOCATE( delta_z_ )
        ENDIF


        WHERE( theta > theta(pml-1,0,0,lpd,0,0) .AND. &
               theta < theta(nx_loc-pml+1,0,0,0,0,0) )
            prof(:,:,:,:,:,:) = 0.0
        ELSE WHERE
            prof(:,:,:,:,:,:) = ABS(theta - (MAXVAL(theta) + MINVAL(theta))/2.0)
            prof(:,:,:,:,:,:) = prof - MINVAL(PACK(prof,prof/=0.0))
            prof(:,:,:,:,:,:) = alpha*(prof/MAXVAL(prof))**2
        END WHERE


        ! Normalise corners
        prof(:,:,:,:,:,:) = prof_x + prof_y + prof_z

        prof_z(:,:,:,:,:,:) = 2.0*prof_z/(1.0+prof/alpha)
        prof_y(:,:,:,:,:,:) = 2.0*prof_y/(1.0+prof/alpha)
        prof_x(:,:,:,:,:,:) = 2.0*prof_x/(1.0+prof/alpha)

        prof(:,:,:,:,:,:) = prof_x + prof_y + prof_z

        taper(:,:,:,:,:,:) = EXP( -0.1 * prof**2 )

        ! Clean up
        DEALLOCATE( x_ )
        DEALLOCATE( y_ )
        DEALLOCATE( z_ )

        END ASSOCIATE

    END SUBROUTINE init_geometric_parameters


    !> Returns a Gaussian taper used as absorbing boundary
    !!
    !! The idea is to multiply the wave field inside a narrow strip along the
    !! artificial boundaries by a Gaussian taper
    !! \f[
    !!   t(\mathbf x) \approx \exp^{-\alpha x^2} \;\;\;
    !!      \text{for} \; x \in \partial\Omega
    !! \f]
    !! The actual implementation works as follows:
    !!  1. Find the precise coordinate values where the taper shall act
    !!  2. Make profiles for all theta, phi and r directions
    !!  3. Normalize corners
    !!
    !! @param width Width in elements of the absorbing boundary
    !! @param alpha Strength ...
    !! @result A rank six field ...
    FUNCTION init_taper(width, alpha) RESULT(t)
        USE coordinates_mod, ONLY : theta_coordinates, phi_coordinates, &
                                    r_coordinates
        USE configuration_mod, ONLY : conf => configuration
        INTEGER, INTENT(IN) :: width
        REAL(real_kind), INTENT(IN) :: alpha
        REAL(real_kind) :: t!(0:conf%nx_glb(),0:conf%ny_glb(),0:conf%nz_glb(), &
                            ! 0:conf%lpd(),0:conf%lpd(),0:conf%lpd())
        REAL(real_kind) :: taper_theta_min, taper_theta_max, &
                           taper_phi_min, taper_phi_max, taper_r_min
        REAL(real_kind) :: theta_glb(0:conf%nx_glb(),0:conf%lpd()), &
                           phi_glb(0:conf%ny_glb(),0:conf%lpd()), &
                           r_glb(0:conf%nz_glb(),0:conf%lpd())

        ! Make global coordinate vectors
        ! Global coordinates are needed here just to figure out at which
        ! coordinate value tapering boundaries start and stop
        theta_glb(:,:) = theta_coordinates( elms=conf%nx_glb(), &
                                            lpd=conf%lpd()+1, &
                                            min=conf%x_glb_min(), &
                                            max=conf%x_glb_max() )
        phi_glb(:,:) = phi_coordinates( elms=conf%ny_glb(), &
                                        lpd=conf%lpd()+1, &
                                        min=conf%y_glb_min(), &
                                        max=conf%y_glb_max() )
        r_glb(:,:) = r_coordinates( elms=conf%nz_glb(), &
                                    lpd=conf%lpd()+1, &
                                    min=conf%z_glb_min(), &
                                    max=conf%z_glb_max() )

        ! Find at which coordinate value the tapering boundaries stops
        ! The width of the taper is defined in elements
        !taper_theta_min = theta_glb(conf%pml()-1,conf%lpd())
        !taper_theta_max = theta_glb(conf%nx_glb()-conf%pml()+1,0)
        !taper_phi_min = phi_glb(conf%pml()-1,conf%lpd())
        !taper_phi_max = phi_glb(conf%ny_glb()-conf%pml()+1,0)
        !taper_r_min = r_glb(conf%pml()-1,conf%lpd())

        ! Make damping profiles for each direction
        !WHERE( theta > theta(pml-1,0,0,lpd,0,0) .AND. &
        !       theta < theta(nx_loc-pml+1,0,0,0,0,0) )
        !    prof(:,:,:,:,:,:) = 0.0
        !ELSE WHERE
        !    prof(:,:,:,:,:,:) = ABS(theta - (MAXVAL(theta) + MINVAL(theta))/2.0)
        !    prof(:,:,:,:,:,:) = prof - MINVAL(PACK(prof,prof/=0.0))
        !    prof(:,:,:,:,:,:) = alpha*(prof/MAXVAL(prof))**2
        !END WHERE

        ! Normalize corners

        ! Make actual taper
        taper(:,:,:,:,:,:) = exp( -alpha * taper**2 )

    END FUNCTION


END MODULE geometric_paras_mod
