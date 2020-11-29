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
!! $Date: 2013-11-25 20:20:46 +0100 (Mon, 25 Nov 2013) $
!! $Author: mauerberger $
!! $Revision: 814 $
!! @copyright GNU General Public License version 3 or later 



!> This module is dedicated to generate Ses3d-NT's coordinate fields
!!
!! A spherical coordinate system is used with \f$ \theta \in [0,\pi]\f$, 
!! \f$ \phi \in [-2\pi:2\pi] \f$ and \f$ r \in [0:r_{earth}]\f$ . 
!! 
!! Discretization is based on a regular, curvilinear grid approach. It splits 
!! into two parts:
!!  1. The computational domain becomes divided into elements. Each coordinate 
!!     direction is cut into \f$n_{theta}\f$, \f$n_{phi}\f$, \f$n_r\f$ equally 
!!     spaced pieces.
!!  2. Each element will be sampled at the so called N+1 GLL collocation points 
!!     (see gll_mod). 
!!
!! @todo Tests coordinates_test.f90
!!
!! This module expects to find a parameter parameters_mod::real_kind in
!! a module named parameters_mod which defines the accuracy of floating point 
!! numbers. <br>
!! In addition it relies on a function gll_mod::get_knots(N) - returning the 
!! N+1 GLL points - which has to be provided by a module named gll_mod.
MODULE coordinates_mod
    USE parameters_mod, ONLY : real_kind

    IMPLICIT NONE
    
    PRIVATE

    !> Generic Interface to distinguish between 1d and 3d case. In the first 
    !! case a rank 2 array is returned in the latter an array of rank 6.
    !! See coordinates_1d() and theta_coordinates_3d() for details.
    INTERFACE theta_coordinates
        MODULE PROCEDURE coordinates_1d
        MODULE PROCEDURE theta_coordinates_3d
    END INTERFACE theta_coordinates 

    !> Generic Interface to distinguish between 1d and 3d case. In the first 
    !! case a rank 2 array is returned in the latter an array of rank 6.
    !! See coordinates_1d() and theta_coordinates_3d() for details.
    INTERFACE phi_coordinates
        MODULE PROCEDURE coordinates_1d
        MODULE PROCEDURE phi_coordinates_3d
    END INTERFACE phi_coordinates 

    !> Generic Interface to distinguish between 1d and 3d case. In the first 
    !! case a rank 2 array is returned in the latter an array of rank 6.
    !! See coordinates_1d() and theta_coordinates_3d() for details.
    !!
    !! @warning In Ses3d-NT r coordinates are upside down. Thus, min and max 
    !!          have to be interchanged 
    INTERFACE r_coordinates
        MODULE PROCEDURE coordinates_1d
        MODULE PROCEDURE r_coordinates_3d
    END INTERFACE r_coordinates 

    PUBLIC :: theta_coordinates, phi_coordinates, r_coordinates, &
              pack_r2tor1

CONTAINS

    !> Returns a 1d coordinate field, an array of rank 2
    !! 
    !! The coordinate values \f$ x_{ij} \f$ of the returned array are covering 
    !! the interval starting from the passed argument min to max 
    !! \f$ [x_{min} : x_{max}] \f$. That are the lower and upper bounds. <br>
    !! The first index \f$ i \f$ is counting the elements. It ranges from 1 to
    !! the passed argument elms, \f$ i \in \{1 \dots N\} \f$. The element width
    !! \c delta is equally spaced \f$ \Delta = (x_{max} - x_{min})/N \f$. <br>
    !! The second index \f$ j \f$ is counting the actual collocation points 
    !! inside of an element \f$ j \f$. It ranges from 0 to the passed argument 
    !! lpd, \f$ j \in \{0 \dots lpd\} \f$. The points inside of an element are 
    !! arranged acording to the GLL collocation points. For details see
    !! gll_mod::get_knots(). <br>
    !! As in all FEM approaches, this means that the upper and lower coordinate
    !! values in adjacent elements are identical \f$ x_{i,lpd} = x_{i+1,0} \f$.
    !! 
    !! @bug How to ensure that boundaries of neighboring elements are matching 
    !!      seamlessly i.e.  coords(i,lpd) == coords(i+1,0)? In the case of
    !!      gaps points may slip through. In particular this crucial amongst
    !!      MPI partitions.
    !!
    !! @warning In Ses3d-NT r coordinates are upside down. Thus, min and max 
    !!          have to be interchanged 
    !!
    !! @param elms Number of elements 
    !! @param lpd Lagrange polynomial degree 
    !! @param min Lower coordinate bound 
    !! @param max Upper coordinate bound 
    !! @return 1d coordinate array of rank 2 and shape \f$ (elms,lpd+1) \f$ 
    FUNCTION coordinates_1d(elms, lpd, min, max) RESULT(coords)
        USE gll_mod, ONLY : get_knots
        ! Passed dummy arguments 
        INTEGER, INTENT(IN) :: elms, lpd
        REAL(real_kind), INTENT(IN) :: min, max
        ! Result
        REAL(real_kind) :: coords(elms, lpd+1)
        ! Local variables
        REAL(real_kind) :: knots(lpd+1), delta
        INTEGER :: i

        ! Retrieve GLL Points 
        knots(:) = get_knots(lpd)
        ! Scale GLL Points to be in the interval [0,1]
        knots(:) = 0.5_real_kind*(1.0_real_kind + knots(:))

        ! Calculate element width
        delta = (max - min)/elms

        ! Make coordinate vector 
        DO i=0, elms - 1
            coords(i+1,:) = min + (i + knots(:))*delta
        END DO

        ! XXX hack: Make sure the upper bound is equal to max
        coords(elms,lpd+1) = max

    END FUNCTION 


    !> Returns a 3d field filled with theta-coordinates
    !! 
    !! Actually this is not doing much more than a reshape. It translates the 
    !! array of rank 2 representing the 1d coordinate vector - retrieved 
    !! through coordinates_1d() - into an array of rank 6 which corresponds to 
    !! to the 3d coordinate field. 
    !! \f[ 
    !!   \theta^{3d}_{ikmjln} = \theta^{1d}_{ij} \;\;\; \text{for} \; 
    !!      i \in \{ 1 \dots e_\theta \}, \;
    !!      k \in \{ 1 \dots e_\phi \}, \;
    !!      m \in \{ 1 \dots e_r \}, \;
    !!      j,l,n \in \{0\dots lpd \} 
    !! \f] 
    !!
    !! @param elms_theta Number of elements in \f$ \theta \f$ direction 
    !! @param elms_phi Number of elements in \f$ \phi \f$ direction 
    !! @param elms_r Number of elements in \f$ r \f$ direction 
    !! @param lpd Lagrange polynomial degree 
    !! @param min Lower bound of \f$ \theta \f$ coordinates
    !! @param max Upper bound of \f$ \theta \f$ coordinates
    !! @return 3d theta coordinates array of rank 6 and shape 
    !!         \f$ (e_\theta,e_\phi,e_r,lpd+1,lpd+1,lpd+1) \f$ 
    FUNCTION theta_coordinates_3d(elms_theta, elms_phi, elms_r, lpd, &
                                  min, max) RESULT(coords)
        ! Passed dummy arguments 
        INTEGER, INTENT(IN) :: elms_theta, elms_phi, elms_r, lpd
        REAL(real_kind), INTENT(IN) :: min, max
        ! Result
        REAL(real_kind) :: coords(elms_theta, elms_phi, elms_r, &
                                  lpd+1, lpd+1, lpd+1)
        ! Local variables
        REAL(real_kind) :: coords_1d(elms_theta, lpd+1)
        INTEGER :: i, j

        ! Get 1d theta coordinates 
        coords_1d = coordinates_1d(elms=elms_theta, lpd=lpd, min=min, max=max)

        ! Translate theta coordinates to 3d field 
        DO i=1, elms_theta
            DO j=1, lpd+1
                coords(i,:,:,j,:,:) = coords_1d(i,j)
            END DO
        END DO

    END FUNCTION 


    !> Returns a 3d field filled with phi coordinates
    !! 
    !! Actually this is not doing much more than a reshape. It translates the 
    !! array of rank 2 representing the 1d coordinate vector - retrieved 
    !! through coordinates_1d() - into an array of rank 6 which corresponds to 
    !! to the 3d coordinate field. 
    !! \f[ 
    !!   \phi^{3d}_{ikmjln} = \phi^{1d}_{ij} \;\;\; \text{for} \; 
    !!      i \in \{ 1 \dots e_\theta \}, \;
    !!      k \in \{ 1 \dots e_\phi \}, \;
    !!      m \in \{ 1 \dots e_r \}, \;
    !!      j,l,n \in \{0\dots lpd \} 
    !! \f] 
    !!
    !! @param elms_theta Number of elements in \f$ \theta \f$ direction 
    !! @param elms_phi Number of elements in \f$ \phi \f$ direction 
    !! @param elms_r Number of elements in \f$ r \f$ direction 
    !! @param lpd Lagrange polynomial degree 
    !! @param min Lower bound of \f$ \phi \f$ coordinates
    !! @param max Upper bound of \f$ \phi \f$ coordinates
    !! @return 3d theta coordinates array of rank 6 and shape 
    !!         \f$ (e_\theta,e_\phi,e_r,lpd+1,lpd+1,lpd+1) \f$ 
    FUNCTION phi_coordinates_3d(elms_theta, elms_phi, elms_r, lpd, &
                                min, max) RESULT(coords)
        ! Passed dummy arguments 
        INTEGER, INTENT(IN) :: elms_theta, elms_phi, elms_r, lpd
        REAL(real_kind), INTENT(IN) :: min, max
        ! Result
        REAL(real_kind) :: coords(elms_theta, elms_phi, elms_r, &
                                  lpd+1, lpd+1, lpd+1)
        ! Local variables
        REAL(real_kind) :: coords_1d(elms_phi, lpd+1)
        INTEGER :: i, j

        ! Get 1d phi coordinates 
        coords_1d = coordinates_1d(elms=elms_phi, lpd=lpd, min=min, max=max)

        ! Translate phi coordinates to 3d field 
        DO i=1, elms_phi
            DO j=1, lpd+1
                coords(:,i,:,:,j,:) = coords_1d(i,j)
            END DO
        END DO

    END FUNCTION 


    !> Returns a 3d field filled with r-coordinates
    !! 
    !! Actually this is not doing much more than a reshape. It translates the 
    !! array of rank 2 representing the 1d coordinate vector - retrieved 
    !! through coordinates_1d() - into an array of rank 6 which corresponds to 
    !! to the 3d coordinate field. 
    !! \f[ 
    !!   r^{3d}_{ikmjln} = r^{1d}_{ij} \;\;\; \text{with} \; 
    !!      i \in \{ 1 \dots e_\theta \}, \;
    !!      k \in \{ 1 \dots e_\phi \}, \;
    !!      m \in \{ 1 \dots e_r \}, \;
    !!      j,l,n \in \{0\dots lpd \} 
    !! \f] 
    !!
    !! @warning In Ses3d-NT r coordinates are upside down. Thus, min and max 
    !!          have to be interchanged 
    !!
    !! @param elms_theta Number of elements in \f$ \theta \f$ direction 
    !! @param elms_phi Number of elements in \f$ \phi \f$ direction 
    !! @param elms_r Number of elements in \f$ r \f$ direction 
    !! @param lpd Lagrange polynomial degree 
    !! @param min Lower bound of \f$ r \f$ coordinates
    !! @param max Upper bound of \f$ r \f$ coordinates
    !! @return 3d theta coordinates array of rank 6 and shape 
    !!         \f$ (e_\theta,e_\phi,e_r,lpd+1,lpd+1,lpd+1) \f$ 
    FUNCTION r_coordinates_3d(elms_theta, elms_phi, elms_r, lpd, &
                              min, max) RESULT(coords)
        ! Passed dummy arguments 
        INTEGER, INTENT(IN) :: elms_theta, elms_phi, elms_r, lpd
        REAL(real_kind), INTENT(IN) :: min, max
        ! Result
        REAL(real_kind) :: coords(elms_theta, elms_phi, elms_r, &
                                  lpd+1, lpd+1, lpd+1)
        ! Local variables
        REAL(real_kind) :: coords_1d(elms_r, lpd+1)
        INTEGER :: i, j

        ! Get 1d r coordinates 
        coords_1d = coordinates_1d(elms=elms_r, lpd=lpd, min=min, max=max)

        ! Translate r coordinates to 3d field 
        DO i=1, elms_r
            DO j=1, lpd+1
                coords(:,:,i,:,:,j) = coords_1d(i,j)
            END DO
        END DO

    END FUNCTION 

    !> Packs Ses3d-NT's rank 2 arrays into an array of rank one
    !! with dropping duplicates
    PURE FUNCTION pack_r2tor1(r2) RESULT(r1)
        REAL(real_kind), INTENT(IN) :: r2(:,:)
        REAL(real_kind) :: r1(SIZE(r2, DIM=1)*(SIZE(r2, DIM=2)-1)+1)
        INTEGER :: nx, lpd

        lpd = SIZE( r2, DIM=2 )
        nx = SIZE( r2, DIM=1 )

        r1 = [ TRANSPOSE( r2(:,:lpd-1) ), r2(nx,lpd) ]

    END FUNCTION 


END MODULE coordinates_mod
