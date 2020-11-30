!> @file
!! Ses3d-NT - simulation of elastic wave propagation in spherical sections
!!
!! (c) by Stefan Mauerberger
!!    and Lion Krischer
!!    and Maksym Melnyk
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
!! $Date: 2013-11-02 15:21:34 +0100 (Sat, 02 Nov 2013) $
!! $Author: mauerberger $
!! $Revision: 757 $
!! @copyright GNU General Public License version 3 or later



!> This module provides few low-level procedures to convert coordinates among
!! spherical and Cartesian coordinate systems.
!!
!! The module expects to find the parameters parameters_mod::real_kind,
!! parameters_mod::earth_radius and parameters_mod::pi in a module named
!! parameters_mod.
MODULE coordinate_utilities_mod
    USE parameters_mod, ONLY : real_kind, earth_radius, pi

    IMPLICIT NONE

    PRIVATE

    PUBLIC :: deg2rad, rad2deg, &
              lat2colat, colat2lat, &
              xyz2tpr, tpr2xyz, &
              depth2radius, radius2depth

CONTAINS

    !> Converting latitude to colatitude
    !!
    !! @param lat latitude [deg]
    !! @return colatitude [deg]
    REAL(real_kind) ELEMENTAL FUNCTION lat2colat(lat) RESULT(colat)
        REAL(real_kind), INTENT(IN) :: lat
        colat = 90.0 - lat
    END FUNCTION

    !> Converting colatitude to latitude
    !!
    !! @param colat colatitude [deg]
    !! @return latitude [deg]
    REAL(real_kind) ELEMENTAL FUNCTION colat2lat(colat) RESULT(lat)
        REAL(real_kind), INTENT(IN) :: colat
        lat = 90.0 - colat
    END FUNCTION

    !> Converting angles from radians to degrees
    !!
    !! @param rad Angle in radians
    !! @return Angle in degrees
    REAL(real_kind) ELEMENTAL FUNCTION rad2deg(rad) RESULT(deg)
        REAL(real_kind), INTENT(IN) :: rad
        deg = rad * 180.0 / pi
    END FUNCTION

    !> Converting angles from degrees to radians
    !!
    !! @param deg Angle in degrees
    !! @return Angle in radians
    REAL(real_kind) ELEMENTAL FUNCTION deg2rad(deg) RESULT(rad)
        REAL(real_kind), INTENT(IN) :: deg
        rad = deg * pi / 180.0
    END FUNCTION

    !> Converts a depth value (measured form the Earth's surface)
    !! into its corresponding radius
    !!
    !! @param depth (positive values; starting from the Earth's surface) [meters]
    !! @return radius counting from the Earth's center [meters]
    REAL(real_kind) ELEMENTAL FUNCTION depth2radius( depth ) RESULT( radius )
        REAL(real_kind), INTENT(IN) :: depth
        radius = earth_radius - depth
    END FUNCTION

    !> Converts a radius into its depth value with respect to the Earth's radius
    !!
    !> @param radius radius [meters]
    !> @return depth from the Earth's surface; positive values [meters]
    REAL(real_kind) ELEMENTAL FUNCTION radius2depth( radius ) RESULT( depth )
        REAL(real_kind), INTENT(IN) :: radius
        depth = earth_radius - radius
    END FUNCTION

    !> Converts Cartesian coordinates (x,y,z) into spherical
    !! coordinates (theta,phi,r)
    !!
    !! In terms of the standard arctan function, whose range
    !! is (-pi/2, pi/2), it can be expressed as follows:
    !!
    !! \f[ \operatorname{atan2}(y, x) = \begin{cases}
    !! \arctan\left(\frac y x\right) & \qquad x > 0 \\
    !! \arctan\left(\frac y x\right) + \pi& \qquad y \ge 0 , x < 0 \\
    !! \arctan\left(\frac y x\right) - \pi& \qquad y < 0 , x < 0 \\
    !! +\frac{\pi}{2} & \qquad y > 0 , x = 0 \\
    !! -\frac{\pi}{2} & \qquad y < 0 , x = 0 \\
    !! \text{undefined} & \qquad y = 0, x = 0
    !! \end{cases} \f]
    !!
    !! @param x x-component in Cartesian coordinates
    !! @param y y-component in Cartesian coordinates
    !! @param z z-component in Cartesian coordinates
    !! @param theta inclination in spherical coordinates [rad]
    !! @param phi azimuth in spherical coordinates [rad]
    !! @param r radial component in spherical coordinates
    !! @return spherical coordinates vector: [theta, phi, r]
    ELEMENTAL SUBROUTINE xyz2tpr( x, y, z, theta, phi, r )
        REAL(real_kind), INTENT(IN) :: x, y, z
        REAL(real_kind), INTENT(OUT) :: theta, phi, r
        r     = SQRT( x**2 + y**2 + z**2 )
        theta = ACOS( z / r )
        phi   = ATAN2( y, x )
    END SUBROUTINE

    !> Converts spherical coordinates (theta,phi,r) into
    !! Cartesian coordinates (x,y,z)
    !!
    !! \f[
    !!  \begin{pmatrix} x\\ y\\ z \end{pmatrix} =
    !!  \begin{pmatrix} r \, \sin \theta \, \cos \phi \\
    !!                  r \, \sin \theta \, \sin \phi \\
    !!                  r \, \cos \theta \end{pmatrix}
    !! \f]
    !!
    !! @param theta inclination in spherical coordinates [rad]
    !! @param phi azimuth in spherical coordinates [rad]
    !! @param r radial-component in spherical coordinates
    !! @param x x-component in Cartesian coordinates
    !! @param y y-component in Cartesian coordinates
    !! @param z z-component in Cartesian coordinates
    !! @return Cartesian coordinates vector: [x, y, z]
    ELEMENTAL SUBROUTINE tpr2xyz( theta, phi, r, x, y, z )
        REAL(real_kind), INTENT(IN) :: theta, phi, r
        REAL(real_kind), INTENT(OUT) :: x, y, z
        x = r * sin( theta ) * cos( phi )
        y = r * sin( theta ) * sin( phi )
        z = r * cos( theta )
    END SUBROUTINE

END MODULE coordinate_utilities_mod
