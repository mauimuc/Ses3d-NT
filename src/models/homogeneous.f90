! SES3D - simulation of elastic wave propagation in spherical sections

! (c) by Stefan Mauerberger <mauerberger@geophysik.uni-muenchen.de>
!        Stefan Wenk <wenk@geophysik.uni-muenchen.de>

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

! Last changed: $Date: 2013-07-29 15:37:16 +0200 (Mon, 29 Jul 2013) $
! By: $Author: mauerberger $
! Revision: $Revision: 543 $


!------------------------------------------------------------------------------
! homogeneous model
!------------------------------------------------------------------------------
MODULE homogeneous_mod
    USE parameters_mod, ONLY : real_kind

    IMPLICIT NONE

    PRIVATE

    PUBLIC :: homogeneous

CONTAINS

    ! homogeneous model
    ELEMENTAL SUBROUTINE homogeneous( rho, mu, lambda, a, b, c, &
                                      qkappa, qmu )
        REAL(real_kind), INTENT(OUT), OPTIONAL :: rho, mu, lambda, a, b, c, &
                                           qkappa, qmu

        REAL(real_kind), PARAMETER :: rho_    = 3000.0
        REAL(real_kind), PARAMETER :: mu_     = 7.0e10
        REAL(real_kind), PARAMETER :: lambda_ = 1.1e11

        REAL(real_kind), PARAMETER :: a_ = 0.0
        REAL(real_kind), PARAMETER :: b_ = 0.0
        REAL(real_kind), PARAMETER :: c_ = 0.0

        REAL(real_kind), PARAMETER :: qkappa_ = 0.0
        REAL(real_kind), PARAMETER :: qmu_    = 0.0


        IF ( PRESENT( rho ) ) &
            rho    = rho_
        IF ( PRESENT( mu ) ) &
            mu     = mu_
        IF ( PRESENT( lambda ) ) &
            lambda = lambda_
        IF ( PRESENT( qkappa ) ) &
            qkappa = qkappa_
        IF ( PRESENT( qmu ) ) &
            qmu = qmu_
        IF ( PRESENT( a ) ) &
            a      = a_
        IF ( PRESENT( b ) ) &
            b      = b_
        IF ( PRESENT( c ) ) &
            c      = c_

    END SUBROUTINE homogeneous

END MODULE homogeneous_mod


