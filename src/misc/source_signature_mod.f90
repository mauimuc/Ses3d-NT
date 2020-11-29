!> @file
!! Ses3d-NT - simulation of elastic wave propagation in spherical sections
!!
!! (c) by Stefan Mauerberger <mauerberger@geophysik.uni-muenchen.de>
!!    and Maksym Melnyk <mmelnyk@geophysik.uni-muenchen.de>
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
!! $Date: 2013-11-01 18:19:11 +0100 (Fri, 01 Nov 2013) $
!! $Author: mmelnyk $
!! $Revision: 754 $
!> @copyright GNU General Public License version 3 or later



!> This module is dedicated to generate source wavelets. It provides
!! functions for the Ricker wavelet, sin**3 and  t*exp(t**2), a peak and
!! Gaussian pulse.
MODULE source_signature_mod
    USE parameters_mod, ONLY : real_kind, pi

    IMPLICIT NONE

    PRIVATE

    PUBLIC :: ricker, sin_3, delta, brune, gauss

CONTAINS


    !> @brief A Ricker wavelet
    !> @details Returns the value of a Ricker wavelet at time \f$t\f$ as
    !! \f[ s(t) = (1-2\left(\frac{\pi(t-T)}{\sigma}\right)^2) \, \exp( -\left(\frac{\pi(t-T)}{\sigma}\right)^2)\f]
    !! with
    !> @param time time \f$t\f$ in seconds,
    !> @param onset time of source onset \f$T\f$ in seconds
    !> @param width and width \f$\sigma\f$.
    !> @result s value of the Ricker-wavelet at time \f$t\f$
    ELEMENTAL FUNCTION ricker(time, onset, width) RESULT(s)
        REAL(real_kind), INTENT(IN) :: time, onset, width
        REAL(real_kind) :: a2, s

        a2 = -(pi/width)**2

        s = 2.0 * ( 0.5 + a2 * (time-onset)**2 ) * EXP( a2 * (time-onset)**2 )

    END FUNCTION ricker


    !> @brief sin(.)-cubed wavelet
    !> @details Returns the value of a sin()-cubed wavelet at time \f$t\f$
    !! \f[ s(t) = \begin{cases}
    !!      0& \text{for}~t \leq T\\
    !!      \sin^3\left(\frac{\pi (t-T)}{\sigma}\right)& \text{for}~T < t < T+\sigma\\
    !!      0& \text{for}~t \geq T+\sigma
    !!      \end{cases}  \f]
    !! with
    !> @param time time \f$t\f$ in seconds,
    !> @param onset time of source onset \f$T\f$ in seconds
    !> @param width and duration of source signal \f$\sigma\f$ seconds
    !> @result s Value of the sin()-cubed wavelet at time \f$t\f$
    ELEMENTAL FUNCTION sin_3(time, onset, width) RESULT(s)
        REAL(real_kind), INTENT(IN) :: time, onset, width
        REAL(real_kind) :: s

        IF ( onset < time .AND. &
             time < onset+width ) THEN
            s = sin(pi*(time-onset)/width)**3
        ELSE
            s = 0.0
        END IF

    END FUNCTION sin_3


    !> @brief Peak at time T
    !> @details Returns 1.0 if time == onset else 0.0
    !! \f[ s(t) = \delta(t-T)\f] with
    !> @param time time \f$t\f$ in seconds and
    !> @param onset time of source onset \f$T\f$ in seconds
    !> @result s Returns 1.0 if time == onset else 0.0
    !> @todo Make it more forgiving with respect to floating-point
    !! representations. It is somewhat fragile because time has to be
    !! precisely the same than onset
    ELEMENTAL FUNCTION delta(time, onset) RESULT(s)
        REAL(real_kind), INTENT(IN) :: time, onset
        REAL(real_kind) :: s

        IF ( time == onset ) THEN
            s = 1.0
        ELSE
            s = 0.0
        END IF

    END FUNCTION delta



    !> @brief Brune Wavelet.
    !> @details Source-time-function used in SISMOWINE as
    !! \f[ s(t) = \begin{cases}
    !!      0& \text{for}~t < T\\
    !!      \frac{t-T}{\sigma^2} \exp(-\frac{t-T}{\sigma})& \text{for}~t \geq T
    !!      \end{cases}  \f]
    !! with
    !> @param time Time \f$t\f$ in seconds,
    !> @param onset time of source onset \f$T\f$ in seconds
    !> @param width width \f$\sigma\f$ in seconds
    !> @result s Returns 1.0 if time == onset else 0.0
    ELEMENTAL FUNCTION brune(time, onset, width) RESULT(s)
        REAL(real_kind), INTENT(IN) :: time, onset, width
        REAL(real_kind) :: s

        IF ( time .ge. onset ) THEN
            s = (time-onset)/width**2*exp(-(time-onset)/width)
        ELSE
            s=0.0
        END IF

    END FUNCTION brune


    !> @brief Gaussian pulse.
    !> @details Returns the value of a Gaussian pulse at time \f$t\f$ as
    !! \f[ s(t) =  \exp(-\frac{(t-T)^2}{\sigma^2})\f]
    !! with
    !> @param time time \f$t\f$ in seconds,
    !> @param onset onset time \f$T\f$ in seconds
    !> @param width and width \f$\sigma\f$ in seconds
    ELEMENTAL FUNCTION gauss(time, onset, width) RESULT(s)
        REAL(real_kind), INTENT(IN) :: time, width, onset
        REAL(real_kind) :: s

        s = exp(-((time-onset)/width)**2)

    END FUNCTION gauss


END MODULE source_signature_mod
