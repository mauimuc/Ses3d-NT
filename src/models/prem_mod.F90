!> @file
!! Ses3d-NT - simulation of elastic wave propagation in spherical sections
!!
!! (c) by Stefan Mauerberger
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
!! $Date: 2013-11-24 14:42:01 +0100 (Sun, 24 Nov 2013) $
!! $Author: mmelnyk $
!! $Revision: 806 $
!> @copyright GNU General Public License version 3 or later



!> Module Preliminary Reference Earth Model (PREM) by Dziewonski & Anderson, 1981
!!
!! Is an average Earth model that incorporates anelastic dispersion and
!! anisotropy and therefore it is frequency-dependent and transversely
!! isotropic for the upper mantle.
!!
!! For further details please see the model websites:
!!  * http://www.iris.edu/dms/products/emc-prem/
MODULE prem_mod
#ifndef STAND_ALONE
    USE parameters_mod, ONLY : real_kind
#endif

    IMPLICIT NONE

    PRIVATE

#ifdef STAND_ALONE
    INTEGER, PARAMETER :: real_kind = 8 ! real_kind
#endif

    ! Defines stratified earth structure
    REAL(real_kind), PARAMETER :: R_EARTH          = 6371000.0, &
                                  ROCEAN           = 6368000.0, &
                                  RMIDDLE_CRUST    = 6356000.0, &
                                  RMOHO            = 6346600.0, &
                                  R80              = 6291000.0, &
                                  R220             = 6151000.0, &
! Andreas uses a different value: R220             = 6140200.0, &
                                  R400             = 5971000.0, &
                                  R600             = 5771000.0, &
                                  R670             = 5701000.0, &
                                  R771             = 5600000.0, &
                                  RTOPDDOUBLEPRIME = 3630000.0, &
                                  RCMB             = 3480000.0, &
!                                 RICB             = 1221000.0
                                  RICB             = 1221500.0

    PUBLIC :: prem, prem_iso

CONTAINS

!> Subroutine Isotropic Preliminary Reference Earth Model
!!
!> @param[in] r Radius corresponding to the depth [km]
!> @param[out] rho Density [Mg/km3]
!> @param[out] drhodr ....
!> @param[out] vs Velocity of the S-waves [km/s]
!> @param[out] vp Velocity of the P-waves [km/s]
!> @param[out] qkappa Bulk quality factor
!> @param[out] qmu Shear quality factor
!> @param[out] nocrust logical if TRUE crustal structure becomes simplified ( default FALSE )
!> @param[out] noocean logical if TRUE water-layer is replaced by crustal material ( default TRUE )
ELEMENTAL SUBROUTINE prem_iso( r, rho, drhodr, vs, vp, qkappa, qmu, &
                               nocrust, noocean )
    REAL(real_kind), INTENT(IN) :: r

    LOGICAL, INTENT(IN), OPTIONAL :: nocrust, noocean
    REAL(real_kind), INTENT(OUT), OPTIONAL :: rho, drhodr, vs, vp, qkappa, qmu

    ! Local variables
    LOGICAL :: nocrust_, noocean_
    REAL(real_kind) :: x, rho_, drhodr_, vs_, vp_, qkappa_, qmu_


    ! TODO move to subroutine prem_no_crust
    nocrust_ = .FALSE.
    IF ( PRESENT(nocrust) ) &
        nocrust_ = nocrust

    ! TODO move to subroutine prem_no_oecan
    noocean_ = .TRUE.
    IF ( PRESENT(noocean) ) &
        noocean_ = noocean

    ! Normalize radius
    x = ABS( r / R_EARTH )

    ! --- inner core
    IF ( r >= 0.0 .AND. r <= RICB ) THEN
        drhodr_ =-2.00*8.83810*x
        rho_    = 13.08850 - 8.83810*x**2
        vp_     = 11.26220 - 6.36400*x**2
        vs_     = 3.66780 - 4.44750*x**2
        qmu_    = 84.60
        qkappa_ = 1327.70

    ! --- outer core
    ELSE IF ( r > RICB .and. r <= RCMB ) THEN
        drhodr_ =-1.263800 - 2.0*3.64260*x - 3.00*5.52810*x**2
        rho_    = 12.58150 - 1.26380*x - 3.64260*x**2 -  5.52810*x**3
        vp_     = 11.04870 - 4.03620*x + 4.80230*x**2 - 13.57320*x**3
        vs_     = 0.0 ! TINY(0.0_real_kind)
        qmu_    = HUGE(0.0_real_kind) ! TODO should be +Inv
        qkappa_ = 57823.00

    ! --- D" at the base of the mantle
    ELSE IF ( r > RCMB .AND. r <= RTOPDDOUBLEPRIME ) THEN
        drhodr_ =-6.47610+2.00*5.52830*x-3.00*3.08070*x*x
        rho_    = 7.95650-6.47610*x+5.52830*x*x-3.08070*x*x*x
        vp_     = 15.38910-5.31810*x+5.52420*x*x-2.55140*x*x*x
        vs_     = 6.92540+1.46720*x-2.08340*x*x+0.97830*x*x*x
        qmu_    = 312.00
        qkappa_ = 57823.00

    ! --- mantle: from top of D" to d670
    ELSE IF ( r > RTOPDDOUBLEPRIME .AND. r <= R771 ) THEN
        drhodr_ =-6.47610+2.00*5.52830*x-3.00*3.08070*x*x
        rho_    = 7.95650-6.47610*x+5.52830*x*x-3.08070*x*x*x
        vp_     = 24.95200-40.46730*x+51.48320*x*x-26.64190*x*x*x
        vs_     = 11.16710-13.78180*x+17.45750*x*x-9.27770*x*x*x
        qmu_    = 312.00
        qkappa_ = 57823.00

    !
    ELSE IF ( r > R771 .AND. r <= R670 ) THEN
        drhodr_ =-6.47610+2.00*5.52830*x-3.00*3.08070*x*x
        rho_    = 7.95650-6.47610*x+5.52830*x*x-3.08070*x*x*x
        vp_     = 29.27660-23.60270*x+5.52420*x*x-2.55140*x*x*x
        vs_     = 22.34590-17.24730*x-2.08340*x*x+0.97830*x*x*x
        qmu_    = 312.00
        qkappa_ = 57823.00

    ! --- mantle: above d670
    ELSE IF ( r > R670 .AND. r <= R600 ) THEN
        drhodr_ =-1.48360
        rho_    = 5.31970-1.48360*x
        vp_     = 19.09570-9.86720*x
        vs_     = 9.98390-4.93240*x
        qmu_    = 143.00
        qkappa_ = 57823.00

    ! ---
    ELSE IF ( r > R600 .AND. r <= R400 ) THEN
        drhodr_ =-8.02980
        rho_    = 11.24940-8.02980*x
        vp_     = 39.70270-32.61660*x
        vs_     = 22.35120-18.58560*x
        qmu_    = 143.00
        qkappa_ = 57823.00

    ! ---
    ELSE IF ( r > R400 .AND. r <= R220 ) THEN
        drhodr_ =-3.80450
        rho_    = 7.10890-3.80450*x
        vp_     = 20.39260-12.25690*x
        vs_     = 8.94960-4.45970*x
        qmu_    = 143.00
        qkappa_ = 57823.00

    ! --- LVZ
    ELSE IF ( R > R220 .AND. R <= R80 ) THEN
        drhodr_ = 0.69240
        rho_    = 2.69100 + 0.69240*x
        vp_     = 4.18750 + 3.93820*x
        vs_     = 2.15190 + 2.34810*x
        qmu_    = 80.00
        qkappa_ = 57823.00

    ! --- LID
    ELSE IF ( r > R80 .AND. r <= RMOHO ) THEN
        drhodr_ = 0.69240
        rho_    = 2.69100 + 0.69240*x
        vp_     = 4.18750 + 3.93820*x
        vs_     = 2.15190 + 2.34810*x
        qmu_    = 600.00
        qkappa_ = 57823.00

    ! --- Mohorovicic discontinuity
    ELSE IF ( r > RMOHO .AND. r <= RMIDDLE_CRUST ) THEN
        drhodr_ = 0.00
        rho_    = 2.90
        vp_     = 6.80
        vs_     = 3.90
        qmu_    = 600.00
        qkappa_ = 57823.00

    ! --- Middle crust
    ELSE IF ( r > RMIDDLE_CRUST .AND. r <= ROCEAN ) THEN
        drhodr_ = 0.00
        rho_    = 2.60
        vp_     = 5.80
        vs_     = 3.20
        qmu_    = 600.00
        qkappa_ = 57823.00

    ! --- Crust
    ELSE IF ( r > ROCEAN .AND. r <= R_EARTH ) THEN
        drhodr_ = 0.00
        rho_    = 1.020
        vp_     = 1.450
        vs_     = 0.00 ! TINY(0.0_real_kind)
        qmu_    = HUGE(0.0_real_kind) ! Should be +Inv
        qkappa_ = 57823.00

    ELSE
        drhodr_ = 0.00
        rho_    = 0.00
        vp_     = 0.00
        vs_     = 0.00 ! TINY(0.0_real_kind)
        qmu_    = HUGE(0.0_real_kind) ! Should be +Inv
        qkappa_ = 0.00

    ENDIF

    ! TODO move to subroutine prem_no_oecan
    if ( noocean_ .AND. ( r > ROCEAN ) ) then
        drhodr_=0.00
        rho_=2.60
        vp_=5.80
        vs_=3.20
        qmu_=600.00
        qkappa_=57823.00
    end if

    ! TODO move to subroutine prem_no_crust
    if ( nocrust_ .AND. ( r > R80 ) ) then
        drhodr_=0.69240
        rho_=2.69100+0.69240*(RMOHO / R_EARTH)
        vp_=4.18750+3.93820*(RMOHO / R_EARTH)
        vs_=2.15190+2.34810*(RMOHO / R_EARTH)
        qmu_=600.00
        qkappa_=57823.00
    end if

    IF ( PRESENT( rho ) ) &
        rho    = rho_
    IF ( PRESENT( drhodr ) ) &
        drhodr = drhodr_
    IF ( PRESENT( vp ) ) &
        vp     = vp_
    IF ( PRESENT( vs ) ) &
        vs     = vs_
    IF ( PRESENT( qkappa ) ) &
        qkappa = qkappa_
    IF ( PRESENT( qmu ) ) &
        qmu = qmu_

  END SUBROUTINE prem_iso



    !> Subroutine Preliminary Reference Earth Model (PREM)
    !!
    !> @param[in] r Radius corresponding to the depth [km]
    !> @param[out] rho Density [Mg/km3]
    !> @param[out] drhodr ....
    !> @param[out] vpv Vertical velocity of the P-waves [km/s]
    !> @param[out] vph Horizontal velocity of the P-waves [km/s]
    !> @param[out] vsv Vertical velocity of the S-waves [km/s]
    !> @param[out] vsh Horizontal velocity of the S-waves [km/s]
    !> @param[out] vp Velocity of the P-waves [km/s]
    !> @param[out] eta [dimensionless]
    !> @param[out] qkappa Bulk quality factor
    !> @param[out] qmu Shear quality factor
    ELEMENTAL SUBROUTINE prem( r, drhodr, rho, vpv, vph, vsv, vsh, &
                               eta, qkappa, qmu )
        REAL(real_kind), INTENT(IN) :: r
        REAL(real_kind), INTENT(OUT), OPTIONAL :: drhodr, rho, vpv, vph, vsv, &
                                                  vsh, eta, qkappa, qmu
        ! Local variables
        REAL(real_kind) :: x, drhodr_, rho_, vsh_, vsv_, vph_, vpv_, &
                    eta_, qkappa_, qmu_
        ! Normalize radius
        x = ABS( r / R_EARTH )

        ! --- LVZ (Low-velocity zone)
        IF ( r > R220 .AND. r <= R80 ) THEN
            drhodr_ = 0.69240
            rho_    = 2.69100 + 0.69240*x
            vph_    = 3.59080 + 4.61720*x
            vpv_    = 0.83710 + 7.21800*x
            vsh_    =-1.08390 + 5.71760*x
            vsv_    = 5.85820 - 1.46780*x
            eta_    = 3.36870 - 2.47780*x
            qmu_    = 80.00
            qkappa_ = 57823.00

        ! --- LID (main part of seismic lithosphere)
        ELSE IF ( r > R80 .AND. r <= RMOHO ) THEN
            drhodr_ = 0.69240
            rho_    = 2.69100 + 0.69240*x
            vph_    = 3.59080 + 4.61720*x
            vpv_    = 0.83710 + 7.21800*x
            vsh_    =-1.08390 + 5.71760*x
            vsv_    = 5.85820 - 1.46780*x
            eta_    = 3.36870 - 2.47780*x
            qmu_    = 600.00
            qkappa_ = 57823.00

        ELSE
            CALL prem_iso( r, rho=rho_, drhodr=drhodr_, vp=vpv_, vs=vsv_, &
                           qmu=qmu_, qkappa=qkappa_ )
            vph_ = vpv_
            vsh_ = vsv_
            eta_ = 1.0

        END IF

        IF ( PRESENT( rho ) ) &
            rho    = rho_
        IF ( PRESENT( drhodr ) ) &
            drhodr = drhodr_
        IF ( PRESENT( vpv ) ) &
            vpv    = vpv_
        IF ( PRESENT( vph ) ) &
            vph    = vph_
        IF ( PRESENT( vsv ) ) &
            vsv    = vsv_
        IF ( PRESENT( vsh ) ) &
            vsh    = vsh_
        IF ( PRESENT( eta ) ) &
            eta    = eta_
        IF ( PRESENT( qmu ) ) &
            qmu    = qmu_
        IF ( PRESENT( qkappa ) ) &
            qkappa = qkappa_

    END SUBROUTINE prem


END MODULE prem_mod

