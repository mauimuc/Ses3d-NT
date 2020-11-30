!> @file
!! Ses3d-NT - simulation of elastic wave propagation in spherical sections
!!
!! (c) by Stefan Mauerberger
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
!! $Date: 2013-10-30 22:07:04 +0100 (Wed, 30 Oct 2013) $
!! $Author: mauerberger $
!! $Revision: 746 $
!! @copyright GNU General Public License version 3 or later



!> This module provides the GLL collocation points and corresponding weights
!! for GLL-integration of degree 1 to 7.
!!
!! The \f$N+1\f$ GLL nodes are given by the \f$N-1\f$ roots of the Legendre
!! polynomials \f$P_N(x)\f$ of degree \f$N\f$. The boundaries of the interval
!! \f$[-1:1]\f$ are additionally included as node points
!! \f[
!!     (1-x_i^2) \, \partial_x \text P_N(x_i) = 0 \;\;\; i \in 0 \dots N
!! \f]
!!
!! For the GLL Quadrature the weights are defined by the integral over the
!! Lagrange polynomials. Including the definition of the \f$N+1\f$ GLL
!! collocation points \f$x_i\f$ the corresponding weights \f$w_i\f$ can be
!! written in a closed form
!! \f[
!!    w_i = \int_{-1}^1 \ell^N_i(x) \; \mathrm dx =
!!        \frac{2}{N (N + 1)} {P_N(x_i)}^{-2} \;\;\; i \in 0 \dots N
!! \f]
!! where \f$\ell^N_i(x)\f$ denotes the \f$i\f$-th Lagrange polynomial of degree
!! \f$N\f$. For Lagrange interpolation see interp_mod.
!!
!! @todo Tests gll_test.f90
!!
!! The archived accuracy is at machine precision. The module expects to find
!! the parameter parameters_mod::real_kind in a module named parameters_mod
!! which defines the accuracy of floating point numbers.
MODULE gll_mod
    USE parameters_mod, ONLY : rk => real_kind

    IMPLICIT NONE

    PRIVATE

    !> GLL collocation points of degree 2
    !! @todo check if those dummy values work
    REAL(rk), PARAMETER, DIMENSION(*) :: knots_1 = &
        [-1.0_rk, &
          1.0_rk ]
    !> Weights for GLL integration of degree 2
    !! @todo check if those dummy values work
    REAL(rk), PARAMETER, DIMENSION(*) :: weights_1 = &
        [ 0.5_rk, &
          0.5_rk ]

    !> GLL collocation points of degree 2
    REAL(rk), PARAMETER, DIMENSION(*) :: knots_2 = &
        [-1.0_rk, &
          0.0_rk, &
          1.0_rk ]
    !> Weights for GLL integration of degree 2
    REAL(rk), PARAMETER, DIMENSION(*) :: weights_2 = &
        [ 1.0_rk / 3.0_rk, &
          4.0_rk / 3.0_rk, &
          1.0_rk / 3.0_rk ]

    !> GLL collocation points of degree 3
    REAL(rk), PARAMETER, DIMENSION(*) :: knots_3 = &
        [-1.0_rk, &
         -SQRT( 1.0_rk / 5.0_rk ), &
          SQRT( 1.0_rk / 5.0_rk ), &
          1.0_rk ]
    !> Weights for GLL integration of degree 3
    REAL(rk), PARAMETER, DIMENSION(*) :: weights_3 = &
        [ 1.0_rk / 6.0_rk, &
          5.0_rk / 6.0_rk, &
          5.0_rk / 6.0_rk, &
          1.0_rk / 6.0_rk ]

    !> GLL collocation points of degree 4
    REAL(rk), PARAMETER, DIMENSION(*) :: knots_4 = &
        [-1.0_rk, &
         -SQRT( 3.0_rk / 7.0_rk ), &
          0.0_rk, &
          SQRT( 3.0_rk / 7.0_rk ), &
          1.0_rk ]
    !> Weights for GLL integration of degree 4
    REAL(rk), PARAMETER, DIMENSION(*) :: weights_4 = &
        [  1.0_rk / 10.0_rk, &
          49.0_rk / 90.0_rk, &
          32.0_rk / 45.0_rk, &
          49.0_rk / 90.0_rk, &
           1.0_rk / 10.0_rk ]

    !> GLL collocation points of degree 5
    REAL(rk), PARAMETER, DIMENSION(*) :: knots_5 = &
        [-1.0_rk, &
         -SQRT( 2.0_rk * SQRT( 7.0_rk ) + 7.0_rk ) / SQRT( 21.0_rk ), &
         -SQRT( 7.0_rk - 2.0_rk * SQRT( 7.0_rk ) ) / SQRT( 21.0_rk ), &
          SQRT( 7.0_rk - 2.0_rk * SQRT( 7.0_rk ) ) / SQRT( 21.0_rk ), &
          SQRT( 2.0_rk * SQRT( 7.0_rk) + 7.0_rk  ) / SQRT( 21.0_rk ), &
          1.0_rk ]
    !> Weights for GLL integration of degree 5
    REAL(rk), PARAMETER, DIMENSION(*) :: weights_5 = &
        [  1.0_rk / 15.0_rk, &
          63.0_rk / (10.0_rk * SQRT( 7.0_rk ) + 140.0_rk ), &
         -63.0_rk / (10.0_rk * SQRT( 7.0_rk ) - 140.0_rk ), &
         -63.0_rk / (10.0_rk * SQRT( 7.0_rk ) - 140.0_rk ), &
          63.0_rk / (10.0_rk * SQRT( 7.0_rk ) + 140.0_rk ), &
           1.0_rk / 15.0_rk ]

    !> GLL collocation points of degree 6
    REAL(rk), PARAMETER, DIMENSION(*) :: knots_6 = &
        [-1.0_rk, &
         -SQRT( 2.0_rk * SQRT( 15.0_rk ) + 15.0_rk ) / SQRT( 33.0_rk ), &
         -SQRT( 15.0_rk - 2.0_rk * SQRT( 15.0_rk ) ) / SQRT( 33.0_rk ), &
          0.0_rk, &
          SQRT( 15.0_rk - 2.0_rk * SQRT( 15.0_rk ) ) / SQRT( 33.0_rk ), &
          SQRT( 2.0_rk * SQRT( 15.0_rk ) + 15.0_rk ) / SQRT( 33.0_rk ), &
          1.0_rk ]
    !> Weights for GLL integration of degree 6
    REAL(rk), PARAMETER, DIMENSION(*) :: weights_6 = &
        [ 1.0_rk / 21.0_rk, &
          14641.0_rk / ( 2450.0_rk * SQRT( 15.0_rk ) + 43400.0_rk ), &
         -14641.0_rk / ( 2450.0_rk * SQRT( 15.0_rk ) - 43400.0_rk ), &
          256.0_rk / 525.0_rk, &
         -14641.0_rk / ( 2450.0_rk * SQRT( 15.0_rk ) - 43400.0_rk ), &
          14641.0_rk / ( 2450.0_rk * SQRT( 15.0_rk ) + 43400.0_rk ), &
          1.0_rk / 21.0_rk ]

    !> GLL collocation points of degree 7
    !! @todo must be replaced by its analytic expressions
    REAL(rk), PARAMETER, DIMENSION(*) :: knots_7 = &
        [-1.0_rk, &
         -0.8717401485096066_rk, &
         -0.5917001814331423_rk, &
         -0.2092992179024789_rk, &
          0.2092992179024789_rk, &
          0.5917001814331423_rk, &
          0.8717401485096066_rk, &
          1.0_rk ]
    !> Weights for GLL integration of degree 7
    REAL(rk), PARAMETER, DIMENSION(*) :: weights_7 = &
        [ 0.0357142857142857_rk, &
          0.2107042271435061_rk, &
          0.3411226924835044_rk, &
          0.4124587946587038_rk, &
          0.4124587946587038_rk, &
          0.3411226924835044_rk, &
          0.2107042271435061_rk, &
          0.0357142857142857_rk ]

    PUBLIC :: get_knots, get_weights, get_dlgll

CONTAINS

    !> Returns the GLL collocation points from degree 1 to 7
    !!
    !! * If lpd is smaller than 1 an zero-size array is returned
    !! * If lpd is bigger tan 7 an empty array of size lpd+1
    !!
    !! @param lpd Lagrange polynomial degree
    !! @result The GLL collocation points; An array of size lpd+1
    PURE FUNCTION get_knots(lpd) RESULT(knots)
        INTEGER, INTENT(IN) :: lpd
        REAL(rk) :: knots(lpd+1)

        SELECT CASE (lpd)
            CASE (1)
                knots(:) = knots_1
            CASE (2)
                knots(:) = knots_2
            CASE (3)
                knots(:) = knots_3
            CASE (4)
                knots(:) = knots_4
            CASE (5)
                knots(:) = knots_5
            CASE (6)
                knots(:) = knots_6
            CASE (7)
                knots(:) = knots_7
            CASE DEFAULT
                knots(:) = 0.0
            CASE (:0)
                ! Array of size 0
        END SELECT

    END FUNCTION get_knots


    !> Returns weights used for the so called GLL integration for degree 1 to 7
    !!
    !! * If lpd is smaller than 1 an zero-size array is returned
    !! * If lpd is bigger tan 7 an empty array of size lpd+1
    !!
    !! @param lpd Lagrange polynomial degree
    !! @result The GLL collocation points; An array of size lpd+1
    PURE FUNCTION get_weights(lpd) RESULT(weights)
        INTEGER, INTENT(IN) :: lpd
        REAL(rk) :: weights(lpd+1)

        SELECT CASE ( lpd )
            CASE (1)
                weights(:) = weights_1
            CASE (2)
                weights(:) = weights_2
            CASE (3)
                weights(:) = weights_3
            CASE (4)
                weights(:) = weights_4
            CASE (5)
                weights(:) = weights_5
            CASE (6)
                weights(:) = weights_6
            CASE (7)
                weights(:) = weights_7
            CASE DEFAULT
                weights(:) = 0.0
            CASE (:0)
                ! Array of size 0
        END SELECT

    END FUNCTION get_weights


    !> Returns the derivatives of GLL polynomials at collocation point j
    !!
    !! @todo Document how the derivative is obtained
    !! @todo This is super cryptic it must become nicer
    !!
    !! @param n Degree of the Lagrange polynomial
    !! @param i Index of the Lagrange polynomial
    !! @param j Index of the collocation point where the derivative is evaluated
    !! @returns 1st derivative of the Lagrange polynomial at the j-th
    !! collocation point
    ELEMENTAL FUNCTION get_dlgll(n,i,j) RESULT(y)
        INTEGER, INTENT(IN) :: n, i, j
        REAL(rk) :: y
        REAL(rk), ALLOCATABLE :: knots(:), x_k(:)
        REAL(rk) :: x_i, x_j

        ! Determine collocation points (knots)
        ALLOCATE( knots(0:n) )
        knots(:) = get_knots( n )

        x_i = knots(i)
        x_j = knots(j)

        ! Pop x_i
        x_k = PACK( knots, knots/=x_i )

        ! Compute derivative of LAGRANGE polynomial i at collocation point j
        IF ( i==j ) THEN
            y = SUM( 1.0 / ( x_i-x_k(:) ) )
        ELSE
            y = PRODUCT( x_j-x_k(:), MASK=( x_k(:)/=x_j ) ) / &
                PRODUCT( x_i-x_k(:) )
        ENDIF

    END FUNCTION get_dlgll

END MODULE gll_mod
