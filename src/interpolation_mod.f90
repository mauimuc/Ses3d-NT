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
!! $Date: 2013-11-03 22:24:38 +0100 (Sun, 03 Nov 2013) $
!! $Author: mauerberger $
!! $Revision: 764 $
!! @copyright GNU General Public License version 3 or later


!> This module provides several routines for Lagrange-Interpolation 
!! 
!! \f[ f(x) \approx \sum_{i=0}^N f_i \ell^N_i(x) \f]
!!
!! The module can be used as follows:
!! 1. Start with a set of N+1 disjoint collocation points named `knots(:)`
!! 2. Pre-calculate the N+1 Lagrange-polynomial's denominators 
!!    `weights(:) = calc_weights(knots)`
!! 3. Pre-calculate the N+1 Lagrange-polynomials for a pint x
!!    `l(:) = calc_lagrange_polynomial(x, knots, weights)`
!! 4. To carry out the actual interpolation a set of N+1 function values 
!! \f$f_i=f(x_i)\f$ is needed. An implementation might be `SUM( l(:)*f(:) )`.
!!
!! In addition there is a convenience function create_3d_vandermonde_matrix() 
!! which simplifies 3d Lagrange interpolation for an orthogonal set of
!! Lagrange polynomials e.g. l_x(:), l_y(:) and l_z(:). 
!!
!! @todo interpolation_test.f90
!! 
!! To specify the precision of floats this module expects to find a parameter 
!! parameters_mod::real_kind in a module named parameters_mod. 
MODULE interp_mod
    USE parameters_mod, ONLY : rk => real_kind

    IMPLICIT NONE

    PRIVATE

    !> Distinguishes if a pre-processed denominator `weights(:)` for Lagrange 
    !! interpolation is provided or not. 
    !!
    !! It may be called by 
    !!  * calc_lagrange_polynomial(x, knots(:))
    !!  * calc_lagrange_polynomial(x, knots(:), weights(:))
    !!
    !! For a detailed description see lagrange_p_kw() or lagrange_p_k()
    INTERFACE calc_lagrange_polynomial 
        MODULE PROCEDURE lagrange_p_kw
        MODULE PROCEDURE lagrange_p_k
    END INTERFACE calc_lagrange_polynomial 
    
    PUBLIC :: calc_lagrange_polynomial, calc_weights, create_3d_vandermonde_matrix

CONTAINS 

    !> Pre-calculates the Lagrange polynomial's denominators `weights(:)`
    !!
    !! This usually is a preprocessing step for optimization.
    !! For degree N Lagrange polynomials the denominator is defined as
    !! \f[
    !!  w_i = \left ( \prod_{i=0}^{N} (x_j-x_i) \right )^{-1} 
    !!      \;\;\; \text{for} \; i \neq j
    !! \f]
    !! with \f$x_i\f$ the N+1 collocation points
    !!
    !! @param knots The N+1 collocation points \f$x_i\f$. Array of size N+1.
    !! @result The Lagrange polynomial's denominators \f$w_i\f$. 
    !!         Array of size N+1. 
    PURE FUNCTION calc_weights(knots) RESULT(weights)
        REAL(rk), INTENT(IN) :: knots(:)
        REAL(rk) :: weights(SIZE(knots))
        INTEGER :: j

        weights(:) = 1.0 / [ ( ( PRODUCT( knots(j)-knots(:j-1) ) * &
                                 PRODUCT( knots(j)-knots(j+1:) ) ), &
                                 j=1, SIZE( knots ) ) ]

    END FUNCTION calc_weights


    !> This is just a convenience function taking together calc_weights() and
    !! lagrange_p_kw(). 
    !!
    !! @waring Do not use this function for evaluating the same polynomial for 
    !!         multiple x. In the interest of performance provide pre-
    !!         calculated denominators i.e. `weights(:)`. 
    !!
    !! To perform an actual interpolation N+1 function values \f$f_i=f(x_i)\f$ 
    !! are needed. An implementation might be `SUM( f(:)*l(:) )`
    !! 
    !! @param x Interpolation point; Should be within the range of collocation 
    !!          points \f$[x_0 \dots x_N]\f$.
    !! @param knots The N+1 disjoint collocation points \f$x_i\f$. Array 
    !!              of size N+1.
    !! @result The Lagrange polynomial \f$(\ell_0^N(x), \dots, \ell_N^N(x)) \f$ 
    !!         evaluated for x. Array of size N+1.
    PURE FUNCTION lagrange_p_k(x, knots) RESULT(l)
        REAL(rk), INTENT(IN) :: x, knots(:)
        REAL(rk) :: l(SIZE(knots))
        
        l(:) = lagrange_p_kw(x=x, knots=knots, weights=calc_weights(knots) )

    END FUNCTION 

    !> Evaluates the Lagrange-polynomial for a single point x. To obtain the 
    !! `weigts(:)` - the Lagrange polynomials denominator - invoke 
    !! calc_weights(). 
    !!
    !! In terms of efficiency the following approach makes more sense than 
    !! implementing the pure definition:
    !! \f{align*}{
    !!  \ell^N_j(x) 
    !!   :&= \prod_{i=0}^{N} \frac{x-x_i}{x_j-x_i}
    !!    \;\;\;\text{for}\; i \neq j \\
    !!   & = \frac{1}{w_j} \prod_{i=0}^{N} (x-x_i)
    !!    \;\;\;\text{for}\; i \neq j \\
    !!   & = \frac{1}{w_j} \; \frac{1}{x-x_j} \prod_{i=0}^{N} (x-x_i) 
    !!    \;\;\;\text{for}\; x \neq x_j \\
    !! \f} 
    !! 
    !! For actual computation the property \f$ \ell^N_i(x_j) = \delta_{ij} \f$
    !! is used. 
    !!
    !! @note If just a single collocation point is passed `[1.0]` is returned. 
    !!
    !! To perform an actual interpolation N+1 function values \f$f_i=f(x_i)\f$ 
    !! are needed. An implementation might be `SUM( f(:)*l(:) )`
    !! 
    !! @param x Interpolation point. Should lie within the range of collocation 
    !!          points \f$[x_0 \dots x_N]\f$.
    !! @param knots The N+1 disjoint collocation points \f$x_i\f$. Array 
    !!              of size N+1.
    !! @param weights Lagrange polynomial's denominators \f$w_i\f$. Array 
    !!        of size N+1. See calc_weights().
    !! @result The Lagrange polynomial \f$(\ell_0^N(x), \dots, \ell_N^N(x)) \f$ 
    !!         evaluated for x. Array of size N+1.
    PURE FUNCTION lagrange_p_kw(x, knots, weights) RESULT(l)
        REAL(rk), INTENT(IN) :: x, knots(:), weights(:)
        REAL(rk) :: l(SIZE(knots))
        LOGICAL :: mask(SIZE(knots))

        IF ( SIZE( knots ) .LE. 1 ) THEN ! only one sampling point
            l(:) = 1.0
            RETURN
        END IF

        mask(:) = x .NE. knots
        IF ( ALL( mask ) ) THEN ! if collocation-p. != interpolation-p.
            l(:) = PRODUCT( x - knots(:) ) / ( x - knots(:) ) * weights
            RETURN
        ELSE ! if interpolation-p. matches with one collocation-p.
        ! uses the property l_j(x_i) = delta_ij
            l(:) = MERGE( 0.0_rk, 1.0_rk, mask ) 
            RETURN
        END IF

    END FUNCTION 


    !> Convenience function which calculates the so called Vandermonde-matrix 
    !! corresponding to a single interpolation point (x,y,z). Actually, this is
    !! not doing much more than a reshape.
    !!
    !! \f[
    !! v_{ijk} = \ell^{N_x}_i(x) \ell^{N_y}_j(y) \ell^{N_z}_k(z)
    !! \f]
    !!
    !! To carry out the actual interpolation for a set of \f$f_{ijk} = 
    !! f(x_i,y_j,z_k)\f$ function values - which is \f$ f(x,y,z) \approx 
    !! \sum_{i,j,k=0}^{N_x,N_y,N_z} f_{ijk} v_{ijk} \f$ - 
    !! say `SUM( vm(:,:,:)*f(:,:,:) )`. 
    !!
    !! @param l_x Lagrange polynomial in x-direction evaluated for x. 
    !!!           Size N_x + 1.
    !! @param l_y Lagrange polynomial in y-direction evaluated for y. 
    !!            Size N_y + 1.
    !! @param l_z Lagrange polynomial in z-direction evaluated for z. 
    !!            Size N_z + 1.
    !! @result The Vandermonde-matrix of shape (N_x+1, N_y+1, N_z+1)
    PURE FUNCTION create_3d_vandermonde_matrix(l_x, l_y, l_z) RESULT(vm)
        REAL(rk), INTENT(IN) :: l_x(:), l_y(:), l_z(:)
        REAL(rk) :: vm(SIZE(l_x), SIZE(l_y), SIZE(l_z))
        INTEGER :: shp(3)

        shp(:) = [ SIZE(l_x), SIZE(l_y), SIZE(l_z) ]

        vm(:,:,:) = RESHAPE(l_x, shp, PAD=l_x, ORDER=[1,2,3])* &
                    RESHAPE(l_y, shp, PAD=l_y, ORDER=[2,1,3])* &
                    RESHAPE(l_z, shp, PAD=l_z, ORDER=[3,2,1])

    END FUNCTION create_3d_vandermonde_matrix

END MODULE interp_mod
