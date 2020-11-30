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
!! $Date: 2013-11-15 16:19:14 +0100 (Fri, 15 Nov 2013) $
!! $Author: mauerberger $
!! $Revision: 781 $
!! @copyright GNU General Public License version 3 or later



!> Model module
MODULE model_mod
    USE error_mod
    USE parameters_mod
    USE general_mod

    IMPLICIT NONE

    PRIVATE

    TYPE, EXTENDS(general_typ) ::  model_typ
    CONTAINS
        PROCEDURE, NON_OVERRIDABLE :: init_model
        PROCEDURE :: init => init_model
        PROCEDURE :: print => print_model

        PROCEDURE, NON_OVERRIDABLE :: x_glb_min
        PROCEDURE, NON_OVERRIDABLE :: x_glb_max
        PROCEDURE, NON_OVERRIDABLE :: x_loc_min
        PROCEDURE, NON_OVERRIDABLE :: x_loc_max
        PROCEDURE, NON_OVERRIDABLE :: dx_glb
        PROCEDURE, NON_OVERRIDABLE :: dx_loc
        PROCEDURE, NON_OVERRIDABLE :: dx_elm

        PROCEDURE, NON_OVERRIDABLE :: y_glb_min
        PROCEDURE, NON_OVERRIDABLE :: y_glb_max
        PROCEDURE, NON_OVERRIDABLE :: y_loc_min
        PROCEDURE, NON_OVERRIDABLE :: y_loc_max
        PROCEDURE, NON_OVERRIDABLE :: dy_glb
        PROCEDURE, NON_OVERRIDABLE :: dy_loc
        PROCEDURE, NON_OVERRIDABLE :: dy_elm

        PROCEDURE, NON_OVERRIDABLE :: z_glb_min
        PROCEDURE, NON_OVERRIDABLE :: z_glb_max
        PROCEDURE, NON_OVERRIDABLE :: dz_glb
        PROCEDURE, NON_OVERRIDABLE :: dz_loc
        PROCEDURE, NON_OVERRIDABLE :: dz_elm
        PROCEDURE, NON_OVERRIDABLE :: z_loc_min
        PROCEDURE, NON_OVERRIDABLE :: z_loc_max
        PROCEDURE                  :: x_coord_loc
        PROCEDURE                  :: x_coord_glb_r1
        PROCEDURE                  :: y_coord_loc
        PROCEDURE                  :: y_coord_glb_r1
        PROCEDURE                  :: z_coord_loc
        PROCEDURE                  :: z_coord_glb_r1
        PROCEDURE, NON_OVERRIDABLE :: model_type

        PROCEDURE, NON_OVERRIDABLE :: file_name_rhoinv
        PROCEDURE, NON_OVERRIDABLE :: file_name_lambda
        PROCEDURE, NON_OVERRIDABLE :: file_name_mu
        PROCEDURE, NON_OVERRIDABLE :: file_name_a
        PROCEDURE, NON_OVERRIDABLE :: file_name_b
        PROCEDURE, NON_OVERRIDABLE :: file_name_c
        PROCEDURE, NON_OVERRIDABLE :: file_name_q
        PROCEDURE, NON_OVERRIDABLE :: file_name_check
        PROCEDURE, NON_OVERRIDABLE :: override
    END TYPE model_typ

    PUBLIC :: model_typ

CONTAINS


    !> @brief This function returns global minimum/maximum coordinate value in x-direction
    REAL(real_kind) ELEMENTAL FUNCTION x_glb_min( self )
        USE coordinate_utilities_mod, ONLY : deg2rad, lat2colat
        CLASS( model_typ ), INTENT( IN ) :: self
        ! Transform into Ses3d's coordinate system
        x_glb_min = deg2rad( lat2colat( self%lat_max_ ) )
    END FUNCTION x_glb_min
    REAL(real_kind) ELEMENTAL FUNCTION x_glb_max( self )
        USE coordinate_utilities_mod, ONLY : deg2rad, lat2colat
        CLASS( model_typ ), INTENT( IN ) :: self
        ! Transform into Ses3d's coordinate system
        x_glb_max = deg2rad( lat2colat( self%lat_min_ ) )
    END FUNCTION x_glb_max


    !> @brief This function calculates global model extent in x-direction
    REAL(real_kind) ELEMENTAL FUNCTION dx_glb( self )
        CLASS( model_typ ), INTENT( IN ) :: self
        dx_glb = self%x_glb_max() - self%x_glb_min()
    END FUNCTION dx_glb


    !> @brief This function calculates local model extent in x-direction
    REAL(real_kind) ELEMENTAL FUNCTION dx_loc( self )
        CLASS( model_typ ), INTENT( IN ) :: self
        dx_loc = self%dx_glb() / self%px()
    END FUNCTION dx_loc


    !> @brief This function returns element extent in x-direction
    REAL(real_kind) ELEMENTAL FUNCTION dx_elm( self )
        CLASS( model_typ ), INTENT( IN ) :: self
        dx_elm = self%dx_glb() / self%nx_glb()
    END FUNCTION dx_elm



    !> @brief This function returns global minimum/maximum coordinate value in y-direction
    REAL(real_kind) ELEMENTAL FUNCTION y_glb_min( self )
        USE coordinate_utilities_mod, ONLY : deg2rad
        CLASS( model_typ ), INTENT( IN ) :: self
        ! Transform into Ses3d's coordinate system
        y_glb_min = deg2rad( self%lon_min_ )
    END FUNCTION y_glb_min
    REAL(real_kind) ELEMENTAL FUNCTION y_glb_max( self )
        USE coordinate_utilities_mod, ONLY : deg2rad
        CLASS( model_typ ), INTENT( IN ) :: self
        ! Transform into Ses3d's coordinate system
        y_glb_max = deg2rad( self%lon_max_ )
    END FUNCTION y_glb_max


    !> @brief This function calculates lower local x-coordinate
    REAL(real_kind) ELEMENTAL FUNCTION x_loc_min( self )
        CLASS( model_typ ), INTENT( IN ) :: self
        x_loc_min = self%x_glb_min() + self%mpi_x_coord() * self%dx_loc()
    END FUNCTION x_loc_min


    !> @brief This function calculates upper local x-coordinate
    REAL(real_kind) ELEMENTAL FUNCTION x_loc_max( self )
        CLASS( model_typ ), INTENT( IN ) :: self
        x_loc_max = self%x_glb_min() + (self%mpi_x_coord()+1) * self%dx_loc()
    END FUNCTION x_loc_max



    !> @brief This function returns global minimum/maximum coordinate value in z-direction
    REAL(real_kind) ELEMENTAL FUNCTION z_glb_min( self )
        CLASS( model_typ ), INTENT( IN ) :: self
        ! Transform into Ses3d's coordinate system
        z_glb_min = self%z_glb_min_
    END FUNCTION z_glb_min
    REAL(real_kind) ELEMENTAL FUNCTION z_glb_max( self )
        CLASS( model_typ ), INTENT( IN ) :: self
        ! Transform into Ses3d's coordinate system
        z_glb_max = self%z_glb_max_
    END FUNCTION z_glb_max


    !> @brief This function calculates global model extent in y-direction
    REAL(real_kind) ELEMENTAL FUNCTION dy_glb( self )
        CLASS( model_typ ), INTENT( IN ) :: self
        dy_glb = self%y_glb_max() - self%y_glb_min()
    END FUNCTION dy_glb


    !> @brief This function calculates local model extent in y-direction
    REAL(real_kind) ELEMENTAL FUNCTION dy_loc( self )
        CLASS( model_typ ), INTENT( IN ) :: self
        dy_loc = self%dy_glb() / self%py()
    END FUNCTION dy_loc


    !> @brief This function returns element extent in y-direction
    REAL(real_kind) ELEMENTAL FUNCTION dy_elm( self )
        CLASS( model_typ ), INTENT( IN ) :: self
        dy_elm = self%dy_glb() / self%ny_glb()
    END FUNCTION dy_elm


    !> @brief This functionCalculates lower local y-coordinate
    REAL(real_kind) ELEMENTAL FUNCTION y_loc_min( self )
        CLASS( model_typ ), INTENT( IN ) :: self
        y_loc_min = self%y_glb_min() + self%mpi_y_coord() * self%dy_loc()
    END FUNCTION y_loc_min


    !> @brief This function calculates upper local y-coordinate
    REAL(real_kind) ELEMENTAL FUNCTION y_loc_max( self )
        CLASS( model_typ ), INTENT( IN ) :: self
        y_loc_max = self%y_glb_min() + (self%mpi_y_coord()+1) * self%dy_loc()
    END FUNCTION y_loc_max




    !> @brief This function calculates global model extent in z-direction
    REAL(real_kind) ELEMENTAL FUNCTION dz_glb( self )
        CLASS( model_typ ), INTENT( IN ) :: self
        dz_glb = self%z_glb_max_ - self%z_glb_min_
    END FUNCTION dz_glb


    !> @brief This function calculates local model extent in z-direction
    REAL(real_kind) ELEMENTAL FUNCTION dz_loc( self )
        CLASS( model_typ ), INTENT( IN ) :: self
        dz_loc = self%dz_glb() / self%pz()
    END FUNCTION dz_loc


    !> @brief This function returns element extent in z-direction
    REAL(real_kind) ELEMENTAL FUNCTION dz_elm( self )
        CLASS( model_typ ), INTENT( IN ) :: self
        dz_elm = self%dz_glb() / self%nz_glb()
    END FUNCTION dz_elm


    !> @brief This function calculates lower local z-coordinate
    REAL(real_kind) ELEMENTAL FUNCTION z_loc_min( self )
        CLASS( model_typ ), INTENT( IN ) :: self
        z_loc_min = self%z_glb_min_ + self%mpi_z_coord() * self%dz_loc()
    END FUNCTION z_loc_min


    !> @brief This function calculates upper local z-coordinate
    REAL(real_kind) ELEMENTAL FUNCTION z_loc_max( self )
        CLASS( model_typ ), INTENT( IN ) :: self
        z_loc_max = self%z_glb_min_ + (self%mpi_z_coord()+1) * self%dz_loc()
    END FUNCTION z_loc_max


    !> @brief This function calculates local x coordinate
    PURE FUNCTION x_coord_loc( self ) RESULT ( x )
        USE gll_mod, ONLY : get_knots
        CLASS( model_typ ), INTENT( IN ) :: self
        INTEGER :: i
        REAL(real_kind) :: x(0:self%nx_loc(),0:self%lpd()), &
                           knots(0:self%lpd())
        knots(:) = get_knots(self%lpd())
        DO i=0, self%nx_loc()
            x(i,:) = self%x_loc_min() &
                + (i)*self%dx_elm() + 0.5*(1.0+knots)*self%dx_elm()
        ENDDO
    END FUNCTION x_coord_loc


    !> @brief This function reshape local x coordinate to r1
    PURE FUNCTION x_coord_glb_r1( self ) RESULT ( x_r1 )
        USE gll_mod, ONLY : get_knots
        CLASS( model_typ ), INTENT( IN ) :: self
        INTEGER :: i
        REAL(real_kind) :: x(0:self%nx_glb()-1,0:self%lpd()), &
                           knots(0:self%lpd()), &
                           x_r1( self%nx_glb()*self%lpd()+1 )
        knots(:) = get_knots(self%lpd())
        DO i=0, self%nx_glb() - 1
            x(i,:) = self%x_glb_min() &
                + i*self%dx_elm() + 0.5*(1.0+knots(:))*self%dx_elm()
        ENDDO

        x_r1(:) = pack_r2tor1( x(:,:) )

    END FUNCTION x_coord_glb_r1



    !> @brief This function calculates local y coordinate
    PURE FUNCTION y_coord_loc( self ) RESULT ( y )
        USE gll_mod, ONLY : get_knots
        CLASS( model_typ ), INTENT( IN ) :: self
        INTEGER :: i
        REAL(real_kind) :: y(0:self%ny_loc(),0:self%lpd()), &
                           knots(0:self%lpd())
        knots(:) = get_knots(self%lpd())
        DO i=0, self%ny_loc()
            y(i,:) = self%y_loc_min() &
                + i*self%dy_elm() + 0.5*(1.0+knots(:))*self%dy_elm()
        ENDDO
    END FUNCTION y_coord_loc


    !> @brief This function  reshape local y coordinate to r1
    PURE FUNCTION y_coord_glb_r1( self ) RESULT ( y_r1 )
        USE gll_mod, ONLY : get_knots
        CLASS( model_typ ), INTENT( IN ) :: self
        INTEGER :: i
        REAL(real_kind) :: y(0:self%ny_glb()-1,0:self%lpd()), &
                           knots(0:self%lpd()), &
                           y_r1( self%ny_glb()*self%lpd()+1 )
        knots(:) = get_knots(self%lpd())
        DO i=0, self%ny_glb()-1
            y(i,:) = self%y_glb_min() &
                + i*self%dy_elm() + 0.5*(1.0+knots(:))*self%dy_elm()
        ENDDO

        y_r1(:) = pack_r2tor1( y(:,:) )

    END FUNCTION y_coord_glb_r1


    !> @brief This function calculates local z coordinate
    PURE FUNCTION z_coord_loc( self ) RESULT ( z )
        USE gll_mod, ONLY : get_knots
        CLASS( model_typ ), INTENT( IN ) :: self
        INTEGER :: i
        REAL(real_kind) :: z(0:self%nz_loc(),0:self%lpd()), &
                           knots(0:self%lpd())
        knots(:) = get_knots(self%lpd())
        ! Why the hell are z-coordinates upside down?
        DO i=0, self%nz_loc()
            z(i,:) = self%z_loc_max() &
                - i*self%dz_elm() - 0.5*(1.0+knots(:))*self%dz_elm()
        ENDDO
    END FUNCTION z_coord_loc


    !> @brief This function reshape local z coordinate to r1
    PURE FUNCTION z_coord_glb_r1( self ) RESULT ( z_r1 )
        USE gll_mod, ONLY : get_knots
        CLASS( model_typ ), INTENT( IN ) :: self
        INTEGER :: i
        REAL(real_kind) :: z(0:self%nz_glb()-1,0:self%lpd()), &
                           knots(0:self%lpd()), &
                           z_r1( self%nz_glb()*self%lpd()+1 )
        knots(:) = get_knots(self%lpd())
        ! Why the hell are z-coordinates upside down?
        DO i=0, self%nz_glb() - 1
            z(i,:) = self%z_glb_max() &
                - i*self%dz_elm() - 0.5*(1.0+knots(:))*self%dz_elm()
        ENDDO

        z_r1(:) = pack_r2tor1( z(:,:) )

    END FUNCTION z_coord_glb_r1


    !> @brief This function returns model-type
    ELEMENTAL FUNCTION model_type( self )
        CLASS( model_typ ), INTENT( IN ) :: self
        INTEGER :: model_type
        model_type = self%model_type_
    END FUNCTION model_type



    !> @brief This function trimmed model-parameter file-names are returned
    PURE FUNCTION file_name_rhoinv( self )
        CLASS( model_typ ), INTENT( IN ) :: self
        CHARACTER( LEN=LEN_TRIM( self%rhoinv_ ) ) :: file_name_rhoinv
        file_name_rhoinv = self%rhoinv_
    END FUNCTION file_name_rhoinv
    PURE FUNCTION file_name_lambda( self )
        CLASS( model_typ ), INTENT( IN ) :: self
        CHARACTER( LEN=LEN_TRIM( self%lambda_ ) ) :: file_name_lambda
        file_name_lambda = self%lambda_
    END FUNCTION file_name_lambda
    PURE FUNCTION file_name_mu( self )
        CLASS( model_typ ), INTENT( IN ) :: self
        CHARACTER( LEN=LEN_TRIM( self%mu_ ) ) :: file_name_mu
        file_name_mu = self%mu_
    END FUNCTION file_name_mu
    PURE FUNCTION file_name_a( self )
        CLASS( model_typ ), INTENT( IN ) :: self
        CHARACTER( LEN=LEN_TRIM( self%a_ ) ) :: file_name_a
        file_name_a = self%a_
    END FUNCTION file_name_a
    PURE FUNCTION file_name_b( self )
        CLASS( model_typ ), INTENT( IN ) :: self
        CHARACTER( LEN=LEN_TRIM( self%b_ ) ) :: file_name_b
        file_name_b = self%b_
    END FUNCTION file_name_b
    PURE FUNCTION file_name_c( self )
        CLASS( model_typ ), INTENT( IN ) :: self
        CHARACTER( LEN=LEN_TRIM( self%c_ ) ) :: file_name_c
        file_name_c = self%c_
    END FUNCTION file_name_c
    PURE FUNCTION file_name_q( self )
        CLASS( model_typ ), INTENT( IN ) :: self
        CHARACTER( LEN=LEN_TRIM( self%a_ ) ) :: file_name_q
        file_name_q = self%q_
    END FUNCTION file_name_q

    PURE LOGICAL FUNCTION override( self )
        CLASS( model_typ ), INTENT( IN ) :: self
        override = self%override_
    END FUNCTION override


    !> @brief This subroutine for printing model
    SUBROUTINE print_model( self, unit )
        CLASS(model_typ), INTENT(IN) :: self
        INTEGER, INTENT(IN) :: unit
        CHARACTER(LEN=fsl) :: fmt

        fmt = '(A)'
        WRITE( UNIT=unit, FMT=fmt ) 'Model:'

        fmt = '( 4x, 2((A," = ",G0):,",  " ) )'
        WRITE( UNIT=unit, FMT=fmt ) 'x_glb_min', self%x_glb_min(),&
                                    'x_glb_max', self%x_glb_max()
        WRITE( UNIT=unit, FMT=fmt ) 'y_glb_min', self%y_glb_min(),&
                                    'y_glb_max', self%y_glb_max()
        WRITE( UNIT=unit, FMT=fmt ) 'z_glb_min', self%z_glb_min_,&
                                    'z_glb_max', self%z_glb_max_

        fmt = '(4x,A,A)'
        WRITE( UNIT=unit, FMT=fmt ) 'rhoinv = ', self%file_name_rhoinv( )
        WRITE( UNIT=unit, FMT=fmt ) 'lambda = ', self%file_name_lambda( )
        WRITE( UNIT=unit, FMT=fmt ) 'mu     = ', self%file_name_mu( )
        WRITE( UNIT=unit, FMT=fmt ) 'a      = ', self%file_name_a( )
        WRITE( UNIT=unit, FMT=fmt ) 'b      = ', self%file_name_b( )
        WRITE( UNIT=unit, FMT=fmt ) 'c      = ', self%file_name_c( )
        WRITE( UNIT=unit, FMT=fmt ) 'q      = ', self%file_name_q( )

        fmt = '( 4x, 2((A," = ",G0):,",  " ) )'
        WRITE( UNIT=unit, FMT=fmt ) 'x_loc_min', self%x_loc_min(), &
                                    'x_loc_max', self%x_loc_max()
        WRITE( UNIT=unit, FMT=fmt ) 'y_loc_min', self%y_loc_min(), &
                                    'y_loc_max', self%y_loc_max()
        WRITE( UNIT=unit, FMT=fmt ) 'z_loc_min', self%z_loc_min(), &
                                    'z_loc_max', self%z_loc_max()

        fmt = '( 4x, 3((A," = ",G0):,",  " ) )'
        WRITE( UNIT=unit, FMT=fmt ) 'dx_elm', self%dx_elm(), &
                                    'dy_elm', self%dy_elm(), &
                                    'dz_elm', self%dz_elm()

    END SUBROUTINE print_model




    !> @brief This subroutine  initializes to Ses3d's coordinate system
    SUBROUTINE init_model( self )
        CLASS(model_typ), INTENT(INOUT) :: self

        ! Call parent init procedure first
        CALL self%init_general()

        IF ( self%workflow() == 'model' ) THEN
            IF ( self%model_type_ == -1 ) &
                CALL abort( "In workflow 'model' a model_type must different&
                           & from -1 must be specified! " )
            IF ( self%file_name_rhoinv() == 'N/A' .AND. &
                 self%file_name_lambda() == 'N/A' .AND. &
                 self%file_name_mu()     == 'N/A' .AND. &
                 self%file_name_a()      == 'N/A' .AND. &
                 self%file_name_b()      == 'N/A' .AND. &
                 self%file_name_c()      == 'N/A' .AND. &
                 self%file_name_q()      == 'N/A' ) &
                CALL abort( "In workflow 'model' at least one model-parameter&
                           &-filename must be specified! " )
        ELSE
            IF ( self%model_type_ == -1 ) THEN
                ! Check if rho, lambda and mu are specified
                IF ( self%file_name_rhoinv() == 'N/A' ) &
                    CALL abort( '"rhoinv" model-file must be given!' )
                IF ( self%file_name_lambda() == 'N/A' ) &
                    CALL abort( '"lambda" model-file must be given!' )
                IF ( self%file_name_mu() == 'N/A' ) &
                    CALL abort( '"mu" model-file must be given!' )
                ! TODO the following should be done at MPI_FILE_OPEN()
                CALL self%file_name_check()
            ELSE
                IF ( self%file_name_rhoinv() /= 'N/A' ) &
                    CALL abort( "'rhoinv' and a 'model_type' can't be specified&
                                & at the same time")
                IF ( self%file_name_lambda() /= 'N/A' ) &
                    CALL abort( "'lambda' and a 'model_type' can't be specified&
                                &at the same time")
                IF ( self%file_name_mu() /= 'N/A' ) &
                    CALL abort( "'mu' and a 'model_type' can't be specified&
                                &at the same time")
                ! TODO how to set a,b,c,q as optional
                !IF ( self%file_name_a() == 'N/A' ) &
                !    CALL abort( "'a' and a 'model_type' can't be specified&
                !                &at the same time")
                !IF ( self%file_name_b() == 'N/A' ) &
                !    CALL abort( "'b' and a 'model_type' can't be specified&
                !                &at the same time")
                !IF ( self%file_name_c() == 'N/A' ) &
                !    CALL abort( "'c' and a 'model_type' can't be specified&
                !                &at the same time")
                !IF ( self%file_name_q() == 'N/A' ) &
                !    CALL abort( "'q' and a 'model_type' can't be specified&
                !                &at the same time")

            END IF
        END IF

    END SUBROUTINE init_model



    !> @brief This subroutine checks if model-parameter files exist
    SUBROUTINE file_name_check( self )
        USE checks_mod, ONLY : file_exists
        CLASS(model_typ), INTENT(IN) :: self

        ! Check if rohinv-file exists
        IF ( .NOT. file_exists( self%rhoinv_ ) ) &
            CALL abort( 'Model-file "' //self%file_name_rhoinv()// &
                '" does not exist!' )

        ! Check if lambda-file exists
        IF ( .NOT. file_exists( self%lambda_ ) ) &
            CALL abort( 'Model-file "' //self%file_name_lambda()// &
                '" does not exist!' )

        ! Check if mu-file exists
        IF ( .NOT. file_exists( self%mu_ ) ) &
            CALL abort( 'Model-file "' //self%file_name_mu()// &
                '" does not exist!' )

        ! In case a model-file is specified
        IF ( self%a_ /= 'N/A' ) THEN
            ! Check if a model-file exists
            IF ( .NOT. file_exists( self%a_ ) ) &
                CALL abort( 'Model-file "' //self%file_name_a()// &
                    '" does not exist!' )
        END IF

        ! In case b model-file is specified
        IF ( self%b_ /= 'N/A' ) THEN
            ! Check if b model-file exists
            IF ( .NOT. file_exists( self%b_ ) ) &
                CALL abort( 'Model-file "' //self%file_name_b()// &
                    '" does not exist!' )
        END IF

        ! In case c model-file is specified
        IF ( self%c_ /= 'N/A' ) THEN
            ! Check if a model-file exists
            IF ( .NOT. file_exists( self%c_ ) ) &
                CALL abort( 'Model-file "' //self%file_name_c()// &
                    '" does not exist!' )
        END IF

        ! In case q model-file is specified
        IF ( self%q_ /= 'N/A' ) THEN
            ! Check if q model-file exists
            IF ( .NOT. file_exists( self%q_ ) ) &
                CALL abort( 'Model-file "' //self%file_name_q()// &
                    '" does not exist!' )
        END IF

    END SUBROUTINE file_name_check



    !> @brief This function packs array of rank two into an array of rank one
    !! with dropping duplicates
    PURE FUNCTION pack_r2tor1( r2 ) RESULT( r1 )
        REAL(real_kind), INTENT(IN) :: r2(:,:)
        REAL(real_kind), ALLOCATABLE :: r1(:)
        INTEGER :: nx, lpd

        lpd = SIZE( r2, DIM=2 )
        nx = SIZE( r2, DIM=1 )

        r1 = [ TRANSPOSE( r2(:,:lpd-1) ), r2(nx,lpd) ]

    END FUNCTION pack_r2tor1


END MODULE model_mod
