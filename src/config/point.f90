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
!! $Date: 2013-12-04 08:50:52 +0100 (Wed, 04 Dec 2013) $
!! $Author: mauerberger $
!! $Revision: 838 $
!! @copyright GNU General Public License version 3 or later



!> Module providing a class (point_cls) to handle points (e.g. receivers)
MODULE point_mod
    USE parameters_mod, ONLY : real_kind
    USE configuration_mod, ONLY : conf => configuration

    IMPLICIT NONE

    PRIVATE

    !> A class which handles points in Ses3d-NT. It is used as a base-class
    !! for receivers and sources.
    !!
    !! It provides the following procedures:
    !!
    !!  * Local element indices: ex, ey and ez
    !!  * Physical coordinates: latitude, longitude and dept
    !!  * The Vandermonde-matrix dedicated to interpolation
    !!  * A print procedure
    !!
    !! There are three structure constructor functions:
    !!
    !!  * Fortran's default constructor expecting arguments ex, ey, ez and vm
    !!  * A constructor expecting no arguments. It returns an uninitialized
    !!    object.
    !!  * The actual used constructor expecting lat, lon, depth as arguments.
    !!    If something goes wrong an uninitialized point_cls object is returned.
    TYPE :: point_cls
        !> Local element index in x/theta direction
        INTEGER, PRIVATE :: ex
        !> Local element index in y/phi direction
        INTEGER, PRIVATE :: ey
        !> Local element index in z/radial direction
        INTEGER, PRIVATE :: ez
        !> Vandermonde matrix used for making interpolation more efficient
        REAL(real_kind), ALLOCATABLE, PRIVATE :: vm(:,:,:)
    CONTAINS
        !> Returns the latitude, obtained through interpolation.
        !! For details see point_mod::get_lat().
        PROCEDURE, PUBLIC :: lat => get_lat
        !> Returns the longitude, obtained through interpolation.
        !! For details see point_mod::get_lon().
        PROCEDURE, PUBLIC :: lon => get_lon
        !> Returns the depth, obtained through interpolation.
        !! For details see point_mod::get_depth().
        PROCEDURE, PUBLIC :: depth => get_depth
        !> Returns the Vandermond-matrix.
        !! For details see point_mod::get_vanderm().
        PROCEDURE, PUBLIC :: vanderm => get_vanderm
        !> Returns the local x/theta element index the point resides in.
        !! For details see point_mod::get_ex().
        PROCEDURE, PUBLIC :: elm_x_loc  => get_ex
        !> Returns the local y/phi element index the point resides in.
        !! For details see point_mod::get_ey().
        PROCEDURE, PUBLIC :: elm_y_loc  => get_ey
        !> Returns the local z/radial element index the point resides in.
        !! For details see point_mod::get_ez().
        PROCEDURE, PUBLIC :: elm_z_loc  => get_ez
        !> Returns .TRUE. if the point object is not initialized.
        !! @warning 'Uninitialized' just says that the object is precisely the
        !!          same as point_mod::uninitialized_point_cls().
        !! For details see point_mod::uninitialized()
        PROCEDURE, PUBLIC :: uninitialized
        !> Writes some information about the object into a unit specifier.
        !! For details see point_mod::print_point().
        PROCEDURE, PUBLIC :: print => print_point
        !> For comparing two point_cls objects. Do not use this bound procedure
        !! explicitly. Use the comparison operator '==' instead.
        !! For details see point_mod::equal_points().
        PROCEDURE, PRIVATE :: equal => equal_points
        !> Defines the comparison operator '==' for point_cls class and any
        !! extension. For details see point_mod::equal_points().
        GENERIC, PUBLIC :: OPERATOR(==) => equal
    END TYPE

    ! Extends Fortran's default structure constructor
    INTERFACE point_cls
        !> Constructs a point_cls object from passed values lat, lon and depth.
        !! For details see point_mod::construct_point_cls(lat, lon, depth)
        MODULE PROCEDURE :: construct_point_cls
        !> A constructor expecting no arguments and returning an uninitialized
        !! point_cls object. That constructor is used if something goes wrong.
        !! For details see point_mod::uninitialized_point_cls()
        MODULE PROCEDURE :: uninitialized_point_cls
    END INTERFACE

    ! Expose the class with its constructors only
    PUBLIC :: point_cls

CONTAINS

!------------------------------------------------------------------------------
! Structure constructors

    !> Structure constructor returning an uninitialized/empty point_cls object
    !!
    !! @note An object of class point_cls is called uninitialized if it is
    !!       instantiated using that very constructor.
    !!
    !! @return Uninitialized object of class point_cls
    ELEMENTAL FUNCTION uninitialized_point_cls() RESULT(new)
        ! Result variable
        TYPE(point_cls) :: new
        ! Local variables
        REAL(real_kind) :: vm(0,0,0)

        ! Invoke Fortran's default constructor
        new = point_cls(ex=HUGE(0), ey=HUGE(0), ez=HUGE(0), vm=vm)

    END FUNCTION


    !> Constructs a point_cls object by passing latitude, longitude and dept.
    !!
    !! It covers a whole bunch of consistency checks as well as several pre-
    !! processing steps for making interpolation more efficient.
    !!
    !! @note If something goes wrong an uninitialized point_cls object is
    !!       returned (see point_mod::uninitialized_point_cls() ).
    !!
    !! @todo Write warnings to log-files
    !!
    !! @param lat Latitude in [deg]
    !! @param lat Longitude in [deg]
    !! @param depth Depth in [m] (counting from parameters_mod::EARTH_RADIUS)
    FUNCTION construct_point_cls(lat, lon, depth) RESULT(new)
        ! Use associated entities
        USE parameters_mod, ONLY : is_mpi_master_rank
        USE error_mod, ONLY : warn
        USE checks_mod, ONLY : valid_lat, valid_lon, valid_depth
        USE interp_mod, ONLY : calc_weights, calc_lagrange_polynomial, &
                               create_3d_vandermonde_matrix
        USE gll_mod, ONLY : get_knots, get_weights
        USE coordinate_utilities_mod, ONLY : deg2rad, lat2colat, depth2radius
        ! Passed dummy arguments
        REAL(real_kind), INTENT(IN) :: lat, lon, depth
        ! Result variable
        TYPE(point_cls) :: new
        ! Local variables
        REAL(real_kind), DIMENSION(0:conf%lpd()) :: l_x_, l_y_, l_z_, knots_
        REAL(real_kind) :: x_, y_, z_
        CHARACTER(LEN=80) :: str_

        ! Check latitude
        IF ( .NOT. valid_lat(lat) ) THEN
            WRITE(str_, FMT='("Got invalid latitude! lat=", G0 )') lat
            ! Print to stdout
            IF ( is_mpi_master_rank() ) &
                CALL warn( TRIM(str_) )
            ! todo log-file
            ! Return uninitialized point_cls object
            new = point_cls()
            RETURN
        END IF
        ! Check longitude
        IF ( .NOT. valid_lon(lon) ) THEN
            WRITE(str_, FMT='("Got invalid longitue! lon=", G0 )') lon
            ! Print to stdout
            IF ( is_mpi_master_rank() ) &
                CALL warn( TRIM(str_) )
            ! todo log-file
            ! Return uninitialized point_cls object
            new = point_cls()
            RETURN
        END IF
        ! Check depth
        IF ( .NOT. valid_depth(depth) ) THEN
            WRITE(str_, FMT='("Got invalid depth! dept=", G0 )') depth
            ! Print to stdout
            IF ( is_mpi_master_rank() ) &
                CALL warn( TRIM(str_) )
            ! todo log-file
            ! Return uninitialized point_cls object
            new = point_cls()
            RETURN
        END IF

        ! Transform into Ses3d-NT's coordinate system
        x_ = deg2rad(lat2colat(lat))
        y_ = deg2rad(lon)
        z_ = depth2radius(depth)

        ! Check if in domain
        IF ( .NOT. in_domain(x=x_, y=y_, z=z_) ) THEN
            WRITE(str_, FMT='("Point (", G0, ", ", G0, ", ", G0, ") not in &
                &computational domain!" )') lat, lon, depth
            ! Print to stdout
            IF ( is_mpi_master_rank() ) &
                CALL warn( TRIM(str_) )
            ! todo log-file
            ! Return uninitialized point_cls object
            new = point_cls()
            RETURN
        END IF

        ! Check is point is in domain
        IF ( .NOT. in_rank(x=x_, y=y_, z=z_) ) THEN
            ! Return uninitialized point_cls object
            new = point_cls()
            RETURN
        END IF

        ! From here on only one MPI rank remains all the others are sorted out

        ! Check if in tapering boundaries
        IF ( in_boundary(x=x_, y=y_, z=z_) ) THEN
            WRITE(str_, FMT='("Point (", G0, ", ", G0, ", ", G0, ") located in &
                &tapering boundaries!")') lat, lon, depth
            ! Print to stdout
            CALL warn( TRIM(str_) )
            ! todo log-file
        END IF

        ! The following pre-processing is due to optimizing Lagrange-
        ! interpolation which is necessary at every time-step
        ! Retrieve GLL-collocation points
        knots_(:) = get_knots( conf%lpd() )
        ! Calculate Lagrange-polynomials at (x,y,z)
        l_x_(:) = calc_lagrange_polynomial( std_x(x=x_), knots_ )
        l_y_(:) = calc_lagrange_polynomial( std_y(y=y_), knots_ )
        l_z_(:) = calc_lagrange_polynomial( std_z(z=z_), knots_ )


        ! Finally instantiate new point_cls object. This covers:
        !  * Determine element indices
        !  * Pre-process Vandermonde-matrix for latter interpolation
        new = point_cls(ex=elm_x_loc(x_), ey=elm_y_loc(y_), ez=elm_z_loc(z_), &
                        vm=create_3d_vandermonde_matrix(l_x_, l_y_, l_z_) )

    END FUNCTION


!------------------------------------------------------------------------------
! source_cls type bound procedures


    !> Writes various information about a point - a point_cls object
    !!
    !! If unit is not specified parameters_mod::output_unit is assumed by default
    !!
    !! @todo get rid of the optional attribute
    !!
    !! @param self Passed-object dummy of class point_cls
    !! @param unit A unit specifier
    SUBROUTINE print_point(self, unit)
        ! Passed dummy arguments
        CLASS(point_cls), INTENT(IN) :: self
        INTEGER, INTENT(IN) :: unit
        ! Local variables
        CHARACTER(LEN=*), PARAMETER :: fmt = &
            '( 4x,  3( A, " = ", (G0,:,3x) )  )' ! Formatter string

        ! Source location
        WRITE(UNIT=unit, FMT=fmt) 'lat',   get_lat(self), &
                                  'lon',   get_lon(self), &
                                  'depth', get_depth(self)
        ! Ses3d's coordinates
        WRITE(UNIT=unit, FMT=fmt) 'x', interpolate_x(self), &
                                  'y', interpolate_y(self), &
                                  'z', interpolate_z(self)
        ! Coordinates in unit-element
        WRITE(UNIT=unit, FMT=fmt) 'std_x', std_x(interpolate_x(self)), &
                                  'std_y', std_y(interpolate_y(self)), &
                                  'std_z', std_z(interpolate_z(self))
        ! Local element indices
        WRITE(UNIT=unit, FMT=fmt) 'elm_x_loc', self%ex, &
                                  'elm_y_loc', self%ey, &
                                  'elm_z_loc', self%ez

    END SUBROUTINE print_point


    !> Returns the x/theta element index of that point
    !!
    !! To access the element a point is located in all three element indices
    !! \f$e_\vartheta\f$, \f$e_\varphi\f$ and \f$e_r\f$ are needed. To actually
    !! select values from one of Ses3d-NT's rank 6 fields say e.g.
    !!
    !!     f_loc(p%elm_x_loc(),p%elm_y_loc(),p%elm_z_loc(),:,:,:)
    !!
    !! where p is an object of class point_cls. This returns an array of shape
    !! [0:lpd,0:lpd,0:lpd].
    !!
    !! @param self Passed-object dummy of class point_cls
    !! @return The x/theta element index
    ELEMENTAL INTEGER FUNCTION get_ex(self)
        ! Passed dummy arguments
        CLASS(point_cls), INTENT(IN) :: self

        get_ex = self%ex

    END FUNCTION


    !> Returns the y/phi element index of that point
    !!
    !! To access the element a point is located in all three element indices
    !! \f$e_\vartheta\f$, \f$e_\varphi\f$ and \f$e_r\f$ are needed. To actually
    !! select values from one of Ses3d-NT's rank 6 fields say e.g.
    !!
    !!     f_loc(p%elm_x_loc(),p%elm_y_loc(),p%elm_z_loc(),:,:,:)
    !!
    !! where p is an object of class point_cls. This returns an array of shape
    !! [0:lpd,0:lpd,0:lpd].
    !!
    !! @param self Passed-object dummy of class point_cls
    !! @return The y/phi element index
    ELEMENTAL INTEGER FUNCTION get_ey(self)
        ! Passed dummy arguments
        CLASS(point_cls), INTENT(IN) :: self

        get_ey = self%ey

    END FUNCTION


    !> Returns the z/radial element index of that point
    !!
    !! To access the element a point is located in all three element indices
    !! \f$e_\vartheta\f$, \f$e_\varphi\f$ and \f$e_r\f$ are needed. To actually
    !! select values from one of Ses3d-NT's rank 6 fields say e.g.
    !!
    !!     f_loc(p%elm_x_loc(),p%elm_y_loc(),p%elm_z_loc(),:,:,:)
    !!
    !! where p is an object of class point_cls. This returns an array of shape
    !! [0:lpd,0:lpd,0:lpd].
    !!
    !! @param self Passed-object dummy of class point_cls
    !! @return The z/radial element index
    ELEMENTAL INTEGER FUNCTION get_ez(self)
        ! Passed dummy arguments
        CLASS(point_cls), INTENT(IN) :: self

        get_ez = self%ez

    END FUNCTION


    !> Compares two objects p1 and p2 of class point_cls. In case these objects
    !! are equal `.TRUE.` is returned otherwise `.FALSE.`.
    !!
    !! @note Two objects p1 and p2 are called identical or equal if all their
    !! components are equal:
    !!
    !!  * ex: x/theta element index
    !!  * ey: y/phi element index
    !!  * ez: z/radial element index
    !!  * vm: Vandermonde Maxrix
    !!
    !! @param p1 Passed-object dummy of class point_cls
    !! @param p2 Passed-object dummy of class point_cls
    !! @return .TRUE. if p1 and p2 are identical otherwise .FALSE.
    ELEMENTAL LOGICAL FUNCTION equal_points(p1, p2)
        ! Passed dummy arguments
        CLASS(point_cls), INTENT(IN) :: p1
        CLASS(point_cls), INTENT(IN) :: p2

        equal_points = .FALSE.

        ! Because p2 may hold any extension of point_cls
        IF ( .NOT. SAME_TYPE_AS(p1, p2) ) RETURN
        ! Check components
        IF ( p1%ex/=p2%ex ) RETURN
        IF ( p1%ey/=p2%ey ) RETURN
        IF ( p1%ez/=p2%ez ) RETURN
        ! See if both vms are allocated
        IF ( ALLOCATED(p1%vm) .NEQV. ALLOCATED(p2%vm) ) RETURN
        ! Check vm's shape
        IF ( .NOT. ALL(SHAPE(p1%vm) == SHAPE(p2%vm)) ) RETURN
        ! Finally see if the content is the same
        IF ( .NOT. ALL(p1%vm==p2%vm)) RETURN

        equal_points = .TRUE.

    END FUNCTION


    !> Checks if an object of point_cls class is initialized or not
    !!
    !! @note A point - an object of class point_cls - is called uninitialized
    !! if it  was instantiated without passing any arguments (i.e. `... =
    !! point_cls()`). For details see: point_mod::uninitialized_point_cls().
    !!
    !! @param self Passed-object dummy of class point_cls
    !! @return .TRUE. if the passed object is uninitialized otherwise .FALSE.
    ELEMENTAL LOGICAL FUNCTION uninitialized(self)
        ! Passed dummy arguments
        CLASS(point_cls), INTENT(IN) :: self

        IF ( equal_points(self, point_cls()) ) THEN
            uninitialized = .TRUE.
        ELSE
            uninitialized = .FALSE.
        END IF

    END FUNCTION


    !> Converts and returns the x/theta coordinate value into the corresponding
    !! latitude
    !!
    !! This function is implemented using:
    !!
    !!  * coordinate_utilities_mod::colat2lat()
    !!  * coordinate_utilities_mod::rad2deg()
    !!  * point_mod::interpolate_x()
    !!
    !! @note The theta coordinate value is interpolated in before hand to
    !!       get precisely that very latitude - affected by FP rounding
    !!       errors - which corresponds to the one which is uses during the
    !!       simulation.
    !!
    !! @param self Passed-object dummy of class point_cls
    !! @return Latitude in [deg]
    ELEMENTAL REAL(real_kind) FUNCTION get_lat(self)
        ! Use associated entities
        USE coordinate_utilities_mod, ONLY : colat2lat, rad2deg
        ! Passed dummy arguments
        CLASS(point_cls), INTENT(IN) :: self

        get_lat = colat2lat(rad2deg( interpolate_x(self) ))

    END FUNCTION


    !> Converts and returns the y/phi coordinate value into the corresponding
    !! longitude
    !!
    !! This function is implemented using:
    !!
    !!  * coordinate_utilities_mod::rad2deg()
    !!  * point_mod::interpolate_y()
    !!
    !! @note The phi coordinate value is interpolated in before hand to
    !!       get precisely that very longitude - affected by FP rounding
    !!       errors - which corresponds to the one which is uses during the
    !!       simulation.
    !!
    !! @param self Passed-object dummy of class point_cls
    !! @return Longitude in [deg]
    ELEMENTAL REAL(real_kind) FUNCTION get_lon(self)
        ! Use associated entities
        USE coordinate_utilities_mod, ONLY : rad2deg
        ! Passed dummy arguments
        CLASS(point_cls), INTENT(IN) :: self

        get_lon = rad2deg( interpolate_y(self) )

    END FUNCTION


    !> Converts and returns the z/radial coordinate value into a depth value
    !! counting positive from the Earth's surface
    !!
    !! This function is implemented using:
    !!
    !!  * coordinate_utilities_mod::radius2depth()
    !!  * point_mod::interpolate_z()
    !!
    !! @note The radial coordinate value is interpolated in before hand to
    !!       get precisely that very depth value - affected by FP rounding
    !!       errors - which corresponds to the one which is uses during the
    !!       simulation.
    !!
    !! @param self Passed-object dummy of class point_cls
    !! @return Depth value (counting from parameters_mod::earth_radius)
    ELEMENTAL REAL(real_kind) FUNCTION get_depth(self)
        ! Use associated entities
        USE coordinate_utilities_mod, ONLY : radius2depth
        ! Passed dummy arguments
        CLASS(point_cls), INTENT(IN) :: self

        get_depth = radius2depth(interpolate_z(self))

    END FUNCTION


    !> Returns the Vandermonde-matrix
    !!
    !! The Vandermonde matrix may be used for interpolation purposes as:
    !! \f[f(\vartheta,\varphi,r) = \sum_{ijk} v_{ijk}*f_{ijk} \f]
    !! With \f$v\f$ the Vandermonde matrix and \f$f\f$ an array of shape
    !! [0:lpd,0:lpd,0:lpd] (lpd: Lagrange polynomial degree)
    !!
    !! @param self Passed-object dummy of class point_cls
    !! @return Vandermode-matrix; An array of shape [0:lpd,0:lpd,0:lpd]
    PURE FUNCTION get_vanderm(self) RESULT(vm)
        ! Passed dummy arguments
        CLASS(point_cls), INTENT(IN) :: self
        ! Result variable
        REAL(real_kind) :: vm(0:conf%lpd(),0:conf%lpd(),0:conf%lpd())

        vm(:,:,:) = self%vm

    END FUNCTION


!------------------------------------------------------------------------------
! Private module procedures


    !> Mapping from the passed local x/theta coordinate value onto the
    !! so called unit interval [-1,1].
    !!
    !! In order to apply the same numerical quadrature scheme to all elements,
    !! coordinates \f$(\vartheta, \varphi, r)\f$ need to be mapped onto the
    !! unit cube \f$ [-1,1]\times [-1,1]\times [-1,1]\f$. For each coordinat
    !! value this splits into two steps:
    !!
    !!  1. Retrieve local element index `x_elm`. See point_mod::elm_x_loc().
    !!  2. Determine the reference position `std_x` inside the actual element.
    !!
    !! Once the local element-index \f$i_{\vartheta,loc}\f$ for the coordinate
    !! value \f$\vartheta\f$ is known, the reference position calculates:
    !! \f{eqnarray*}
    !!  std_{\vartheta}: [\vartheta_{loc,min},\vartheta_{loc,max}] & \rightarrow & [-1,1] \\
    !!  \vartheta & \mapsto &
    !!  2 \big ( (\vartheta-\vartheta_{loc,min})/ \Delta\vartheta_{elm} - i_{\vartheta,loc} \big ) - 1
    !! \f}
    !! Where the subscript \f$_{loc}\f$ refers to the local MPI rank and
    !! \f$\Delta\vartheta_{elm}\f$ is denoting the width of a x/theta element.
    !!
    !! @warning If \f$\vartheta\f$ is not a coordinate value in the local MPI
    !!          rank - which is the interval \f$[\vartheta_{loc,min},
    !!          \vartheta_{loc,max}]\f$ - returned values are meaningless.
    !!
    !! @param x Co-latitude \f$\vartheta\f$ in [rad]
    !! @returns The x-coordinate reference value in [-1,1]
    REAL(real_kind) ELEMENTAL FUNCTION std_x(x)
        ! Passed dummy arguments
        REAL(real_kind), INTENT(IN) :: x

        ! The ASSOCIATE is just to keep the actual code concise and readable
        ASSOCIATE ( x_loc_min => conf%x_loc_min(), &
                    dx_elm    => conf%dx_elm(), &
                    x_elm     => elm_x_loc(x) )

        ! XXX might it be better to multiply by 2.0 at first and than divide?
        std_x = 2.0*( (x - x_loc_min)/dx_elm - x_elm ) - 1.0

        END ASSOCIATE

    END FUNCTION std_x


    !> Mapping from the passed local y/phi coordinate value onto the
    !! so called unit interval [-1,1].
    !!
    !! In order to apply the same numerical quadrature scheme to all elements,
    !! coordinates \f$(\vartheta, \varphi, r)\f$ need to be mapped onto the
    !! unit cube \f$ [-1,1]\times [-1,1]\times [-1,1]\f$. For each coordinat
    !! value this splits into two steps:
    !!
    !!  1. Retrieve local element index `y_elm`. See point_mod::elm_y_loc().
    !!  2. Determine the reference position `std_y` inside the actual element.
    !!
    !! Once the local element-index \f$i_{\varphi,loc}\f$ for the coordinate
    !! value \f$\varphi\f$ is known, the reference position calculates:
    !! \f{eqnarray*}
    !!  std_{\varphi}: [\varphi_{loc,min},\varphi_{loc,max}] & \rightarrow & [-1,1] \\
    !!  \varphi & \mapsto &
    !!  2 \big ( (\varphi-\varphi_{loc,min})/ \Delta\varphi_{elm} - i_{\varphi,loc} \big ) - 1
    !! \f}
    !! Where the subscript \f$_{loc}\f$ refers to the local MPI rank and
    !! \f$\Delta\varphi_{elm}\f$ is denoting the width of a y/phi element.
    !!
    !! @warning If \f$\varphi\f$ is not a coordinate value in the local MPI
    !!          rank - which is the interval \f$[\varphi_{loc,min},
    !!          \varphi_{loc,max}]\f$ - returned values are meaningless.
    !!
    !! @param y Longitude \f$\varphi\f$ in [rad]
    !! @returns The y-coordinate reference value in [-1,1]
    REAL(real_kind) ELEMENTAL FUNCTION std_y(y)
        ! Passed dummy arguments
        REAL(real_kind), INTENT(IN) :: y

        ! XXX might it be better to multiply by 2.0 at first and than divide?
        ASSOCIATE ( y_loc_min => conf%y_loc_min(), &
                    dy_elm    => conf%dy_elm(), &
                    y_elm     => elm_y_loc(y) )

        ! XXX might it be better to multiply first by 2.0 and than divide?
        std_y = 2.0*( (y - y_loc_min)/dy_elm - y_elm ) - 1.0

        END ASSOCIATE

    END FUNCTION std_y


    !> Mapping from the passed local z/radial coordinate value onto the
    !! so called unit interval [-1,1].
    !!
    !! In order to apply the same numerical quadrature scheme to all elements,
    !! coordinates \f$(\vartheta,  r , r)\f$ need to be mapped onto the
    !! unit cube \f$ [-1,1]\times [-1,1]\times [-1,1]\f$. For each coordinat
    !! value this splits into two steps:
    !!
    !!  1. Retrieve local element index `z_elm`. See point_mod::elm_z_loc().
    !!  2. Determine the reference position `std_z` inside the actual element.
    !!
    !! Once the local element-index \f$i_{r,loc}\f$ for the coordinate
    !! value \f$r\f$ is known, the reference position calculates:
    !! \f{eqnarray*}
    !!  std_{r}: [ r_{loc,min}, r_{loc,max}] & \rightarrow & [-1,1] \\
    !!   r  & \mapsto &
    !!  2 \big ( (r_{loc,max}-r)/ \Delta r_{elm} - i_{r,loc} \big ) - 1
    !! \f}
    !! Where the subscript \f$_{loc}\f$ refers to the local MPI rank and
    !! \f$\Delta r _{elm}\f$ is denoting the width of a z/radial element.
    !!
    !! @warning If \f$ r \f$ is not a coordinate value in the local MPI
    !!          rank - which is the interval \f$[ r_{loc,min},
    !!          r_{loc,max}]\f$ - returned values are meaningless.
    !! @warning z-coordinates are upside down
    !!
    !! @param z Longitude \f$ r \f$ in [rad]
    !! @returns The z-coordinate reference value in [-1,1]
    REAL(real_kind) ELEMENTAL FUNCTION std_z(z)
        ! Passed dummy arguments
        REAL(real_kind), INTENT(IN) :: z

        ! The ASSOCIATE is just to keep the actual code concise and readable
        ASSOCIATE ( z_loc_max => conf%z_loc_max(), &
                    dz_elm    => conf%dz_elm(), &
                    z_elm     => elm_z_loc(z) )

        ! XXX might it be better to multiply by 2.0 at first and than divide?
        ! Why the hell are z-coordinates upside down?
        std_z = 2.0*( (z_loc_max - z)/dz_elm - z_elm ) - 1.0

        END ASSOCIATE

    END FUNCTION std_z


    !> Returns .TRUE. if the passed point (x,y,z) is located inside the
    !! relaxing boundary elements
    !!
    !! To avoid unphysical reflections the lateral and lower boundaries are
    !! enclosed by so called relaxing boundary elements which attenuate the
    !! wave-field. A point \f$(\vartheta, \varphi, r) \f$ is located inside
    !! the boundary layers if:
    !!
    !! - \f$ \vartheta_{glb,min} + w_{elm} \Delta\vartheta_{elm} \leq \vartheta \leq \vartheta_{glb,max} - w_{elm} \Delta\vartheta_{elm} \f$
    !! - \f$ \varphi_{glb,min} + w_{elm} \Delta\varphi_{elm} \leq \varphi \leq \varphi_{glb,max} - w_{elm} \Delta\varphi_{elm} \f$
    !! - \f$  r_{glb,min} + w_{elm} \Delta r_{elm} \leq  r \f$
    !!
    !! Where \f$w_{elm}\f$ denotes the boundary width (counting in number of
    !! elements) and the subscripts \f$_{glb}\f$, \f$_{elm}\f$ are
    !! abbreviations for 'global' and 'element'. Here the term 'global' refers
    !! to the entire domain, not split into partitions for MPI parallelization.
    !!
    !! @note Actually, this approach is not a 100% bullet proof. Floating point
    !!       representation may cause rounding errors.
    !! @note Someone should neither place receiver nor sources within the
    !!       boundary layer.
    !!
    !! @param x Co-latitude \f$\vartheta\f$ in [rad]
    !! @param y Longitude \f$\varphi\f$ in [rad]
    !! @param z Radius \f$r\f$ counting from the Earth's center in [m]
    !! @returns .TRUE. if the point (x,y,z) is located in PML otherwise .FALSE.
    ELEMENTAL LOGICAL FUNCTION in_boundary(x, y, z)
        ! Passed dummy arguments
        REAL(real_kind), INTENT(IN) :: x, y, z

        ! The ASSOCIATE is just to keep the actual code concise and readable
        ASSOCIATE ( taper_width => conf%taper_width, &
                    x_glb_min   => conf%x_glb_min(), &
                    x_glb_max   => conf%x_glb_max(), &
                    y_glb_min   => conf%y_glb_min(), &
                    y_glb_max   => conf%y_glb_max(), &
                    z_glb_min   => conf%z_glb_min_, &
                    dx_elm      => conf%dx_elm(), &
                    dy_elm      => conf%dy_elm(), &
                    dz_elm      => conf%dz_elm() )

        IF ( ( x_glb_min + taper_width*dx_elm ) > x .OR. &
             x > ( x_glb_max - taper_width*dx_elm ) .OR. &
             ( y_glb_min + taper_width*dy_elm ) > y .OR. &
             y > ( y_glb_max - taper_width*dy_elm ) .OR. &
             ( z_glb_min + taper_width*dz_elm ) > z ) THEN
            in_boundary = .TRUE.
        ELSE
            in_boundary = .FALSE.
        END IF

        END ASSOCIATE

    END FUNCTION


    !> Returns .TRUE. if point (x,y,z) is located inside of the computational
    !! domain otherwise .FALSE. is the result.
    !!
    !! A point \f$(\vartheta, \varphi, r) \f$ is located inside of the
    !! computational domain if:
    !! * \f$ \vartheta_{glb,min} < \vartheta < \vartheta_{glb,max}) \f$
    !! * \f$ \varphi_{glb,min}   < \varphi   < \varphi_{glb,max}) \f$
    !! * \f$ r_{glb,min}         < r      \leq r_{glb,max}) \f$ (surface)
    !!
    !! Where the subscript \f$_{glb}\f$ is an abbreviation for 'global'.
    !! Here the term 'global' refers to the entire domain, not split into
    !! partitions for MPI parallelization.
    !!
    !! @note Mind the less or equal at the surface.
    !!
    !! @param x Co-latitude \f${\vartheta}\f$ in[rad]
    !! @param y Longitude \f${\varphi}\f$ in [rad]
    !! @param z Radius \f${r}\f$ from the Earth's center in [m]
    !! @returns True if point is located in computational domain otherwise false
    ELEMENTAL LOGICAL FUNCTION in_domain(x, y, z)
        ! Passed dummy arguments
        REAL(real_kind), INTENT(IN) :: x, y, z

        ! The ASSOCIATE is just to keep the actual code concise and readable
        ASSOCIATE ( x_glb_min   => conf%x_glb_min(), &
                    x_glb_max   => conf%x_glb_max(), &
                    y_glb_min   => conf%y_glb_min(), &
                    y_glb_max   => conf%y_glb_max(), &
                    z_glb_min   => conf%z_glb_min_, &
                    z_glb_max   => conf%z_glb_max_ )

        IF ( x_glb_min < x .AND. x <  x_glb_max .AND. &
             y_glb_min < y .AND. y <  y_glb_max .AND. &
             z_glb_min < z .AND. z <= z_glb_max ) THEN
            in_domain = .TRUE.
        ELSE
            in_domain = .FALSE.
        END IF

        END ASSOCIATE

    END FUNCTION in_domain


    !> Checks if the passed coordinates are inside the local MPI rank. Returns
    !! .TRUE. if the passed point (x,y,z) is located in this very MPI rank
    !! otherwise .FALSE. is the result.
    !!
    !! A point \f$(\vartheta, \varphi, r) \f$ is located inside that very MPI
    !! rank if:
    !!
    !! - \f$ \vartheta_{loc,min} < \vartheta \leq \vartheta_{loc,max} \f$
    !! - \f$ \varphi_{loc,min} < \varphi \leq \varphi_{loc,max} \f$
    !! - \f$ r_{loc,min} < r \leq r_{loc,max} \f$
    !!
    !! Where \f$loc\f$ is an abbreviation for the 'local' MPI rank. Here the
    !! term 'local' refers to a MPI partition. Remember, for parallelization
    !! the entire domain is split into partitions.
    !!
    !! @note That approach circumvents duplicates at partition boundaries.
    !!
    !! @param x Co-latitude \f$\vartheta\f$ in [rad]
    !! @param y Longitude \f$\varphi\f$ in [rad]
    !! @param z Radius \f$r\f$ counting from the Earth's center in [m]
    !! @returns .TRUE. if point is located in local rank otherwise .FALSE.
    ELEMENTAL LOGICAL FUNCTION in_rank(x, y, z)
        ! Passed dummy arguments
        REAL(real_kind), INTENT(IN) :: x, y, z

        ! The ASSOCIATE is just to keep the actual code concise and readable
        ASSOCIATE ( x_loc_min => conf%x_loc_min(), &
                    x_loc_max => conf%x_loc_max(), &
                    y_loc_min => conf%y_loc_min(), &
                    y_loc_max => conf%y_loc_max(), &
                    z_loc_min => conf%z_loc_min(), &
                    z_loc_max => conf%z_loc_max() )

        IF ( x_loc_min < x .AND. x <= x_loc_max .AND. &
             y_loc_min < y .AND. y <= y_loc_max .AND. &
             z_loc_min < z .AND. z <= z_loc_max ) THEN
            in_rank = .TRUE.
        ELSE
            in_rank = .FALSE.
        END IF

        END ASSOCIATE

    END FUNCTION


    !> Returns the local x/theta element index corresponding to a points
    !! x/theta coordinate value
    !!
    !! The index is determined as follows:
    !!   1. Calculate the local offset: \f$ x_{loc,max} - x \f$
    !!   2. Divide it by the width of the element: \f$ \frac{x_{loc,max} - x}{\Delta x_{elm}} \f$
    !!   3. Get the greatest integer less than or equal: \f$ {\mathrm FLOOR}\left( \frac{x_{loc,max} - x}{\Delta x_{elm}}\right) \f$
    !!   4. Subtract if from the local number of elements: \f$ nx_{loc} - {\mathrm FLOOR}\left( \frac{x_{loc,max} - x}{\Delta x_{elm}}\right) \f$
    !!
    !! where \f$loc\f$ indicates the local MPI-rank, \f$nx\f$ denotes the local
    !! number of elements and \f$x\f$ stands for the radial coordinate value.
    !!
    !! @todo I do not understand why step 4 is necessary.
    !!
    !! @warning Assumes array indices counting from 0 to nx_loc_max
    !! @warning If the point is not located in the local MPI rank the returned
    !!          value is meaningless
    !!
    !! @param x The x/theta coordinate value
    !! @return Element index in x-direction with respect to the local MPI rank
    ELEMENTAL INTEGER FUNCTION elm_x_loc(x)
        ! Passed dummy arguments
        REAL(real_kind), INTENT(IN) :: x

        elm_x_loc = conf%nx_loc() - FLOOR((conf%x_loc_max() - x)/conf%dx_elm())

    END FUNCTION


    !> Returns the local y/phi element index corresponding to a points
    !! y/phi coordinate value
    !!
    !! The index is determined as follows:
    !!   1. Calculate the local offset: \f$ y_{loc,max} - y \f$
    !!   2. Divide it by the width of the element: \f$ \frac{y_{loc,max} - y}{\Delta y_{elm}} \f$
    !!   3. Get the greatest integer less than or equal: \f$ {\mathrm FLOOR}\left( \frac{y_{loc,max} - y}{\Delta y_{elm}}\right) \f$
    !!   4. Subtract if from the local number of elements: \f$ ny_{loc} - {\mathrm FLOOR}\left( \frac{y_{loc,max} - y}{\Delta y_{elm}}\right) \f$
    !!
    !! where \f$loc\f$ indicates the local MPI-rank, \f$ny\f$ denotes the local
    !! number of elements and \f$y\f$ stands for the radial coordinate value.
    !!
    !! @todo I do not understand why step 4 is necessary.
    !!
    !! @warning Assumes array indices counting from 0 to ny_loc_max
    !! @warning If the point is not located in the local MPI rank the returned
    !!          value is meaningless
    !!
    !! @param y The y/phi coordinate value
    !! @return Element index in y-direction with respect to the local MPI rank
    ELEMENTAL INTEGER FUNCTION elm_y_loc(y)
        ! Passed dummy arguments
        REAL(real_kind), INTENT(IN) :: y

        elm_y_loc = conf%ny_loc() - FLOOR((conf%y_loc_max() - y)/conf%dy_elm())

    END FUNCTION


    !> Returns the local z/radial-element index corresponding to a points
    !! r/radial coordinate value
    !!
    !! The index is determined as follows:
    !!   1. Calculate the local offset: \f$ z_{loc,max} - z \f$
    !!   2. Divide it by the width of the element: \f$ \frac{z_{loc,max} - z}{\Delta z_{elm}} \f$
    !!   3. Get the greatest integer less than or equal: \f$ {\mathrm FLOOR}\left( \frac{z_{loc,max} - z}{\Delta z_{elm}}\right) \f$
    !!
    !! where \f$ loc \f$ indicates the local MPI-rank and \f$z\f$ stands for
    !! the radial coordinate value.
    !!
    !! @warning Assumes array indices counting from 0 to nz_loc_max
    !! @warning z-coordinates are upside down
    !! @warning If the point is not located in the local MPI rank the returned
    !!          value is meaningless
    !!
    !! @param z The r/radial coordinate value
    !! @return Element index in z-direction with respect to the local MPI rank
    ELEMENTAL INTEGER FUNCTION elm_z_loc(z)
        ! Passed dummy arguments
        REAL(real_kind), INTENT(IN) :: z

        ! Why the hell are z-coordinates upside down?
        elm_z_loc = FLOOR( (conf%z_loc_max() - z)/conf%dz_elm() )

    END FUNCTION


    !> Returns the interpolated x/theta coordinate of a point
    !!
    !! The interpolation scheme is:
    !! \f[\vartheta = \sum_{ijk} v_{ijk}*\vartheta_{e_\vartheta e_\varphi e_r,ijk} \f]
    !! With \f$v\f$ the Vandermonde matrix and \f$e\f$ the point's element
    !! indices.
    !!
    !! @note This is useful to get precisely that very value - affected by FP
    !!       rounding errors - Ses3d-NT actually uses during the simulation
    !!
    !! @param self Passed-object dummy of class point_cls
    !! @return x/theta coordinate value of the passed point
    ELEMENTAL REAL(real_kind) FUNCTION interpolate_x(self)
        ! Use associated entities
        USE geometric_paras_mod, ONLY : theta
        ! Passed dummy arguments
        CLASS(point_cls), INTENT(IN) :: self

        interpolate_x = SUM( self%vm*theta(self%ex,self%ey,self%ez,:,:,:) )

    END FUNCTION


    !> Returns the interpolated y/phi coordinate of a point
    !!
    !! The interpolation scheme is:
    !! \f[\varphi = \sum_{ijk} v_{ijk}*\varphi_{e_\vartheta e_\varphi e_r,ijk} \f]
    !! With \f$v\f$ the Vandermonde matrix and \f$e\f$ the point's element
    !! indices.
    !!
    !! @note This is useful to get precisely that very value - affected by FP
    !!       rounding errors - Ses3d-NT actually uses during the simulation
    !!
    !! @param self Passed-object dummy of class point_cls
    !! @return y/phi coordinate value of the passed point
    ELEMENTAL REAL(real_kind) FUNCTION interpolate_y(self)
        ! Use associated entities
        USE geometric_paras_mod, ONLY : phi
        ! Passed dummy arguments
        CLASS(point_cls), INTENT(IN) :: self

        interpolate_y = SUM( self%vm*phi(self%ex,self%ey,self%ez,:,:,:) )

    END FUNCTION


    !> Returns the interpolated z/radial coordinate of a point
    !!
    !! The interpolation scheme is:
    !! \f[r = \sum_{ijk} v_{ijk}*r_{e_\vartheta e_\varphi e_r,ijk} \f]
    !! With \f$v\f$ the Vandermonde matrix and \f$e\f$ the point's element
    !! indices.
    !!
    !! @note This is useful to get precisely that very value - affected by FP
    !!       rounding errors - Ses3d-NT actually uses during the simulation
    !!
    !! @param self Passed-object dummy of class point_cls
    !! @return z/radial coordinate value of the passed point
    ELEMENTAL REAL(real_kind) FUNCTION interpolate_z(self)
        ! Use associated entities
        USE geometric_paras_mod, ONLY : r
        ! Passed dummy arguments
        CLASS(point_cls), INTENT(IN) :: self

        interpolate_z = SUM( self%vm*r(self%ex,self%ey,self%ez,:,:,:) )

    END FUNCTION


END MODULE point_mod
