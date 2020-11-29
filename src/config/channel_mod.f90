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
!! $Date: 2013-12-02 09:25:30 +0100 (Mon, 02 Dec 2013) $
!! $Author: mauerberger $
!! $Revision: 832 $
!! @copyright GNU General Public License version 3 or later



!> This module provides a fundamental Ses3d-NT class which is used for
!! recording and adding observables. It covers recorders, sources and volume
!! snapshots.
MODULE channel_mod
    USE parameters_mod, ONLY : real_kind
    USE configuration_mod, ONLY : config => configuration
    USE grid_mod, ONLY : lpd ! Workaround for a GCC bug

    IMPLICIT NONE

    PRIVATE

    !> A class which handles different channels in Ses3d-NT to:
    !!  * Record synthetics
    !!  * Add sources
    !!  * Return volume snapshots
    !! To give an examples this might be the radial component of the
    !! displacement velocity.
    TYPE, ABSTRACT :: channel_cls
    PRIVATE
        !> Either records a channel or stores the source signature
        REAL(real_kind), ALLOCATABLE :: data(:)
        !> The channels trivial name (e.g. vx, rot_vx ...)
        CHARACTER(LEN=8) :: attribute
    CONTAINS
        !> Dummy procedure which simply does nothing. Actual implementations
        !! are overriding it by procedures recording seismograms.
        PROCEDURE :: record => record_nothing
        !> Dummy procedure which simply does nothing. Actual implementations
        !! are overriding it with procedures adding an exciting source
        PROCEDURE :: add => add_nothing
        !> Dummy procedure which simply does nothing. Actual implementations
        !! are overriding it with functions returning the corresponding field.
        PROCEDURE, NOPASS :: field => field_nothing
    END TYPE


    !> A container class ... not the very best approach
    !!
    !!@todo split into a recorder container, a source container and a field
    !!      container
    TYPE :: channel_container
    PRIVATE
        !> The actual channel object
        CLASS(channel_cls), ALLOCATABLE :: channel
    CONTAINS
        !> Returns the data array of the component channel. Depending on the
        !! actual instantiation of channel this may either be a synthetic
        !! seismogram or the source wavelet
        PROCEDURE :: data
        !> Returns the channels trivial name (e.g. vx, rot_vx ...)
        PROCEDURE :: attribute
        !> Procedure to record synthetic seismograms
        PROCEDURE :: record
        !> Returns the whole field corresponding to the channel
        PROCEDURE :: field
        !> Adds the source
        PROCEDURE :: add
        !> Returns true if the channel is unknown
        PROCEDURE :: unknown
        !> For comparing channel_container objects. Do not use this bound
        !! procedure explicitly. Use the comparison operator '==' instead.
        !! For details see channel_mod::equal_channels().
        PROCEDURE, PRIVATE :: equal_channels
        !> Defines the comparison operator '==' for channel_container classes.
        !! For details see channel_mod::equal_channels().
        GENERIC :: OPERATOR(==) => equal_channels
    END TYPE


    !> An object of class channel_unknown will be instantiated if something goes
    !! wrong e.g. the channel attribute is unknown.
    TYPE, EXTENDS(channel_cls) :: channel_unknown
    CONTAINS
    END TYPE

    !> To record the theta-component of the displacement velocity
    TYPE, EXTENDS(channel_cls) :: channel_vx
    CONTAINS
        PROCEDURE :: record => record_vx
        PROCEDURE, NOPASS :: field => field_vx
    END TYPE

    !> To record the north-component of the displacement velocity
    TYPE, EXTENDS(channel_cls) :: channel_n
    CONTAINS
        PROCEDURE :: record => record_n
    END TYPE

    !> To record the phi-component of the displacement velocity
    TYPE, EXTENDS(channel_cls) :: channel_vy
    CONTAINS
        PROCEDURE :: record => record_vy
        PROCEDURE, NOPASS :: field => field_vy
    END TYPE

    !> To record the radial-component of the displacement velocity
    TYPE, EXTENDS(channel_cls) :: channel_vz
    CONTAINS
        PROCEDURE :: record => record_vz
        PROCEDURE, NOPASS :: field => field_vz
    END TYPE

    !> To record the curls theta-component of displacement velocity
    TYPE, EXTENDS(channel_cls) :: channel_rot_vx
    CONTAINS
        PROCEDURE :: record => record_rot_vx
        PROCEDURE, NOPASS :: field => field_rot_vx
    END TYPE

    !> To record the curls phi-component of displacement velocity
    TYPE, EXTENDS(channel_cls) :: channel_rot_vy
    CONTAINS
        PROCEDURE :: record => record_rot_vy
        PROCEDURE, NOPASS :: field => field_rot_vy
    END TYPE

    !> To record the curls radial-component of displacement velocity
    TYPE, EXTENDS(channel_cls) :: channel_rot_vz
    CONTAINS
        PROCEDURE :: record => record_rot_vz
        PROCEDURE, NOPASS :: field => field_rot_vz
    END TYPE

    !> To record hydrostatic pressure rate
    TYPE, EXTENDS(channel_cls) :: channel_pressure
    CONTAINS
        PROCEDURE :: record => record_pressure
        PROCEDURE, NOPASS :: field => field_pressure
    END TYPE


    !> To add or record theta-theta-direction dipole force sources
    TYPE, EXTENDS(channel_cls) :: channel_mxx
    CONTAINS
        PROCEDURE :: add => add_mxx
        PROCEDURE :: record => record_mxx
    END TYPE

    !> To add or record theta-phi-direction dipole force sources
    TYPE, EXTENDS(channel_cls) :: channel_mxy
    CONTAINS
        PROCEDURE :: add => add_mxy
        PROCEDURE :: record => record_mxy
    END TYPE

    !> To add or record theta-r-direction dipole force sources
    TYPE, EXTENDS(channel_cls) :: channel_mxz
    CONTAINS
        PROCEDURE :: add => add_mxz
        PROCEDURE :: record => record_mxz
    END TYPE

    !> To add or record phi-phi-direction dipole force sources
    TYPE, EXTENDS(channel_cls) :: channel_myy
    CONTAINS
        PROCEDURE :: add => add_myy
        PROCEDURE :: record => record_myy
    END TYPE

    !> To add or record phi-r-direction dipole force sources
    TYPE, EXTENDS(channel_cls) :: channel_myz
    CONTAINS
        PROCEDURE :: add => add_myz
        PROCEDURE :: record => record_myz
    END TYPE

    !> To add or record rr-direction dipole force sources
    TYPE, EXTENDS(channel_cls) :: channel_mzz
    CONTAINS
        PROCEDURE :: add => add_mzz
        PROCEDURE :: record => record_mzz
    END TYPE

    !> To add or record theta-component monopole force sources
    TYPE, EXTENDS(channel_cls) :: channel_fx
    CONTAINS
        PROCEDURE :: record => record_fx
        PROCEDURE :: add => add_fx
    END TYPE

    !> To add or record phi-component monopole force sources
    TYPE, EXTENDS(channel_cls) :: channel_fy
    CONTAINS
        PROCEDURE :: record => record_fy
        PROCEDURE :: add => add_fy
    END TYPE

    !> To add or record radial-component monopole forces
    TYPE, EXTENDS(channel_cls) :: channel_fz
    CONTAINS
        PROCEDURE :: record => record_fz
        PROCEDURE :: add => add_fz
    END TYPE

    !> To record theta-theta component strain
    TYPE, EXTENDS(channel_cls) :: channel_exx
    CONTAINS
        PROCEDURE, NOPASS :: field => field_exx
    END TYPE

    !> To record theta-phi component strain
    TYPE, EXTENDS(channel_cls) :: channel_exy
    CONTAINS
        PROCEDURE, NOPASS :: field => field_exy
    END TYPE

    !> To record theta-r component strain
    TYPE, EXTENDS(channel_cls) :: channel_exz
    CONTAINS
        PROCEDURE, NOPASS :: field => field_exz
    END TYPE

    !> To record phi-phi component strain
    TYPE, EXTENDS(channel_cls) :: channel_eyy
    CONTAINS
        PROCEDURE, NOPASS :: field => field_eyy
    END TYPE

    !> To record phi-r component strain
    TYPE, EXTENDS(channel_cls) :: channel_eyz
    CONTAINS
        PROCEDURE, NOPASS :: field => field_eyz
    END TYPE

    !> To record r-r component strain
    TYPE, EXTENDS(channel_cls) :: channel_ezz
    CONTAINS
        PROCEDURE, NOPASS :: field => field_ezz
    END TYPE

    !> Structure constructor for the channel_container class. Either used to
    !!  * instantiate an object to record data (receiver)
    !!  * or to setup an exiting source (monopole or dipole)
    INTERFACE channel_container
        MODULE PROCEDURE construct_channel_container_receiver
        MODULE PROCEDURE construct_channel_container_source
    END INTERFACE

    ! Just expose the class. All the rest remains under the hoods.
    PUBLIC :: channel_container

CONTAINS

    ELEMENTAL FUNCTION construct_channel_container_receiver(attribute) &
            RESULT(container)
        CHARACTER(LEN=*), INTENT(IN) :: attribute
        TYPE(channel_container) :: container
        CLASS(channel_cls), ALLOCATABLE :: channel_
        REAL(real_kind) :: zeros_(config%nt())

        zeros_ = 0.0

        SELECT CASE(attribute)
            ! Record displacement velocity
            CASE('vx')
                ALLOCATE( channel_, &
                    source=channel_vx(attribute=attribute, data=zeros_) )
            CASE ('n')
                ALLOCATE( channel_, &
                    source=channel_n(attribute=attribute, data=zeros_) )
            CASE('vy', 'e')
                ALLOCATE( channel_, &
                    source=channel_vy(attribute=attribute, data=zeros_) )
            CASE('vz', 'z')
                ALLOCATE( channel_, &
                    source=channel_vz(attribute=attribute, data=zeros_) )
            ! Record curl of displacement velocity
            CASE ('rot_vx')
                ALLOCATE( channel_, &
                    source=channel_rot_vx(attribute=attribute, data=zeros_) )
            CASE ('rot_vy')
                ALLOCATE( channel_, &
                    source=channel_rot_vy(attribute=attribute, data=zeros_) )
            CASE ('rot_vz')
                ALLOCATE( channel_, &
                    source=channel_rot_vz(attribute=attribute, data=zeros_) )
            ! Record moment-tensor forces
            CASE ('mxx')
                ALLOCATE( channel_, &
                    source=channel_mxx(attribute=attribute, data=zeros_) )
            CASE ('mxy','myx')
                ALLOCATE( channel_, &
                    source=channel_mxy(attribute=attribute, data=zeros_) )
            CASE ('mxz','mzx')
                ALLOCATE( channel_, &
                    source=channel_mxz(attribute=attribute, data=zeros_) )
            CASE ('myy')
                ALLOCATE( channel_, &
                    source=channel_myy(attribute=attribute, data=zeros_) )
            CASE ('myz','mzy')
                ALLOCATE( channel_, &
                    source=channel_myz(attribute=attribute, data=zeros_) )
            CASE ('mzz')
                ALLOCATE( channel_, &
                    source=channel_mzz(attribute=attribute, data=zeros_) )
            ! Record monopole forces
            CASE ('fx')
                ALLOCATE( channel_, &
                    source=channel_fx(attribute=attribute, data=zeros_) )
            CASE ('fy')
                ALLOCATE( channel_, &
                    source=channel_fy(attribute=attribute, data=zeros_) )
            CASE ('fz')
                ALLOCATE( channel_, &
                    source=channel_fz(attribute=attribute, data=zeros_) )
            ! Allocate to 'nothing' if attribute is not known
            CASE DEFAULT
                ALLOCATE( channel_, &
                    source=channel_unknown(attribute=attribute, data=zeros_) )
        END SELECT

        ! XXX isn't it possible to use the default constructor?
        !container = channel_container(channel=channel_)
        ALLOCATE(container%channel, SOURCE=channel_)

    END FUNCTION construct_channel_container_receiver

    PURE FUNCTION construct_channel_container_source(attribute, data) &
            RESULT(container)
        CHARACTER(LEN=*), INTENT(IN) :: attribute
        REAL(real_kind), INTENT(IN) :: data(:)
        TYPE(channel_container) :: container
        CLASS(channel_cls), ALLOCATABLE :: channel_

        SELECT CASE(attribute)
            ! Monopole forces
            CASE('fx')
                ALLOCATE(channel_, &
                    source=channel_fx(attribute=attribute, data=data))
            CASE('fy')
                ALLOCATE(channel_, &
                    source=channel_fy(attribute=attribute, data=data))
            CASE('fz')
                ALLOCATE(channel_, &
                    source=channel_fz(attribute=attribute, data=data))
            ! Moment-tensor forces
            CASE ('mxx')
                ALLOCATE(channel_, &
                    source=channel_mxx(attribute=attribute, data=data))
            CASE ('mxy','myx')
                ALLOCATE( channel_, &
                    source=channel_mxy(attribute=attribute, data=data) )
            CASE ('mxz','mzx')
                ALLOCATE( channel_, &
                    source=channel_mxz(attribute=attribute, data=data) )
            CASE ('myy')
                ALLOCATE( channel_, &
                    source=channel_myy(attribute=attribute, data=data) )
            CASE ('myz','mzy')
                ALLOCATE( channel_, &
                    source=channel_myz(attribute=attribute, data=data) )
            CASE ('mzz')
                ALLOCATE( channel_, &
                    source=channel_mzz(attribute=attribute, data=data) )
            ! Allocate to 'nothing' if attribute is not known
            CASE DEFAULT
                ALLOCATE( channel_, &
                    source=channel_unknown(attribute=attribute, data=data) )
        END SELECT

        ! XXX isn't it possible to use the default constructor?
        !container = channel_container(channel=channel_)
        ALLOCATE(container%channel, SOURCE=channel_)

    END FUNCTION construct_channel_container_source



    !======================!
    ! Container class TBPs !
    !======================!

    !> Returns .TRUE. if channel name is unknown
    ELEMENTAL FUNCTION unknown(self)
        CLASS(channel_container), INTENT(IN) :: self
        LOGICAL :: unknown

        SELECT TYPE (channel => self%channel)
            TYPE IS (channel_unknown)
                unknown = .TRUE.
            CLASS DEFAULT
                unknown = .FALSE.
        END SELECT

    END FUNCTION

    !> Returns the channel name
    ELEMENTAL FUNCTION attribute(self)
        CLASS(channel_container), INTENT(IN) :: self
        CHARACTER(LEN=8) :: attribute
        attribute = TRIM(self%channel%attribute)
    END FUNCTION

    !> Returns the channel's data section.
    PURE FUNCTION data(self)
        CLASS(channel_container), INTENT(IN) :: self
        REAL(real_kind) :: data(SIZE(self%channel%data))
        data = self%channel%data
    END FUNCTION

    !> Records the current value of the desired channel into the data section
    !! at timestep/position 'it'
    !!
    !! @param self An object of class channel_container
    !! @param it Current timestep; Position in data array
    !! @param point A point_mod::point_cls object
    ELEMENTAL SUBROUTINE record(self, it, point)
        USE point_mod, ONLY : point_cls
        CLASS(channel_container), INTENT(INOUT) :: self
        INTEGER, INTENT(IN) :: it
        TYPE(point_cls), INTENT(IN) :: point
        CALL self%channel%record( it=it, &
                                  vanderm=point%vanderm(), &
                                  ex=point%elm_x_loc(), &
                                  ey=point%elm_y_loc(), &
                                  ez=point%elm_z_loc() )
    END SUBROUTINE

    !> Adds to field
    SUBROUTINE add(self, it, point)
        USE point_mod, ONLY : point_cls
        CLASS(channel_container), INTENT(IN) :: self
        INTEGER, INTENT(IN) :: it
        TYPE(point_cls), INTENT(IN) :: point
        CALL self%channel%add( it=it, &
                               vanderm=point%vanderm(), &
                               ex=point%elm_x_loc(), &
                               ey=point%elm_y_loc(), &
                               ez=point%elm_z_loc() )
    END SUBROUTINE

    !> Returns current field
    PURE FUNCTION field(self)
        CLASS(channel_container), INTENT(IN) :: self
        REAL(real_kind) :: field(1+config%nx_loc(), &
                                 1+config%ny_loc(), &
                                 1+config%ny_loc(), &
                                 1+config%lpd(), &
                                 1+config%lpd(), &
                                 1+config%lpd() )
        field(:,:,:,:,:,:) = self%channel%field()
    END FUNCTION


    !> Compares two objects c1 and c2 of class channel_container. In case
    !! these objects are equal `.TRUE.` is returned otherwise `.FALSE.`.
    !!
    !! @note Two channels c1 and c2 are called identical or equal if all their
    !!       components are equal:
    !!
    !!  * data: The data section, array of reals
    !!  * attribute: The observable, 8-chars
    !!
    !! @param c1 Passed-object dummy of class channel_container
    !! @param c2 Passed-object dummy of class channel_container
    !! @return .TRUE. if c1 and c2 are identical otherwise .FALSE.
    ELEMENTAL LOGICAL FUNCTION equal_channels(c1, c2)
        ! Passed dummy arguments
        CLASS(channel_container), INTENT(IN) :: c1, c2

        IF ( attribute(c1)==attribute(c2) .AND. &
             ALL(data(c1)==data(c2)) ) THEN
            equal_channels = .TRUE.
        ELSE
            equal_channels = .FALSE.
        END IF

    END FUNCTION


    !====================!
    ! Channel class TBPs !
    !====================!

!------------------------------------------------------------------------------
! vx

    !> Records the displacement velocities' theta component at timestep it
    PURE SUBROUTINE record_vx(self, it, vanderm, ex, ey, ez)
        USE elastic_vars_mod, ONLY : vx
        CLASS(channel_vx), INTENT(INOUT) :: self
        INTEGER, INTENT(IN) :: it, ex, ey, ez
        REAL(real_kind), INTENT(IN) :: vanderm(lpd(config)+1,lpd(config)+1,lpd(config)+1)

        self%data(it) = SUM( vanderm(:,:,:) * vx(ex,ey,ez,:,:,:) )

    END SUBROUTINE record_vx

    !> Returns the current theta component of the displacement field
    !!
    !! @return theta displacement field, array of rank 6
    PURE FUNCTION field_vx() RESULT(field)
        USE elastic_vars_mod, ONLY : vx
        REAL(real_kind) :: field(0:config%nx_loc(),0:config%ny_loc(),&
                                 0:config%nz_loc(),0:config%lpd(),&
                                 0:config%lpd(),0:config%lpd())

        field(:,:,:,:,:,:) = vx

    END FUNCTION field_vx


!------------------------------------------------------------------------------
! n

    PURE SUBROUTINE record_n(self, it, vanderm, ex, ey, ez)
        USE elastic_vars_mod, ONLY : vx
        CLASS(channel_n), INTENT(INOUT) :: self
        INTEGER, INTENT(IN) :: it, ex, ey, ez
        REAL(real_kind), INTENT(IN) :: vanderm(lpd(config)+1,lpd(config)+1,lpd(config)+1)

        self%data(it) = SUM( -1.0*vanderm(:,:,:)*vx(ex,ey,ez,:,:,:) )

    END SUBROUTINE record_n


!------------------------------------------------------------------------------
! vy

    !> Records the displacement velocities' phi component at timestep it
    PURE SUBROUTINE record_vy(self, it, vanderm, ex, ey, ez)
        USE elastic_vars_mod, ONLY : vy
        CLASS(channel_vy), INTENT(INOUT) :: self
        INTEGER, INTENT(IN) :: it, ex, ey, ez
        REAL(real_kind), INTENT(IN) :: vanderm(lpd(config)+1,lpd(config)+1,lpd(config)+1)

        self%data(it) = SUM( vanderm(:,:,:) * vy(ex,ey,ez,:,:,:) )

    END SUBROUTINE record_vy

    !> Returns the current phi component of the displacement field
    PURE FUNCTION field_vy() RESULT(field)
        USE elastic_vars_mod, ONLY : vy
        REAL(real_kind) :: field(0:config%nx_loc(),0:config%ny_loc(),&
                                 0:config%nz_loc(),0:config%lpd(),&
                                 0:config%lpd(),0:config%lpd())

        field(:,:,:,:,:,:) = vy

    END FUNCTION field_vy


!------------------------------------------------------------------------------
! vz

    !> Records the displacement velocities' r component at timestep it
    PURE SUBROUTINE record_vz(self, it, vanderm, ex, ey, ez)
        USE elastic_vars_mod, ONLY : vz
        CLASS(channel_vz), INTENT(INOUT) :: self
        INTEGER, INTENT(IN) :: it, ex, ey, ez
        REAL(real_kind), INTENT(IN) :: vanderm(lpd(config)+1,lpd(config)+1,lpd(config)+1)

        self%data(it) = SUM( vanderm(:,:,:) * vz(ex,ey,ez,:,:,:) )

    END SUBROUTINE record_vz

    !> Returns the current r component of the displacement field
    PURE FUNCTION field_vz() RESULT(field)
        USE elastic_vars_mod, ONLY : vz
        REAL(real_kind) :: field(0:config%nx_loc(),0:config%ny_loc(),&
                                 0:config%nz_loc(),0:config%lpd(),&
                                 0:config%lpd(),0:config%lpd())

        field(:,:,:,:,:,:) = vz

    END FUNCTION field_vz


!------------------------------------------------------------------------------
! rot_vx

    !> Records the curl of the displacement-velocity in theta direction
    !!
    !!  In spherical coordinates the curl's theta (i.e.) x component of a vector-field v reads as
    !! \f[ \mathrm{rot}(\mathbf v)_\theta =
    !! \frac{1}{r} \left( \frac{1}{\sin\theta} \frac{\partial \mathrm v_r}{\partial \phi} - \frac{\partial}{\partial r} \left( r \mathrm v_\phi \right) \right) =
    !! \left( \frac{1}{r\sin\theta} \frac{\partial \mathrm v_r}{\partial \phi} - \frac{\mathrm v_\phi}{r} - \frac{\partial \mathrm v_\phi}{\partial r} \right)
    !! \f]
    !! with \f$ \dot\epsilon_{r\phi} = \frac{\partial \mathrm v_r}{\partial \phi} \f$
    !! and \f$ \dot\epsilon_{\phi r} = \frac{\partial \mathrm v_\phi}{\partial r} \f$ it follows
    !! \f[ \mathrm{rot}(\mathbf v)_\theta =
    !! \frac{\dot\epsilon_{r \phi}}{r\sin\theta} - \frac{\mathrm v_\phi}{r} - \dot\epsilon_{\phi r}
    !! \f]
    !! @todo Check if strain indices i.e. e?? are positioned properly
    !! @todo Check eyz vs. ezy aren't they the same?
    !! @todo Where does the minus comes from?
    !!
    !! @param self A channel_rot_vx object; usually passed as the dummy object
    !! @param it Current timestep; Position in data array
    !! @param vanderm The Vandermonde matrix used for interpolation
    !! @param ex Element index in theta direction it
    !! @param ex Element index in phi direction it
    !! @param ex Element index in r direction it
    PURE SUBROUTINE record_rot_vx(self, it, vanderm, ex, ey, ez)
        USE elastic_vars_mod, ONLY : vy, ezy, eyz
        USE geometric_paras_mod, ONLY : r, r_sin_theta
        ! Passed dummy arguments
        CLASS(channel_rot_vx), INTENT(INOUT) :: self
        INTEGER, INTENT(IN) :: it, ex, ey, ez
        REAL(real_kind), INTENT(IN) :: vanderm(lpd(config)+1,lpd(config)+1,lpd(config)+1)
        ! Local variables
        REAL(real_kind) :: elm(config%lpd()+1, &
                               config%lpd()+1, &
                               config%lpd()+1)

        elm(:,:,:) = -vy(ex,ey,ez,:,:,:) / r(ex,ey,ez,:,:,:) + &
            eyz(ex,ey,ez,:,:,:) / r_sin_theta(ex,ey,ez,:,:,:) + &
            ezy(ex,ey,ez,:,:,:)

        self%data(it) = SUM( vanderm(:,:,:) * elm(:,:,:) )

    END SUBROUTINE record_rot_vx

    !> Returns the theta component of the curl of the displacement field
    !!
    !! @return field ; array of rank 6
    PURE FUNCTION field_rot_vx() RESULT(field)
        USE elastic_vars_mod, ONLY : vy, eyz, ezy
        USE geometric_paras_mod, ONLY : r, r_sin_theta
        REAL(real_kind) :: field(0:config%nx_loc(),0:config%ny_loc(),&
                                 0:config%nz_loc(),0:config%lpd(),&
                                 0:config%lpd(),0:config%lpd())

        field(:,:,:,:,:,:) = -vy/r + eyz/r_sin_theta + ezy

    END FUNCTION


!------------------------------------------------------------------------------
! rot_vy

    !> @brief Records the curl of the displacement-velocity in phi direction
    !> @details In spherical coordinates the curl's phi (i.e. y) component of a vector-field v reads as
    !! \f[ \mathrm{rot}(\mathbf v)_\phi =
    !! \frac{1}{r}\left( \frac{\partial}{\partial r} \left( r \mathrm v_\theta \right) - \frac{\partial \mathrm v_r}{\partial \theta} \right) =
    !! \frac{\mathrm v_\theta}{r} + \frac{\partial \mathrm v_\theta}{\partial r} - \frac{1}{r} \frac{\partial \mathrm v_r}{\partial \theta}
    !! \f]
    !! with \f$ \dot\epsilon_{r \theta} = \frac{\partial \mathrm v_r}{\partial \theta} \f$
    !! and \f$ \dot\epsilon_{\theta r} = \frac{\partial \mathrm v_\theta}{\partial r} \f$ it follows
    !! \f[ \mathrm{rot}(\mathbf v)_\phi =
    !! \frac{\mathrm v_\theta}{r} + \dot\epsilon_{\theta r} - \frac{\dot\epsilon_{r \theta}}{r}
    !! \f]
    !! @todo Check if strain indices i.e. e?? are positioned properly
    !! @todo Why do signs not match up?
    !! @todo Check exz vs. ezx aren't they the same?
    !!
    !! @param self A channel_rot_vy object; usually passed as the dummy object
    !! @param it Current timestep; Position in data array
    !! @param vanderm The Vandermonde matrix used for interpolation
    !! @param ex Element index in theta direction it
    !! @param ex Element index in phi direction it
    !! @param ex Element index in r direction it
    PURE SUBROUTINE record_rot_vy(self, it, vanderm, ex, ey, ez)
        USE elastic_vars_mod, ONLY : ezx, vx, exz
        USE geometric_paras_mod, ONLY : r
        ! Passed dummy arguments
        CLASS(channel_rot_vy), INTENT(INOUT) :: self
        INTEGER, INTENT(IN) :: it, ex, ey, ez
        REAL(real_kind), INTENT(IN) :: vanderm(lpd(config)+1,lpd(config)+1,lpd(config)+1)
        ! Local variables
        REAL(real_kind) :: elm(config%lpd()+1, &
                               config%lpd()+1, &
                               config%lpd()+1)

        elm(:,:,:) = -ezx(ex,ey,ez,:,:,:) + ( vx(ex,ey,ez,:,:,:) - &
            exz(ex,ey,ez,:,:,:) ) / r(ex,ey,ez,:,:,:)

        self%data(it) = SUM( vanderm(:,:,:) * elm(:,:,:) )

    END SUBROUTINE record_rot_vy

    !> Returns the phi component of the curl of the displacement field
    !!
    !! @return field -ezx + ( vx - exz )/r; array of rank 6
    PURE FUNCTION field_rot_vy() RESULT(field)
        USE elastic_vars_mod, ONLY : vx, exz, ezx
        USE geometric_paras_mod, ONLY : r
        REAL(real_kind) :: field(0:config%nx_loc(),0:config%ny_loc(),&
                                 0:config%nz_loc(),0:config%lpd(),&
                                 0:config%lpd(),0:config%lpd())

        field(:,:,:,:,:,:) = -ezx + ( vx - exz )/r

    END FUNCTION


!------------------------------------------------------------------------------
! rot_vz

    !> Records the curl of the displacement-velocity in r direction
    !!
    !! In spherical coordinates the curl's r (i.e. z) component of a vector-field v reads as
    !! \f[ \mathrm{rot}(\mathbf v)_r =
    !! \frac{1}{r\sin\theta} \left(\frac{\partial}{\partial \theta} \left( \mathrm v_\phi\sin\theta \right) - \frac{\partial \mathrm v_\theta}{\partial \phi}\right) =
    !! \frac{1}{r} \left(
    !! \frac{\partial \mathrm v_\phi}{\partial \theta} + \mathrm v_\phi \, \cot\theta -
    !! \frac{1}{\sin\theta}\frac{\partial \mathrm v_\theta}{\partial \phi} \right)
    !! \f]
    !! with \f$ \dot\epsilon_{\phi \theta} = \frac{\partial \mathrm v_\phi}{\partial \theta} \f$
    !! and \f$ \dot\epsilon_{\theta \phi} = \frac{\partial \mathrm v_\theta}{\partial \phi} \f$ it follows
    !! \f[ \mathrm{rot}(\mathbf v)_r =
    !! \frac{1}{r} \left(
    !! \dot\epsilon_{\phi \theta} + \mathrm v_\phi \, \cot\theta -
    !! \frac{\dot\epsilon_{\theta \phi}}{\sin\theta}  \right)
    !! \f]
    !!
    !! @todo Check if strain indices i.e. e?? are positioned properly
    !! @todo Check exz vs. ezx aren't they the same?
    !!
    !! @param self A channel_rot_vz object; usually passed as the dummy object
    !! @param it Current timestep; Position in data array
    !! @param vanderm The Vandermonde matrix used for interpolation
    !! @param ex Element index in theta direction it
    !! @param ex Element index in phi direction it
    !! @param ex Element index in r direction it
    PURE SUBROUTINE record_rot_vz(self, it, vanderm, ex, ey, ez)
        USE elastic_vars_mod, ONLY : vy, exy, eyx
        USE geometric_paras_mod, ONLY : r, cot_theta, sin_theta
        ! Passed dummy arguments
        CLASS(channel_rot_vz), INTENT(INOUT) :: self
        INTEGER, INTENT(IN) :: it, ex, ey, ez
        REAL(real_kind), INTENT(IN) :: vanderm(lpd(config)+1,lpd(config)+1,lpd(config)+1)
        ! Local variables
        REAL(real_kind) :: elm(config%lpd()+1, &
                               config%lpd()+1, &
                               config%lpd()+1)

        elm(:,:,:) = ( exy(ex,ey,ez,:,:,:) + vy(ex,ey,ez,:,:,:)*&
            cot_theta(ex,ey,ez,:,:,:) - eyx(ex,ey,ez,:,:,:)/&
            sin_theta(ex,ey,ez,:,:,:) ) / r(ex,ey,ez,:,:,:)

        self%data(it) = SUM( vanderm(:,:,:) * elm(:,:,:) )

    END SUBROUTINE record_rot_vz

    !> Returns the r component of the curl of the displacement field
    !!
    !! @return field ( exy + vy*cot_theta - eyx/sin_theta )/r; array of rank 6
    PURE FUNCTION field_rot_vz() RESULT(field)
        USE elastic_vars_mod, ONLY : vy, exy, eyx
        USE geometric_paras_mod, ONLY : r, cot_theta, sin_theta
        REAL(real_kind) :: field(0:config%nx_loc(),0:config%ny_loc(),&
                                 0:config%nz_loc(),0:config%lpd(),&
                                 0:config%lpd(),0:config%lpd())

        field(:,:,:,:,:,:) = ( exy + vy*cot_theta - eyx/sin_theta )/r

    END FUNCTION


!------------------------------------------------------------------------------
! pressure

    !> Records pressure rate
    !!
    !! @param self A channel_pressure object; usually passed as the dummy object
    !! @param it Current timestep; Position in data array
    !! @param vanderm The Vandermonde matrix used for interpolation
    !! @param ex Element index in theta direction it
    !! @param ex Element index in phi direction it
    !! @param ex Element index in r direction it
    PURE SUBROUTINE record_pressure(self, it, vanderm, ex, ey, ez)
        USE elastic_vars_mod, ONLY : sxx=>sxx_pml, syy=>syy_pml, szz=>szz_pml
        ! Passed dummy arguments
        CLASS(channel_pressure), INTENT(INOUT) :: self
        INTEGER, INTENT(IN) :: it, ex, ey, ez
        REAL(real_kind), INTENT(IN) :: vanderm(lpd(config)+1,lpd(config)+1,lpd(config)+1)

        self%data(it) = SUM(-1.0/3.0*(sxx(ex,ey,ez,:,:,:) + &
                                      syy(ex,ey,ez,:,:,:) + &
                                      szz(ex,ey,ez,:,:,:))*vanderm)

    END SUBROUTINE record_pressure

    !> Returns field of pressure rate
    !!
    !! @return field -1.0/3.0*(sxx + syy + szz); array of rank 6
    PURE FUNCTION field_pressure() RESULT(field)
        USE elastic_vars_mod, ONLY : sxx=>sxx_pml, syy=>syy_pml, szz=>szz_pml
        REAL(real_kind) :: field(0:config%nx_loc(),0:config%ny_loc(),&
                                 0:config%nz_loc(),0:config%lpd(),&
                                 0:config%lpd(),0:config%lpd())

        field(:,:,:,:,:,:) = -1.0/3.0*(sxx + syy + szz)

    END FUNCTION field_pressure


!------------------------------------------------------------------------------
! nothing

    !> Does simply nothing
    !!
    !! @param self A channel_cls object; usually passed as the dummy object
    !! @param it Current timestep; Position in data array
    !! @param vanderm The Vandermonde matrix used for interpolation
    !! @param ex Element index in theta direction it
    !! @param ex Element index in phi direction it
    !! @param ex Element index in r direction it
    PURE SUBROUTINE record_nothing(self, it, vanderm, ex, ey, ez)
        CLASS(channel_cls), INTENT(INOUT) :: self
        INTEGER, INTENT(IN) :: it, ex, ey, ez
        REAL(real_kind), INTENT(IN) :: vanderm(lpd(config)+1,lpd(config)+1,lpd(config)+1)

    END SUBROUTINE record_nothing

    !> Does simply nothing
    !!
    !! @param self A channel_cls object; usually passed as the dummy object
    !! @param it Current timestep; Position in data array
    !! @param vanderm The Vandermonde matrix used for interpolation
    !! @param ex Element index in theta direction it
    !! @param ex Element index in phi direction it
    !! @param ex Element index in r direction it
    SUBROUTINE add_nothing(self, it, vanderm, ex, ey, ez)
        CLASS(channel_cls), INTENT(IN) :: self
        INTEGER, INTENT(IN) :: it, ex, ey, ez
        REAL(real_kind), INTENT(IN) :: vanderm(lpd(config)+1,lpd(config)+1,lpd(config)+1)

    END SUBROUTINE add_nothing

    !> Returns an empty field
    PURE FUNCTION field_nothing() RESULT(field)
        REAL(real_kind) :: field(0:config%nx_loc(),0:config%ny_loc(),&
                                 0:config%nz_loc(),0:config%lpd(),&
                                 0:config%lpd(),0:config%lpd())

    END FUNCTION field_nothing

!------------------------------------------------------------------------------
! mxx

    SUBROUTINE add_mxx(self, it, vanderm, ex, ey, ez)
        USE elastic_vars_mod, ONLY : src_xx
        USE geometric_paras_mod, ONLY : jac, rr_sin_theta, wx_wy_wz
        CLASS(channel_mxx), INTENT(IN) :: self
        INTEGER, INTENT(IN) :: it, ex, ey, ez
        REAL(real_kind), INTENT(IN) :: vanderm(lpd(config)+1,lpd(config)+1,lpd(config)+1)

        src_xx(ex,ey,ez,:,:,:) = -self%data(it)*vanderm(:,:,:) / &
            (jac*wx_wy_wz(ex,ey,ez,:,:,:)*rr_sin_theta(ex,ey,ez,:,:,:))

    END SUBROUTINE add_mxx

    PURE SUBROUTINE record_mxx(self, it, vanderm, ex, ey, ez)
        USE elastic_vars_mod, ONLY : src_xx
        CLASS(channel_mxx), INTENT(INOUT) :: self
        INTEGER, INTENT(IN) :: it, ex, ey, ez
        REAL(real_kind), INTENT(IN) :: vanderm(lpd(config)+1,lpd(config)+1,lpd(config)+1)

        self%data(it) = SUM( vanderm(:,:,:) * src_xx(ex,ey,ez,:,:,:) )

    END SUBROUTINE

!------------------------------------------------------------------------------

    SUBROUTINE add_mxy(self, it, vanderm, ex, ey, ez)
        USE elastic_vars_mod, ONLY : src_xy
        USE geometric_paras_mod, ONLY : jac, rr_sin_theta, wx_wy_wz
        CLASS(channel_mxy), INTENT(IN) :: self
        INTEGER, INTENT(IN) :: it, ex, ey, ez
        REAL(real_kind), INTENT(IN) :: vanderm(lpd(config)+1,lpd(config)+1,lpd(config)+1)

        src_xy(ex,ey,ez,:,:,:) = -self%data(it)*vanderm(:,:,:) / &
            (jac*wx_wy_wz(ex,ey,ez,:,:,:)*rr_sin_theta(ex,ey,ez,:,:,:))

    END SUBROUTINE add_mxy

    PURE SUBROUTINE record_mxy(self, it, vanderm, ex, ey, ez)
        USE elastic_vars_mod, ONLY : src_xy
        CLASS(channel_mxy), INTENT(INOUT) :: self
        INTEGER, INTENT(IN) :: it, ex, ey, ez
        REAL(real_kind), INTENT(IN) :: vanderm(lpd(config)+1,lpd(config)+1,lpd(config)+1)

        self%data(it) = SUM( vanderm(:,:,:) * src_xy(ex,ey,ez,:,:,:) )

    END SUBROUTINE

!------------------------------------------------------------------------------

    SUBROUTINE add_mxz(self, it, vanderm, ex, ey, ez)
        USE elastic_vars_mod, ONLY : src_xz
        USE geometric_paras_mod, ONLY : jac, rr_sin_theta, wx_wy_wz
        CLASS(channel_mxz), INTENT(IN) :: self
        INTEGER, INTENT(IN) :: it, ex, ey, ez
        REAL(real_kind), INTENT(IN) :: vanderm(lpd(config)+1,lpd(config)+1,lpd(config)+1)

        src_xz(ex,ey,ez,:,:,:) = -self%data(it)*vanderm(:,:,:) / &
            (jac*wx_wy_wz(ex,ey,ez,:,:,:)*rr_sin_theta(ex,ey,ez,:,:,:))

    END SUBROUTINE add_mxz

    PURE SUBROUTINE record_mxz(self, it, vanderm, ex, ey, ez)
        USE elastic_vars_mod, ONLY : src_xz
        CLASS(channel_mxz), INTENT(INOUT) :: self
        INTEGER, INTENT(IN) :: it, ex, ey, ez
        REAL(real_kind), INTENT(IN) :: vanderm(lpd(config)+1,lpd(config)+1,lpd(config)+1)

        self%data(it) = SUM( vanderm(:,:,:) * src_xz(ex,ey,ez,:,:,:) )

    END SUBROUTINE

!------------------------------------------------------------------------------

    SUBROUTINE add_myy(self, it, vanderm, ex, ey, ez)
        USE elastic_vars_mod, ONLY : src_yy
        USE geometric_paras_mod, ONLY : jac, rr_sin_theta, wx_wy_wz
        CLASS(channel_myy), INTENT(IN) :: self
        INTEGER, INTENT(IN) :: it, ex, ey, ez
        REAL(real_kind), INTENT(IN) :: vanderm(lpd(config)+1,lpd(config)+1,lpd(config)+1)

        src_yy(ex,ey,ez,:,:,:) = -self%data(it)*vanderm(:,:,:) / &
            (jac*wx_wy_wz(ex,ey,ez,:,:,:)*rr_sin_theta(ex,ey,ez,:,:,:))

    END SUBROUTINE add_myy

    PURE SUBROUTINE record_myy(self, it, vanderm, ex, ey, ez)
        USE elastic_vars_mod, ONLY : src_yy
        CLASS(channel_myy), INTENT(INOUT) :: self
        INTEGER, INTENT(IN) :: it, ex, ey, ez
        REAL(real_kind), INTENT(IN) :: vanderm(lpd(config)+1,lpd(config)+1,lpd(config)+1)

        self%data(it) = SUM( vanderm(:,:,:) * src_yy(ex,ey,ez,:,:,:) )

    END SUBROUTINE

!------------------------------------------------------------------------------

    SUBROUTINE add_myz(self, it, vanderm, ex, ey, ez)
        USE elastic_vars_mod, ONLY : src_yz
        USE geometric_paras_mod, ONLY : jac, rr_sin_theta, wx_wy_wz
        CLASS(channel_myz), INTENT(IN) :: self
        INTEGER, INTENT(IN) :: it, ex, ey, ez
        REAL(real_kind), INTENT(IN) :: vanderm(lpd(config)+1,lpd(config)+1,lpd(config)+1)

        src_yz(ex,ey,ez,:,:,:) = -self%data(it)*vanderm(:,:,:) / &
            (jac*wx_wy_wz(ex,ey,ez,:,:,:)*rr_sin_theta(ex,ey,ez,:,:,:))

    END SUBROUTINE

    PURE SUBROUTINE record_myz(self, it, vanderm, ex, ey, ez)
        USE elastic_vars_mod, ONLY : src_yz
        CLASS(channel_myz), INTENT(INOUT) :: self
        INTEGER, INTENT(IN) :: it, ex, ey, ez
        REAL(real_kind), INTENT(IN) :: vanderm(lpd(config)+1,lpd(config)+1,lpd(config)+1)

        self%data(it) = SUM( vanderm(:,:,:) * src_yz(ex,ey,ez,:,:,:) )

    END SUBROUTINE

!------------------------------------------------------------------------------

    SUBROUTINE add_mzz(self, it, vanderm, ex, ey, ez)
        USE elastic_vars_mod, ONLY : src_zz
        USE geometric_paras_mod, ONLY : jac, rr_sin_theta, wx_wy_wz
        CLASS(channel_mzz), INTENT(IN) :: self
        INTEGER, INTENT(IN) :: it, ex, ey, ez
        REAL(real_kind), INTENT(IN) :: vanderm(lpd(config)+1,lpd(config)+1,lpd(config)+1)

        src_zz(ex,ey,ez,:,:,:) = -self%data(it)*vanderm(:,:,:) / &
            (jac*wx_wy_wz(ex,ey,ez,:,:,:)*rr_sin_theta(ex,ey,ez,:,:,:))

    END SUBROUTINE

    PURE SUBROUTINE record_mzz(self, it, vanderm, ex, ey, ez)
        USE elastic_vars_mod, ONLY : src_zz
        CLASS(channel_mzz), INTENT(INOUT) :: self
        INTEGER, INTENT(IN) :: it, ex, ey, ez
        REAL(real_kind), INTENT(IN) :: vanderm(lpd(config)+1,lpd(config)+1,lpd(config)+1)

        self%data(it) = SUM( vanderm(:,:,:) * src_zz(ex,ey,ez,:,:,:) )

    END SUBROUTINE

!------------------------------------------------------------------------------
! fx

    !> Records the field fx which corresponds to the source signature of a
    !! a monopole force acting in theta direction.
    !!
    !! @param self A channel_fz object; usually passed as the dummy object
    !! @param it Current timestep; Position in data array
    !! @param vanderm The Vandermonde matrix used for interpolation
    !! @param ex Element index in theta direction it
    !! @param ex Element index in phi direction it
    !! @param ex Element index in r direction it
    PURE SUBROUTINE record_fx(self, it, vanderm, ex, ey, ez)
        USE elastic_vars_mod, ONLY : src_x
        CLASS(channel_fx), INTENT(INOUT) :: self
        INTEGER, INTENT(IN) :: it, ex, ey, ez
        REAL(real_kind), INTENT(IN) :: vanderm(lpd(config)+1,lpd(config)+1,lpd(config)+1)

        self%data(it) = SUM( vanderm(:,:,:) * src_x(ex,ey,ez,:,:,:) )

    END SUBROUTINE record_fx

    SUBROUTINE add_fx(self, it, vanderm, ex, ey, ez)
        USE elastic_vars_mod, ONLY : src_x
        CLASS(channel_fx), INTENT(IN) :: self
        INTEGER, INTENT(IN) :: it, ex, ey, ez
        REAL(real_kind), INTENT(IN) :: vanderm(lpd(config)+1,lpd(config)+1,lpd(config)+1)

       src_x(ex,ey,ez,:,:,:) = -self%data(it) * vanderm(:,:,:)

    END SUBROUTINE add_fx

!------------------------------------------------------------------------------
! fy

    !> Records the field fy which corresponds to the source signature of a
    !! a monopole force acting in phi direction.
    !!
    !! @param self A channel_fz object; usually passed as the dummy object
    !! @param it Current timestep; Position in data array
    !! @param vanderm The Vandermonde matrix used for interpolation
    !! @param ex Element index in theta direction it
    !! @param ex Element index in phi direction it
    !! @param ex Element index in r direction it
    PURE SUBROUTINE record_fy(self, it, vanderm, ex, ey, ez)
        USE elastic_vars_mod, ONLY : src_y
        CLASS(channel_fy), INTENT(INOUT) :: self
        INTEGER, INTENT(IN) :: it, ex, ey, ez
        REAL(real_kind), INTENT(IN) :: vanderm(lpd(config)+1,lpd(config)+1,lpd(config)+1)

        self%data(it) = SUM( vanderm(:,:,:) * src_y(ex,ey,ez,:,:,:) )

    END SUBROUTINE record_fy

    SUBROUTINE add_fy(self, it, vanderm, ex, ey, ez)
        USE elastic_vars_mod, ONLY : src_y
        CLASS(channel_fy), INTENT(IN) :: self
        INTEGER, INTENT(IN) :: it, ex, ey, ez
        REAL(real_kind), INTENT(IN) :: vanderm(lpd(config)+1,lpd(config)+1,lpd(config)+1)

       src_y(ex,ey,ez,:,:,:) = -self%data(it) * vanderm(:,:,:)

    END SUBROUTINE add_fy

!------------------------------------------------------------------------------
! fz

    !> Records the field fz which corresponds to the source signature of a
    !! a monopole force acting in r direction.
    !!
    !! @param self A channel_fz object; usually passed as the dummy object
    !! @param it Current timestep; Position in data array
    !! @param vanderm The Vandermonde matrix used for interpolation
    !! @param ex Element index in theta direction it
    !! @param ex Element index in phi direction it
    !! @param ex Element index in r direction it
    PURE SUBROUTINE record_fz(self, it, vanderm, ex, ey, ez)
        USE elastic_vars_mod, ONLY : src_z
        CLASS(channel_fz), INTENT(INOUT) :: self
        INTEGER, INTENT(IN) :: it, ex, ey, ez
        REAL(real_kind), INTENT(IN) :: vanderm(lpd(config)+1,lpd(config)+1,lpd(config)+1)

        self%data(it) = SUM( vanderm(:,:,:) * src_z(ex,ey,ez,:,:,:) )

    END SUBROUTINE record_fz

    SUBROUTINE add_fz(self, it, vanderm, ex, ey, ez)
        USE elastic_vars_mod, ONLY : src_z
        CLASS(channel_fz), INTENT(IN) :: self
        INTEGER, INTENT(IN) :: it, ex, ey, ez
        REAL(real_kind), INTENT(IN) :: vanderm(lpd(config)+1,lpd(config)+1,lpd(config)+1)

       src_z(ex,ey,ez,:,:,:) = -self%data(it) * vanderm(:,:,:)

    END SUBROUTINE add_fz


!------------------------------------------------------------------------------
! exx

    !> Returns theta-theta strain
    !!
    !! @return field dxux; array of rank 6
    PURE FUNCTION field_exx() RESULT(field)
        USE elastic_vars_mod, ONLY : dxux
        REAL(real_kind) :: field(0:config%nx_loc(),0:config%ny_loc(),&
                                 0:config%nz_loc(),0:config%lpd(),&
                                 0:config%lpd(),0:config%lpd())
        field(:,:,:,:,:,:) = dxux
    END FUNCTION


!------------------------------------------------------------------------------
! exy

    !> Returns theta-phi strain
    !!
    !! @todo Actually, symmetrization should not be necessary
    !!
    !! @return field ( dxuy + dyux ) / 2.0; array of rank 6
    PURE FUNCTION field_exy() RESULT(field)
        USE elastic_vars_mod, ONLY : dxuy, dyux
        REAL(real_kind) :: field(0:config%nx_loc(),0:config%ny_loc(),&
                                 0:config%nz_loc(),0:config%lpd(),&
                                 0:config%lpd(),0:config%lpd())
        field(:,:,:,:,:,:) = ( dxuy + dyux ) / 2.0
    END FUNCTION

!------------------------------------------------------------------------------
! exz

    !> Returns theta-r strain
    !!
    !! @todo Actually, symmetrization should not be necessary
    !!
    !! @return field ( dxuz + dzux ) / 2.0; array of rank 6
    PURE FUNCTION field_exz() RESULT(field)
        USE elastic_vars_mod, ONLY : dxuz, dzux
        REAL(real_kind) :: field(0:config%nx_loc(),0:config%ny_loc(),&
                                 0:config%nz_loc(),0:config%lpd(),&
                                 0:config%lpd(),0:config%lpd())
        field(:,:,:,:,:,:) = ( dxuz + dzux ) / 2.0
    END FUNCTION

!------------------------------------------------------------------------------
! eyy

    !> Returns phi-phi strain
    !!
    !! @return field dyuy ; array of rank 6
    PURE FUNCTION field_eyy() RESULT(field)
        USE elastic_vars_mod, ONLY : dyuy
        REAL(real_kind) :: field(0:config%nx_loc(),0:config%ny_loc(),&
                                 0:config%nz_loc(),0:config%lpd(),&
                                 0:config%lpd(),0:config%lpd())
        field(:,:,:,:,:,:) = dyuy
    END FUNCTION


!------------------------------------------------------------------------------
! eyz

    !> Returns phi-r strain
    !!
    !! @todo Actually, symmetrization should not be necessary
    !!
    !! @return field ( dyuz + dzuy ) / 2.0; array of rank 6
    PURE FUNCTION field_eyz() RESULT(field)
        USE elastic_vars_mod, ONLY : dyuz, dzuy
        REAL(real_kind) :: field(0:config%nx_loc(),0:config%ny_loc(),&
                                 0:config%nz_loc(),0:config%lpd(),&
                                 0:config%lpd(),0:config%lpd())
        field(:,:,:,:,:,:) = ( dyuz + dzuy ) / 2.0
    END FUNCTION


!------------------------------------------------------------------------------
! ezz

    !> Returns r-r strain
    !!
    !! @return field dzuz; array of rank 6
    PURE FUNCTION field_ezz() RESULT(field)
        USE elastic_vars_mod, ONLY : dzuz
        REAL(real_kind) :: field(0:config%nx_loc(),0:config%ny_loc(),&
                                 0:config%nz_loc(),0:config%lpd(),&
                                 0:config%lpd(),0:config%lpd())
        field(:,:,:,:,:,:) = dzuz
    END FUNCTION


END MODULE channel_mod
