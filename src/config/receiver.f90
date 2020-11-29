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
!! $Date: 2013-12-04 08:50:52 +0100 (Wed, 04 Dec 2013) $
!! $Author: mauerberger $
!! $Revision: 838 $
!! @copyright GNU General Public License version 3 or later



!> This module is dedicated to handle Ses3d-NT's synthetic seismograms. It
!! provides a class receiver_cls which extends point_mod::point_cls.
!!
!! Based upon Fortran's NAMELIST-feature the function parse_receiver_file()
!! is provided which parses multiple '&receiver'-groups specified in a file.
!! In return it gives an array of receiver_cls objects filled with values from
!! the receiver file.
MODULE receiver_mod
    USE error_mod, ONLY : warn, abort
    USE parameters_mod, ONLY : real_kind, fnl, &
                               mpi_master_rank, is_mpi_master_rank
    USE configuration_mod, ONLY : config => configuration
    USE point_mod, ONLY : point_cls
    USE channel_mod, ONLY : channel_container

    IMPLICIT NONE

    PRIVATE

    !> Class handling receivers
    !!
    !! Extends point_mod::point_cls
    !!
    !! @todo ... doc-string missing ...
    TYPE, EXTENDS(point_cls) :: receiver_cls
        !> The filename prefix. A directory where traces are stored in.
        CHARACTER(LEN=fnl), PRIVATE :: prefix
        !> Seismic network name
        CHARACTER(LEN=8), PRIVATE :: network
        !> Station name
        CHARACTER(LEN=8), PRIVATE :: station
        !> Station location
        CHARACTER(LEN=8), PRIVATE :: location
        !> A flag which says whether existing files are overridden or not
        LOGICAL, PRIVATE :: override
        !> An array of channels objects channel_mod::channel_container
        TYPE(channel_container), ALLOCATABLE, PRIVATE :: channels(:)
    CONTAINS
        !> Writes some information about the object into a unit specifier.
        !! For details see receiver_mod::output_header().
        !! Overloads the inherited procedure point_mod::point_cls::print()
        PROCEDURE, PUBLIC :: print => output_header
        !> Records the values of channels at the current timestep
        !! @todo ... doc-string missing ...
        !! For details see receiver_mod::record().
        PROCEDURE, PUBLIC :: record
        !> Records the values of channels at the current timestep
        !! @todo ... doc-string missing ...
        !! For details see receiver_mod::write().
        PROCEDURE, PUBLIC :: write
        !> Records the values of channels at the current timestep
        !! @todo ... doc-string missing ...
        !! For details see receiver_mod::id().
        PROCEDURE, PUBLIC :: id
        !> Returns .TRUE. if the receiver_cls object is not initialized.
        !! @warning 'Uninitialized' just says that the object is precisely the
        !!          same as receiver_mod::uninitialized_receiver_cls().
        !! For details see receiver_mod::uninitialized().
        !! Overloads the inherited procedure point_mod::point_cls::uninitialized()
        PROCEDURE, PUBLIC :: uninitialized
        !> For comparing two receiver_cls objects. Do not use this procedure
        !! explicitly. Use the comparison operator '==' instead.
        !! For details see receiver_mod::equal_receivers().
        !! Overloads the inherited procedure point_mod::point_cls::equal_points()
        PROCEDURE, PRIVATE :: equal => equal_receivers
    END TYPE receiver_cls

    ! An interface ...
    INTERFACE receiver_cls
        !> Constructor function ...
        !! @todo ... doc-string missing ...
        !! For details see receiver_mod::create_receiver().
        MODULE PROCEDURE create_receiver
        !> A constructor expecting no arguments and returning an uninitialized
        !! receiver_cls object. That constructor is used if something goes
        !! wrong. For details see receiver_mod::uninitialized_receiver_cls().
        MODULE PROCEDURE :: uninitialized_receiver_cls
    END INTERFACE

    ! Expose the class and its constructors accompanied by the parser function
    PUBLIC :: receiver_cls, parse_receiver_file

CONTAINS

!------------------------------------------------------------------------------
! receiver_cls class bound procedures


    !> Returns a SEED compatible identifier which contains the network,
    !! station, location and channel code for the reveiver_cls object.
    !!
    !! @param self Passed-object dummy of class receiver_cls
    !! @return A SEED compatible identifier of the trace
    PURE FUNCTION id(self)
        CLASS(receiver_cls), INTENT(IN) :: self
        CHARACTER(LEN=LEN_TRIM(self%network) + 1 + &
                      LEN_TRIM(self%station) + 1 + &
                      LEN_TRIM(self%location)) :: id
        id = TRIM(self%network) // '.' // &
             TRIM(self%station) // '.' //&
             TRIM(self%location)
    END FUNCTION id


    !> Subroutine to record synthetic seismograms
    !!
    !! @todo ... doc-string missing ...
    !!
    !! @param it Tim-step/position
    !! @param self Passed-object dummy of class receiver_cls
    ELEMENTAL SUBROUTINE record(self, it)
        CLASS(receiver_cls), INTENT(INOUT) :: self
        INTEGER, INTENT(IN) :: it

        CALL self%channels(:)%record(it=it, point=self%point_cls )

    END SUBROUTINE record


    !> Writes SAC traces
    !!
    !! @todo ... doc-string missing ...
    !! @todo: Well, we should better have a 'format'-parameter to select between
    !!        different file formats
    !! @todo: Make sure that only one rank is writing at the same time
    !!       not to spread out many writing file accesses in the same directory
    !!
    !! @param self Passed-object dummy of class receiver_cls
    SUBROUTINE write(self)
        USE sac_io_mod, ONLY : sac_cls
        USE parameters_mod, ONLY : isl
        CLASS(receiver_cls), INTENT(INOUT) :: self

        CHARACTER(LEN=isl) :: iomsg
        TYPE(sac_cls) :: sac_obj
        INTEGER :: i, iostat

        ! The actual code writing SAC files
        DO i=1, SIZE(self%channels)
            ! Create SAC-object
            sac_obj = sac_cls(lat=self%lat(), &
                              lon=self%lon(), &
                              depth=self%depth(), &
                              dt=config%dt(), &
                              data=self%channels(i)%data(), &
                              network=self%network, &
                              station=self%station, &
                              location=self%location, &
                              channel=self%channels(i)%attribute(), &
                              event=config%event_name(), &
                              date_time=config%date_time(), &
                              override=self%override)
            ! Write SAC trace
            CALL sac_obj%write(TRIM(self%prefix)//sac_obj%auto_file_name(), &
                               iostat=iostat, iomsg=iomsg)
            ! Write status to screen and log-file
            IF (iostat /= 0) THEN
                CALL warn(iomsg)
                IF (config%log_file_true()) THEN
                    CALL warn(msg=iomsg, unit=config%log_file_unit)
                END IF
            ELSE
                WRITE(UNIT=config%log_file_unit, FMT='(A)') TRIM(iomsg)
            END IF
        END DO

    END SUBROUTINE write


    !> This subroutine writes various information about receivers into a passed
    !! unit-specifier
    !!
    !! @param self Passed-object dummy of class receiver_cls
    !! @param unit An already opened unit specifier to write into
    SUBROUTINE output_header(self, unit)
        CLASS(receiver_cls), INTENT(IN) :: self
        INTEGER, INTENT(IN) :: unit

        WRITE( UNIT=unit, FMT='(A)' ) 'Receiver: '

        ! Naming
        WRITE( UNIT=unit, FMT='(4x,A," = ",A)' ) 'network', self%network
        WRITE( UNIT=unit, FMT='(4x,A," = ",A)' ) 'station', self%station
        WRITE( UNIT=unit, FMT='(4x,A," = ",A)' ) 'location', self%location

        ! Channels
        WRITE( UNIT=unit, FMT=* ) '   channels: ', self%channels(:)%attribute()

        ! Use the point_cls' print procedure to write coordinates
        CALL self%point_cls%print(unit=unit)

    END SUBROUTINE output_header


    !> Compares two objects p1 and p2 of class receiver_cls. In case these
    !! objects are equal `.TRUE.` is returned otherwise `.FALSE.`.
    !!
    !! Invokes the comparison operator '==' of point_mod::point_cls.
    !! In addition it uses the comparison operator '==' for objects of type
    !! channel_mod::channel_container.
    !!
    !! The function is used to overload the inherited comparison operator of
    !! class receiver_cls.
    !!
    !! @note Two objects p1 and p2 are called identical or equal if all their
    !!       components are equal:
    !!
    !!  * The base objects point_mod::point_cls are equal
    !!  * prefix: The fine-name prefix
    !!  * override: The flag which says whether existing files are overridden
    !!  * network: Seismic network name
    !!  * station: Station name
    !!  * location: Station location
    !!  * All channel_mod::channel_container objects are equal
    !!
    !! @param p1 Passed-object dummy of class receiver_cls
    !! @param p2 Passed-object dummy of class receiver_cls
    !! @return .TRUE. if p1 and p2 are identical otherwise .FALSE.
    ELEMENTAL LOGICAL FUNCTION equal_receivers(p1, p2)
        ! Passed dummy arguments
        CLASS(receiver_cls), INTENT(IN) :: p1
        CLASS(point_cls), INTENT(IN) :: p2

        equal_receivers = .FALSE.

        ! If p2 is not of class source_cls return
        ! Because p2 may hold any extension of point_cls
        IF ( .NOT. SAME_TYPE_AS(p1, p2) ) RETURN
        ! Makes p2's type explicit; It gives access to p2's components
        SELECT TYPE(p2)
            CLASS IS (receiver_cls) ! Guaranteed to execute since same type
                IF ( .NOT.(p1%point_cls==p2%point_cls) ) RETURN
                IF ( ALLOCATED(p1%channels) .NEQV. ALLOCATED(p2%channels) ) RETURN
                IF ( SIZE(p1%channels) /= SIZE(p2%channels) ) RETURN
                IF ( .NOT. ALL(p1%channels==p2%channels) ) RETURN
                IF ( p1%prefix/=p2%prefix ) RETURN
                IF ( p1%network/=p2%network ) RETURN
                IF ( p1%station/=p2%station) RETURN
                IF ( p1%location/=p2%location) RETURN
                IF ( p1%override.NEQV.p2%override) RETURN
        END SELECT

        equal_receivers = .TRUE.

    END FUNCTION


    !> Checks if an object of receiver_cls class is initialized or not
    !!
    !! @note A receiver - an object of class receiver_cls - is called
    !!       uninitialized if it was instantiated without passing any arguments
    !!       (i.e. `... = receiver_cls()`). For details see:
    !!       receiver_mod::uninitialized_point_cls().
    !!
    !! @warning The converse 'not uninitialized' does not mean that an object
    !!          is already properly initialized.
    !!
    !! @param self Passed-object dummy of class receiver_cls
    !! @return .TRUE. if the passed object is uninitialized otherwise .FALSE.
    ELEMENTAL LOGICAL FUNCTION uninitialized(self)
        ! Passed dummy arguments
        CLASS(receiver_cls), INTENT(IN) :: self

        IF ( equal_receivers(self, receiver_cls()) ) THEN
            uninitialized = .TRUE.
        ELSE
            uninitialized = .FALSE.
        END IF

    END FUNCTION


!------------------------------------------------------------------------------
! receiver_cls structure constructors

    !> Structure constructor returning an uninitialized/empty
    !! receiver_cls object
    !!
    !! Its instantiated using an empty point_cls object ...  an zero sized
    !! channel_container array.
    !!
    !! @note An object of class receiver_cls is called uninitialized if it is
    !!       instantiated using that very constructor.
    !!
    !! @return Uninitialized object of type receiver_cls
    ELEMENTAL FUNCTION uninitialized_receiver_cls() RESULT(new)
        ! Result variable
        TYPE(receiver_cls) :: new
        ! Local variables
        TYPE(channel_container) :: channels(0)

        ! Initialize object to 'null' invoking Fortran's default constructor
        !new = receiver_cls( point_cls=point_cls(), prefix='', &
        !                    network='', station='', location='', &
        !                    override=.FALSE., channels=channels )
        new%point_cls=point_cls()
        new%prefix=''
        new%network=''
        new%station=''
        new%location=''
        new%override=.FALSE.
        ALLOCATE(new%channels(0), source=channels)

    END FUNCTION


    !> receiver_cls constructor function ...
    !!
    !! @todo ... doc-string missing ...
    !!
    !! @return An object of type receiver_cls
    FUNCTION create_receiver(rec_in) RESULT(new)
        ! Use associated entities
        USE parser_mod, ONLY : receiver_typ
        USE channel_mod, ONLY : channel_container
        USE string_utilities_mod, ONLY : strip, array_to_string, lower_case
        ! Passed dummy arguments
        CLASS(receiver_typ), INTENT(IN) :: rec_in
        ! Result variable
        TYPE(receiver_cls) :: new
        ! Local variables
        CHARACTER(LEN=8) :: network
        CHARACTER(LEN=8), ALLOCATABLE :: attributes_(:)
        TYPE(channel_container), ALLOCATABLE :: channels_(:)
        TYPE(point_cls) :: point

        ! Invoke the point-class structure constructor
        point = point_cls(lat=rec_in%lat, lon=rec_in%lon, depth=rec_in%depth)
        IF ( point%uninitialized() ) THEN
            new = receiver_cls()
            RETURN
        END IF

        ! do not accept empty network or station
        IF ( rec_in%network == '' ) THEN
            CALL warn( "Got invalid value="//TRIM(rec_in%network)//&
                       " as network name" )
            network = 'Ses3d-NT'
        ELSE
            network = rec_in%network
        END IF

        IF ( rec_in%station == '' ) THEN
            CALL warn( "Got invalid value="//TRIM(rec_in%station)//&
                       " as station name" )
            new = receiver_cls()
            ! TODO write warning to log-file
            RETURN
        END IF

        ! TODO prefix has to end with a slash

        ! Strip, allocate and assign local attributes field
        ! Empty entities are dropped; all strings are aligned left
        attributes_ = strip(rec_in%attributes(:))
        ! Make all attributes lower case
        attributes_(:) = lower_case(attributes_(:))

        ! Allocate local channels class
        ALLOCATE( channels_(SIZE(attributes_(:))) )
        ! Initialize and assign channels
        channels_(:) = channel_container(attributes_)

        ! Warn about unknown attributes
        attributes_ = PACK( channels_(:)%attribute(), channels_(:)%unknown() )
        IF ( SIZE(attributes_) /= 0 ) &
            CALL WARN( "Unknown attributes " // array_to_string(attributes_(:)) &
                // " in receiver '" // new%id() // "'."  )

        ! Drop unknown attributes
        !DEALLOCATE( new%channels )
        ALLOCATE( new%channels( COUNT(.NOT. channels_(:)%unknown()) ) )
        new%channels(:) = PACK( channels_(:), .NOT. channels_(:)%unknown())

        new%point_cls = point
        new%prefix = rec_in%prefix
        new%network = network
        new%station = rec_in%station
        new%location = rec_in%location
        new%override = rec_in%override

        !new = receiver_cls( point_cls=point, &
        !                    prefix=rec_in%prefix, &
        !                    network=rec_in%network, &
        !                    station=rec_in%station, &
        !                    location=rec_in%location, &
        !                    override=rec_in%override, &
        !                    channels=channels_ )

        ! Warn if no channels are specified
        IF ( SIZE(new%channels(:)) == 0 .AND. is_mpi_master_rank() ) THEN
            CALL WARN( "No attributes specified in receiver "//new%id() )
            new = receiver_cls()
            RETURN
        END IF

    END FUNCTION create_receiver


!------------------------------------------------------------------------------
! Public module procedures


    !> Subroutine parsing all occurring &receiver NAMELIST groups in a file
    !!
    !! @todo ... doc-string missing ...
    !!
    !! @param receiver_file_name File name
    !! @param receiver An unallocated array of receiver_cls types
    SUBROUTINE parse_receiver_file( receiver_file_name, receivers )
        USE parser_mod, ONLY : count_namelists, parse_namelist, &
                               receiver_typ
        USE bcast_mod, ONLY : bcast_namelist, bcast_integer, reduce_sum
        ! Passed dummy arguments
        CHARACTER(LEN=*), INTENT(IN) :: receiver_file_name
        TYPE( receiver_cls ), INTENT(OUT), ALLOCATABLE :: receivers(:)
        ! Local variables
        INTEGER :: receiver_file_unit
        TYPE( receiver_typ) :: receiver
        INTEGER :: n_recs, i

        ! Open receiver file
        IF ( is_mpi_master_rank() ) &
            OPEN( NEWUNIT=receiver_file_unit, FILE=receiver_file_name, &
                  ACTION='READ', STATUS='OLD' )

        !----------------
        ! Parse receivers
        !----------------
        IF ( is_mpi_master_rank() ) &
            ! Count number of NAMELIST groups
            n_recs = count_namelists( unit=receiver_file_unit, nml_cls=receiver)

        ! Broadcast n_recs amongst all ranks
        CALL bcast_integer(i=n_recs, mr=mpi_master_rank)

        ! Allocate array of receivers
        ALLOCATE( receivers(n_recs) )

        IF ( is_mpi_master_rank() ) &
            REWIND( UNIT=receiver_file_unit )

        ! Iterate through all &receiver groups
        DO i=1, n_recs
            IF ( is_mpi_master_rank() ) THEN
                ! Parse NAMELIST
                CALL parse_namelist(unit=receiver_file_unit, nml_cls=receiver)
            END IF
            ! Broadcast values
            CALL bcast_namelist(receiver, mpi_master_rank)
            ! Invoke structure constructor
            receivers(i) = receiver_cls(receiver)
        END DO

        ! all done -> close file
        IF ( is_mpi_master_rank() ) &
            CLOSE( receiver_file_unit )

        ! Sort out all receivers which are not located in this rank
        receivers = PACK( receivers(:), .NOT. receivers(:)%uninitialized() )

        ! Check if there is at least one receiver station left
        n_recs = reduce_sum(i=SIZE(receivers), mr=mpi_master_rank)
        IF ( n_recs < 1 .AND. is_mpi_master_rank() ) &
            CALL warn( 'No receivers defined' )

    END SUBROUTINE parse_receiver_file


END MODULE receiver_mod
