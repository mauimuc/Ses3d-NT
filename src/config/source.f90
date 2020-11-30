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



!> Module to handle multiple sources
!!
!! This module is dedicated to handle Ses3d-NT's monopole-sources as well as
!! moment-tensor sources.
MODULE source_mod
    USE parameters_mod, ONLY : real_kind, fnl, &
                               mpi_master_rank, is_mpi_master_rank
    USE configuration_mod, ONLY : config => configuration
    USE point_mod, ONLY : point_cls
    USE channel_mod, ONLY : channel_container

    IMPLICIT NONE

    PRIVATE

    !> Class handling multiple sources ...
    !!
    !! Extends point_mod::point_cls
    !!
    !! @todo ... doc-string missing ...
    TYPE, EXTENDS(point_cls) :: source_cls
        !> An array of channels objects channel_mod::channel_container
        TYPE(channel_container), ALLOCATABLE, PRIVATE :: channels(:)
    CONTAINS
        !> Writes some information about the object into a unit specifier.
        !! For details see source_mod::output_header().
        !! Overloads the procedure inherited from point_mod::point_cls::print
        PROCEDURE, PUBLIC :: print => output_header
        !> Adds ...
        !! @todo ... doc-string missing ...
        !! For details see source_mod::add()
        PROCEDURE, PUBLIC :: add
        !> Returns .TRUE. if the source_cls object is not initialized.
        !! @warning 'Uninitialized' just says that the object is precisely the
        !!          same as source_mod::uninitialized_source_cls().
        !! For details see source_mod::uninitialized().
        !! Overloads the inherited procedure point_mod::point_cls::uninitialized
        PROCEDURE, PUBLIC :: uninitialized
        !> For comparing two source_cls objects. Do not use this procedure
        !! explicitly. Use the comparison operator '==' instead.
        !! For details see source_mod::equal_sources().
        !! Overloads the inherited procedure point_mod::point_cls::equal_points
        PROCEDURE, PRIVATE :: equal => equal_sources
    END TYPE source_cls

    ! Interface
    INTERFACE source_cls
        !> Create from monopole force ...
        !! @todo ... doc-string missing ...
        !! For details see source_mod::create_from_source_sf().
        MODULE PROCEDURE create_from_source_sf
        !> Create from dipole force ...
        !! @todo ... doc-string missing ...
        !! For details see source_mod::create_from_source_mt().
        MODULE PROCEDURE create_from_source_mt
        !> Returns a source_cls object instantiated from passing an object
        !! of class sac_io_mod::sac_cls.
        !! For details see source_mod::uninitialized_source.
        MODULE PROCEDURE create_from_sac
        !> A constructor expecting no arguments and returning an uninitialized
        !! receiver_cls object. That constructor is used if something goes
        !! wrong. For details see source_mod::uninitialized_source_cls().
        MODULE PROCEDURE :: uninitialized_source_cls
    END INTERFACE

    ! Expose the class and its constructors accompanied by the parser function
    PUBLIC :: source_cls, parse_source_file

CONTAINS

!------------------------------------------------------------------------------
! source_cls structure constructors

    !> Structure constructor returning an uninitialized/empty souce_cls object
    !!
    !! Its instantiated using an empty point_cls object and an zero sized
    !! channel_container array.
    !!
    !! @note An object of class source_cls is called uninitialized if it is
    !!       instantiated using that very constructor.
    !!
    !! @return Uninitialized object of type source_cls
    ELEMENTAL FUNCTION uninitialized_source_cls() RESULT(new)
        ! Result variable
        TYPE(source_cls) :: new
        ! Local variables
        TYPE(channel_container) :: channels(0)

        ! Initialize object to 'null' invoking Fortran's default constructor
        !new = source_cls(point_cls=point_cls(), channels=channels)
        new%point_cls = point_cls()
        ALLOCATE( new%channels(0), source=channels )

    END FUNCTION


    !> Constructor function which instantiates a source_cls object out of a
    !! SAC-trace. The SAC-filename is passed as a parser_mod::source_sac_typ
    !! object.
    !!
    !! The new source_cls object is created by using the following values:
    !!   * Latitude as entry STLA of the SAC-header
    !!   * Longitude as entry STLO of the SAC-header
    !!   * Depth as entry STDP of the SAC-header
    !!   * Channel as entry KSTCMP of the SAC-header
    !!   * Source signature is the dependent data section
    !!   * The whole trace is multiplied by the SCALE SAC-header field
    !!
    !! For details about the SAC file-format see sac_io_mod.
    !!
    !! @note If the SAC-trace is too short the run is aborted
    !! @note If the channel is not know the source will be sorted out
    !! @note A SAC-trace provides just one attribute/channel
    !!
    !! @param src_sac An object of type parser_mod::source_sac_typ
    !! @return An object of type source_cls holding one channel
    FUNCTION create_from_sac(src_sac) RESULT(new)
        ! Use associated entities
        USE error_mod, ONLY : warn, abort
        USE parser_mod, ONLY : source_sac_typ
        USE sac_io_mod, ONLY : sac_cls
        ! Passed dummy arguments
        TYPE(source_sac_typ), INTENT(IN) :: src_sac
        ! Result variable
        TYPE(source_cls) :: new
        ! Local Variables
        CHARACTER(LEN=8) :: attribute
        REAL(real_kind), ALLOCATABLE :: signature(:)
        REAL(real_kind) :: lat, lon, depth, scale
        TYPE(sac_cls) :: sac_obj
        TYPE(point_cls) :: point
        TYPE(channel_container) :: channel

        ! Read SAC-file
        sac_obj = sac_cls(src_sac%file_name)

        ! Get values from header
        lat = sac_obj%get_header('STLA')
        lon = sac_obj%get_header('STLO')
        depth = sac_obj%get_header('STDP')

        ! Invoke the point-class structure constructor
        point = point_cls(lat=lat, lon=lon, depth=depth)
        !! If source is not in MPI rank or something went wrong sort it out
        IF ( point%uninitialized() ) THEN
            ! Return uninitialized sourc_cls object
            new = source_cls()
            RETURN
        END IF

        ! Get channel
        attribute = sac_obj%get_header('KCMPNM')
        ! Get source signature
        scale = sac_obj%get_header('scale')
        signature = sac_obj%get_data1()
        ! Check length of the trace
        IF ( SIZE(signature) < config%nt() ) THEN
            ! todo write warning to log-file
            CALL abort("SAC-trace too short!")
        END IF

        channel = channel_container(attribute=attribute, data=signature*scale)

        ! Warn if unknown attributes
        IF ( channel%unknown() ) THEN
            CALL WARN( "Unknown attribute '" // &
                TRIM(channel%attribute()) // "' in SAC source." )
            ! todo write warning to log-file
            ! Return uninitialized sourc_cls object
            new = source_cls()
        END IF

        ! Initialize object to 'null' invoking Fortran's default constructor
        !new = source_cls(point_cls=point, channels=[channel])
        new%point_cls = point
        ALLOCATE( new%channels(1), source=channel )

    END FUNCTION


    !> Creates a source_cls object from a monopole force
    !! parser_mod::source_sf_typ
    !!
    !! @todo ... doc-string missing ...
    !! @todo Check all parsed values!!!
    !! @todo Sort out unknown or zero channels
    !!
    !! @param src_sf parser_mod::source_sf_typ
    !! @return An object of type source_cls
    FUNCTION create_from_source_sf(src_sf) RESULT(new)
        ! Use associated entities
        USE error_mod, ONLY : warn
        USE parser_mod, ONLY : source_sf_typ, default_real
        ! Passed dummy variables
        TYPE(source_sf_typ), INTENT(IN) :: src_sf
        ! Result variable
        TYPE(source_cls) :: new
        ! Local Variables
        REAL(real_kind) :: data(config%nt())
        TYPE(channel_container) :: channels(3)
        TYPE(point_cls) :: point

        ! Invoke the point-class structure constructor
        point = point_cls(lat=src_sf%lat, lon=src_sf%lon, depth=src_sf%depth)
        ! If something went wrong while constructing the point return 'null'
        IF ( point%uninitialized() ) THEN
            new = source_cls()
            RETURN
        END IF

        ! Check if all components got values assigned
        IF ( ANY(src_sf%direction == default_real) ) THEN
            IF ( is_mpi_master_rank() ) THEN
                CALL warn('Direction vector of monopole force invalid!')
            END IF
            ! todo log-file
            new = source_cls()
            RETURN
        END IF

        ! Make source signature; Parameters arch double-checked there
        data(:) = make_source_signature_(src_sf%wavelet, &
                                         src_sf%onset, &
                                         src_sf%width)

        channels(:) = [ channel_container('fx', data(:)*src_sf%direction(1)), &
                        channel_container('fy', data(:)*src_sf%direction(2)), &
                        channel_container('fz', data(:)*src_sf%direction(3)) ]

        ! todo sort out channels with ALL(data(:) == 0.0)

        ! Initialize object invoking Fortran's default constructor
        !new = source_cls(point_cls=point, channels=channels)
        new%point_cls = point
        ALLOCATE( new%channels(3), source=channels)

    END FUNCTION


    !> Creates a source_cls object from a moment tensor force
    !! parser_mod::source_mt_typ
    !!
    !! @todo ... doc-string missing ...
    !! @todo Check all parsed values!!!
    !! @todo Sort out unknown or zero channels
    !!
    !! @param src_mt parser_mod::source_mt_typ
    !! @return An object of type source_cls
    FUNCTION create_from_source_mt(src_mt) RESULT(new)
        ! Use associated entities
        USE error_mod, ONLY : warn
        USE parser_mod, ONLY : source_mt_typ, default_real
        ! Passed dummy variables
        TYPE(source_mt_typ), INTENT(IN) :: src_mt
        ! Result variable
        TYPE(source_cls) :: new
        ! Local Variables
        INTEGER :: i
        REAL(real_kind) :: data(config%nt())
        TYPE(channel_container) :: channels(6)
        LOGICAL :: mask(6)
        TYPE(point_cls) :: point

        ! Invoke the point-class structure constructor
        point = point_cls(lat=src_mt%lat, lon=src_mt%lon, depth=src_mt%depth)
        ! If something went wrong while constructing the point return 'null'
        IF ( point%uninitialized() ) THEN
            new = source_cls()
            RETURN
        END IF

        ! Check if all components got values assigned
        IF ( ANY(src_mt%moment_tensor == default_real) ) THEN
            IF ( is_mpi_master_rank() ) THEN
                CALL warn('Moment-tensor of dipol force invalid!')
            END IF
            ! todo log-file
            new = source_cls()
            RETURN
        END IF

        data(:) = make_source_signature_(src_mt%wavelet, &
                                         src_mt%onset, &
                                         src_mt%width)

        channels(:) = [channel_container('mxx', data*src_mt%moment_tensor(1)),&
                       channel_container('mxy', data*src_mt%moment_tensor(2)),&
                       channel_container('mxz', data*src_mt%moment_tensor(3)),&
                       channel_container('myy', data*src_mt%moment_tensor(4)),&
                       channel_container('myz', data*src_mt%moment_tensor(5)),&
                       channel_container('mzz', data*src_mt%moment_tensor(6))]

        !new = source_cls(point_cls=point, channels=channels)
        new%point_cls = point
        DO i=1, SIZE(channels)
            mask(i) = .NOT. ALL(channels(i)%data() == 0.0_real_kind )
        END DO
        ALLOCATE( new%channels(COUNT(mask)), source=PACK(channels, mask))

        ! todo check if components are left

    END FUNCTION


!------------------------------------------------------------------------------
! source_cls class bound procedures


    !> This subroutine writes various information about sources into a passed
    !! unit-specifier
    !!
    !! @param self Passed object-dummy of class source_cls
    !! @param unit An already opened unit specifier to write into
    SUBROUTINE output_header(self, unit)
        ! Passed dummy variables
        CLASS(source_cls), INTENT(IN) :: self
        INTEGER, INTENT(IN) :: unit

        WRITE(UNIT=unit, FMT='(A)' ) 'Source: '
        WRITE(UNIT=unit, FMT=*) '   channels: ', self%channels(:)%attribute()

        ! Use the point_cls' print procedure to write coordinates
        CALL self%point_cls%print(unit=unit)

    END SUBROUTINE output_header


    !> Subroutine adding source contributions at time-step 'it' to src-fields
    !!
    !! @todo ... doc-string missing ...
    !!
    !! @param self Passed object-dummy of class source_cls
    !! @param it Current time-step
    SUBROUTINE add(self, it)
        ! Passed dummy variables
        CLASS(source_cls), INTENT(IN) :: self
        INTEGER, INTENT(IN) :: it
        ! Local Variables
        INTEGER :: i

        DO i=1, SIZE(self%channels)
            CALL self%channels(i)%add(it=it, point=self%point_cls)
        END DO

    END SUBROUTINE add


    !> Compares two objects p1 and p2 of class source_cls. In case these
    !! objects are equal `.TRUE.` is returned otherwise `.FALSE.`.
    !!
    !! Invokes the comparison operator '==' of point_mod::point_cls.
    !! In addition it uses the comparison operator '==' for objects of type
    !! channel_mod::channel_container.
    !!
    !! The function is used to overload the inherited comparison operator of
    !! class source_cls.
    !!
    !! @note Two objects p1 and p2 are called identical or equal if all their
    !!       components are equal:
    !!
    !!  * The base objects point_mod::point_cls are equal
    !!  * All channel_mod::channel_container objects are equal
    !!
    !! @param p1 Passed-object dummy of class source_cls
    !! @param p2 Passed-object dummy of class source_cls
    !! @return .TRUE. if p1 and p2 are identical otherwise .FALSE.
    ELEMENTAL LOGICAL FUNCTION equal_sources(p1, p2)
        ! Passed dummy arguments
        CLASS(source_cls), INTENT(IN) :: p1
        CLASS(point_cls), INTENT(IN) :: p2

        equal_sources = .FALSE.

        ! If p2 is not of class source_cls return
        ! Because p2 may hold any extension of point_cls
        IF ( .NOT. SAME_TYPE_AS(p1, p2) ) RETURN
        ! Makes p2's type explicit; It gives access to p2's components
        SELECT TYPE (p2)
            CLASS IS (source_cls) ! Guaranteed to execute since same type
                IF ( .NOT.(p1%point_cls==p2%point_cls) ) RETURN
                IF ( ALLOCATED(p1%channels) .NEQV. ALLOCATED(p2%channels) ) RETURN
                IF ( SIZE(p1%channels) /= SIZE(p2%channels) ) RETURN
                IF ( .NOT. ALL(p1%channels==p2%channels) ) RETURN
        END SELECT

        equal_sources = .TRUE.

    END FUNCTION


    !> Checks if an object of source_cls class is initialized or not
    !!
    !! Overloads the procedure inherited from point_mod::point_cls.
    !!
    !! @note A source - an object of class source_cls - is called
    !!       uninitialized if it was instantiated without passing any arguments
    !!       (i.e. `... = source_cls()`). For details see:
    !!       source_mod::uninitialized_point_cls().
    !!
    !! @param self Passed-object dummy of class source_cls
    !! @return .TRUE. if the passed object is uninitialized otherwise .FALSE.
    ELEMENTAL LOGICAL FUNCTION uninitialized(self)
        ! Passed dummy arguments
        CLASS(source_cls), INTENT(IN) :: self

        IF ( equal_sources(self, source_cls()) ) THEN
            uninitialized = .TRUE.
        ELSE
            uninitialized = .FALSE.
        END IF

    END FUNCTION


!------------------------------------------------------------------------------
! Public module procedures


    !> Subroutine parsing all occurring &source_sf and &source_mt
    !! NAMELIST groups in a file
    !!
    !! @param source_file_name Configuration file-name
    !! @param sources An allocatable array of type source_cls
    SUBROUTINE parse_source_file(source_file_name, sources)
        ! Use associated entities
        USE error_mod, ONLY : warn, abort
        USE parser_mod, ONLY : count_namelists, parse_namelist, &
                               source_mt_typ, source_sf_typ, source_sac_typ
        USE bcast_mod, ONLY : bcast_namelist, bcast_integer, reduce_sum
        ! Passed dummy variables
        CHARACTER(LEN=*), INTENT(IN) :: source_file_name
        TYPE(source_cls), INTENT(OUT), ALLOCATABLE :: sources(:)
        ! Local variables
        INTEGER :: source_file_unit, n_sf, n_mt, n_sac, i, n_srcs
        TYPE(source_sf_typ) :: source_sf
        TYPE(source_mt_typ) :: source_mt
        TYPE(source_sac_typ) :: source_sac

        ! Open receiver file
        IF ( is_mpi_master_rank() ) &
            OPEN( NEWUNIT=source_file_unit, FILE=source_file_name, &
                ACTION='READ', STATUS='OLD' )

        IF ( is_mpi_master_rank() ) THEN
            n_sf = count_namelists( unit=source_file_unit, nml_cls=source_sf )
            n_mt = count_namelists( unit=source_file_unit, nml_cls=source_mt )
            n_sac= count_namelists( unit=source_file_unit, nml_cls=source_sac)
        END IF

        ! Broadcast values
        CALL bcast_integer(i=n_sf, mr=mpi_master_rank)
        CALL bcast_integer(i=n_mt, mr=mpi_master_rank)
        CALL bcast_integer(i=n_sac,mr=mpi_master_rank)

        ! Allocate array of sources
        ALLOCATE( sources(n_sf+n_mt+n_sac) )

        IF ( is_mpi_master_rank() ) &
            REWIND( UNIT=source_file_unit )
        DO i=1, n_sf
            IF ( is_mpi_master_rank() ) THEN
                ! Parse NAMELIST
                CALL parse_namelist( unit=source_file_unit, nml_cls=source_sf )
            END IF
            ! Broadcast values
            CALL bcast_namelist(source_sf, mpi_master_rank)
            ! Invoke structure constructor
            sources(i) = source_cls(source_sf)
        END DO

        IF ( is_mpi_master_rank() ) &
            REWIND(UNIT=source_file_unit)
        DO i=1, n_mt
            IF ( is_mpi_master_rank() ) THEN
                ! Parse NAMELIST
                CALL parse_namelist( unit=source_file_unit, nml_cls=source_mt )
            END IF
            ! Broadcast values
            CALL bcast_namelist(source_mt, mpi_master_rank)
            ! Invoke structure constructor
            sources(i+n_sf) = source_cls(source_mt)
        END DO

        IF ( is_mpi_master_rank() ) &
            REWIND(UNIT=source_file_unit)
        DO i=1, n_sac
            IF ( is_mpi_master_rank() ) THEN
                ! Parse NAMELIST
                CALL parse_namelist( unit=source_file_unit, nml_cls=source_sac)
            END IF
            ! Broadcast values
            CALL bcast_namelist(source_sac, mpi_master_rank)
            ! Invoke structure constructor
            sources(i+n_sf+n_mt) = source_cls(source_sac)
        END DO


        ! All done -> close file
        IF ( is_mpi_master_rank() ) &
            CLOSE(source_file_unit)

        ! Sort out all source which are not located in this rank
        sources = PACK( sources(:), .NOT. sources(:)%uninitialized() )

        ! Check if there remains at least one source
        n_srcs = reduce_sum( i=SIZE(sources), mr=mpi_master_rank )
        IF ( n_srcs < 1 .AND. is_mpi_master_rank() ) &
            CALL abort( 'No sources defined in config file. &
                        &At least one source must be specified.' )

    END SUBROUTINE parse_source_file


!------------------------------------------------------------------------------
! Private module procedures


    !> Procedure which returns the source signature. Any other string than
    !! DELTA, RICKER, SIN_3, GAUSS, BRUNE is treated as ASCII file with
    !! floating point numbers.
    !!
    !! Built in Source wavelets are:
    !!  * DELTA: see source_signature_mod::delta()
    !!  * RICKER: see source_signature_mod::ricker()
    !!  * GAUSS: see source_signature_mod::gauss()
    !!  * SIN**3: see source_signature_mod::sin_3()
    !!  * BRUNE: see source_signature_mod::brune()
    !!  * Any other string is treated as filename: see read_stf()
    !!
    !! @todo 1st derivative of GAUSS
    !! @todo Built-in source wavelets: Double-check conditions for onset and
    !! width.
    !! @todo Make waring-messages more expressive.
    !!
    !! @param wavelet Either a filename or the name of a source wavelet (DELTA,
    !!                RICKER, GAUSS, SIN**3, BRUNE)
    !! @param onset Parameter for wavelet
    !! @param width Parameter for wavelet
    FUNCTION make_source_signature_(wavelet, onset, width) &
            RESULT(signature)
        ! Use associated entities
        USE error_mod, ONLY : warn, abort
        USE source_signature_mod, ONLY : delta, ricker, sin_3, gauss, brune
        ! Passed dummy variables
        CHARACTER(LEN=*), INTENT(IN) :: wavelet
        REAL(real_kind), INTENT(IN) :: onset, width
        ! Result variable
        REAL(real_kind) :: signature(config%nt())
        ! Local variables
        REAL(real_kind) :: t(config%nt()), dt
        INTEGER :: i, nt

        nt = config%nt()
        dt = config%dt()

        ! Make time vector
        t(:) = REAL([ (i, i=0, nt-1) ], real_kind)*dt

        ! Select the source signature
        SELECT CASE ( TRIM(wavelet) )

            ! Excite with delta-pulse
            CASE ( 'DELTA' )
                ! Check if parameter onset is appropriate
                IF ( .NOT. ANY(t(:)==onset) ) &
                    CALL warn("For a DELTA-pulse the value of 'onset' must pre&
                        &cisely coinside with a time-step (n*dt). This is not &
                        &the case! The entire exciting wavelet is set to zero.")
                ! Make source signature
                signature(:) = delta(time=t(:), onset=onset)

            ! Excite with Ricker-wavelet
            CASE ( 'RICKER' )
                ! Check if parameters onset and width are appropriate
                IF ( t(1) > onset .OR. t(nt) < onset ) &
                    CALL warn("For a RICKER-wavelet the value of 'onset' should&
                        & be in the range of in the simulation time.")
                IF ( width < 3*dt .OR. width > t(nt/2) ) &
                    CALL warn("The RICKER-wavelet's width appears to be either &
                        &too wide or too narrow.")
                ! Make source signature
                signature(:) = ricker( time=t(:), onset=onset, width=width)

            ! Excite with Gaussian pulse
            CASE ( 'GAUSS' )
                ! Check if parameters onset and width are appropriate
                IF ( t(1) > onset .OR. t(nt) < onset ) &
                    CALL warn("For a GAUSS-puls the value of 'onset' should&
                        & be in the range of in the simulation time.")
                IF ( width < 3*dt .OR. width > t(nt/2) ) &
                    CALL warn("The GAUSS-puls' width appears to be either &
                        &too wide or too narrow.")
                ! Make source signature
                signature(:) = gauss(time=t(:), onset=onset, width=width)

            ! Excite with sin()**3 wavelet
            CASE ( 'SIN**3' )
                ! Check if parameters onset and width are appropriate
                IF ( t(1) > onset .OR. t(nt) < onset ) &
                    CALL warn("For a SIN**3-puls the value of 'onset' should&
                        & be in the range of in the simulation time.")
                IF ( width < 3*dt .OR. width > t(nt/2) ) &
                    CALL warn("The SIN**3-puls' width appears to be either &
                        &too wide or too narrow.")
                ! Make source signature
                signature(:) = sin_3(time=t(:), onset=onset, width=width)

            ! Excite with Brune signature
            CASE ( 'BRUNE' )
                ! Check if parameters onset and width are appropriate
                IF ( t(1) > onset .OR. t(nt) < onset ) &
                    CALL warn("For a BRUNE-puls the value of 'onset' should&
                        & be in the range of in the simulation time.")
                IF ( width < 3*dt .OR. width > t(nt/2) ) &
                    CALL warn("The BRUNE-puls' width appears to be either &
                        &too wide or too narrow.")
                ! Make source signature
                signature(:) = brune(time=t(:), onset=onset, width=width)

            ! Any string which is not listed above will be interpreted
            ! as file-name
            CASE DEFAULT
                signature(:) = read_stf(wavelet, config%nt())

        END SELECT

    END FUNCTION


    !> Reads a series of reals line by line from an ASCII file
    !!
    !! @exception If something goes wrong reading the file
    !!
    !! @param fn File name
    !! @param nt Number of floating point numbers
    !! @result stf Source signature an array of nt reals
    FUNCTION read_stf(fn, nt) RESULT(stf)
        ! Use associated entities
        USE error_mod, ONLY : warn, abort
        USE string_utilities_mod, ONLY : i2c
        ! Passed dummy variables
        INTEGER, INTENT(IN) :: nt
        CHARACTER(LEN=*), INTENT(IN) :: fn
        ! Result variable
        REAL(real_kind) :: stf(nt)
        ! Local variables
        INTEGER :: lun, stat, i
        CHARACTER(fnl) :: msg

        OPEN( NEWUNIT=lun, FILE=TRIM(fn), ACTION='read', STATUS='old', &
              IOSTAT=stat, IOMSG=msg )
        IF ( stat /= 0 ) &
            CALL abort( TRIM(msg) )
        DO i=1, nt
            READ( UNIT=lun, FMT=*, IOSTAT=stat, IOMSG=msg ) stf(i)
            IF ( stat /= 0 ) &
                CALL abort("At line "//i2c(i)//" of File '"//TRIM(fn) &
                    // "' - " // TRIM(msg) )
        END DO
        CLOSE( UNIT=lun )

    END FUNCTION read_stf


END MODULE source_mod
