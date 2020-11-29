! Ses3d-NT - simulation of elastic wave propagation in spherical sections

! (c) by Stefan Mauerberger <mauerberger@geophysik.uni-muenchen.de>
!    and Maksym Melnyk <mmelnyk@geophysik.uni-muenchen.de>

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



!> Module for reading, writing and modifying SAC trace files
!!
!! For details see: http://www.iris.edu/software/sac/manual/file_format.html
!!
!! @todo Documentation
!!
!! $Date: 2013-11-27 13:59:40 +0100 (Wed, 27 Nov 2013) $
!! $Author: mauerberger $
!! $Revision: 822 $
!! @copyright GNU General Public License version 3 or later
MODULE sac_io_mod
    USE sac_header_mod, ONLY : rk, ascii, logical_kind, iday, &
                               header, header_cls, enum

    IMPLICIT NONE

    PRIVATE

    INTEGER :: i

    !> Class for reading, writing and modifying SAC traces
    TYPE :: sac_cls
    PRIVATE
        ! SAC header
        REAL(rk) :: floats(1:70)  = -12345.0
        INTEGER :: ints(71:105) = [(-12345,i=1,6), 6,(-12345,i=1,28)]
        LOGICAL :: bools(106:110) = [.TRUE., .FALSE., .FALSE., .FALSE., .FALSE.]
        CHARACTER(LEN=4) :: chars(111:158) = [(['-123','45..'],i=1,24)]
        ! First Data Section
        ! Dependent variable, amplitude or REAL component
        REAL(rk), ALLOCATABLE :: data1(:)
        ! Second Data Section (if present)
        ! Independent variable unevenly spaced, phase or imaginary component
        REAL(rk), ALLOCATABLE :: data2(:)
    CONTAINS
    PRIVATE
        PROCEDURE :: allocate => allocate_data_arrays
        !> Returns "[network].[station].[location].[channel].sac"
        !! (see sac_io_mod::auto_file_name())
        PROCEDURE, PUBLIC :: auto_file_name
        !> Writes trace to 'auto_file_name'
        PROCEDURE, PRIVATE :: write_sac_auto
        !> Writes trace to 'file'
        PROCEDURE, PRIVATE :: write_sac
        !> Writes SAC trace; file, iostat and iomsg are optional
        GENERIC, PUBLIC :: write => write_sac, write_sac_auto

        !> Write the SAC header to stdout (see sac_mod::print_sac_header)
        PROCEDURE :: print_sac_header
        !> Writes the SAC header into a unit identifier
        !! (see sac_mod::print_sac_header_unit)
        PROCEDURE :: print_sac_header_unit
        !> Generic interface for writing the SAC header.
        !! Optional a unit identifier may be passed.
        GENERIC, PUBLIC :: print => print_sac_header, print_sac_header_unit

        PROCEDURE :: get_header_key
        PROCEDURE :: get_header_word
        GENERIC, PUBLIC :: get_header => get_header_key, get_header_word

        PROCEDURE :: set_hdr_int_flt
        PROCEDURE :: set_hdr_int_int
        PROCEDURE :: set_hdr_int_bool
        PROCEDURE :: set_hdr_int_char
        PROCEDURE :: set_hdr_char_flt
        PROCEDURE :: set_hdr_char_int
        PROCEDURE :: set_hdr_char_bool
        PROCEDURE :: set_hdr_char_char
        GENERIC, PUBLIC :: set_header => set_hdr_int_flt,  set_hdr_char_flt, &
                                         set_hdr_int_int,  set_hdr_char_int, &
                                         set_hdr_int_bool, set_hdr_char_bool, &
                                         set_hdr_int_char, set_hdr_char_char

        !> Set the dependent data section (see sac_io_mod::set_data1_)
        PROCEDURE :: set_data1_
        !> Set a dependent data point at certain time-step
        !! see sac_io_mod::set_data1_at)
        PROCEDURE :: set_data1_at
        !> Generic interface for either setting the entire dependent data section
        !! or a single point at a certain time-step
        GENERIC, PUBLIC :: set_data1 => set_data1_, set_data1_at
        !> Set the independent data section (see sac_io_mod::set_data2_)
        PROCEDURE :: set_data2_
        !> Set an independent data point at certain time-step
        !! (see sac_io_mod::set_data2_at)
        PROCEDURE :: set_data2_at
        !> Generic interface for either setting the entire independent data
        !! section or a single point at a certain time-step
        GENERIC, PUBLIC :: set_data2 => set_data2_, set_data2_at
        !> Set dependent and independent data sections
        !! (see sac_io_mod::set_data1_data2_)
        PROCEDURE :: set_data1_data2_
        !> Set a dependent and an independent data point at certain time-step
        !! (see sac_io_mod::set_data1_data2_at)
        PROCEDURE :: set_data1_data2_at
        !> Generic interface for either setting the entire dependent and
        !! independent data section or a single point at a certain time-step
        GENERIC, PUBLIC :: set_data1_data2 => set_data1_data2_, set_data1_data2_at

        !> Returns the whole dependent data section (see sac_mod::get_data1_)
        PROCEDURE :: get_data1_
        !> Returns a point of dependent data at a certain time-step
        !! (see sac_mod::get_data1_at)
        PROCEDURE :: get_data1_at
        !> Generic interface for either returning the whole dependent data or
        !! at at certain time-step
        GENERIC, PUBLIC :: get_data1 => get_data1_, get_data1_at
        !> Returns the whole independent data section (see sac_mod::get_data2_)
        PROCEDURE :: get_data2_
        !> Returns a point of independent data at a certain time-step
        !! (see sac_mod::get_data2_at)
        PROCEDURE :: get_data2_at
        !> Generic interface for either returning the whole dependent data or
        !! at at certain time-step
        GENERIC, PUBLIC :: get_data2 => get_data2_, get_data2_at
        !> Returns both, the whole dependent and independent data sections
        !! (see sac_mod::get_data1_data2_)
        PROCEDURE :: get_data1_data2_
        !> Returns a point of dependent and independent data at a certain
        !! time-step
        !! (see sac_mod::get_data1_data2_at)
        PROCEDURE :: get_data1_data2_at
        !> Generic interface for either returning the whole dependent and
        !! independent data or at at certain time-step
        GENERIC, PUBLIC :: get_data1_data2 => get_data1_data2_, get_data1_data2_at
    END TYPE

    !> Structure constructor
    INTERFACE sac_cls
        !> See sac_mod::read_sac
        MODULE PROCEDURE read_sac
        !> See sac_mod::construct_w_req
        MODULE PROCEDURE construct_w_req
        !> See sac_mod::construct_trace
        MODULE PROCEDURE construct_trace
    END INTERFACE sac_cls

    PUBLIC :: sac_cls

CONTAINS


    !> Function to construct a SAC trace
    !> @param lat Latitude [rad]
    !> @param lon Longitude [rad]
    !> @param depth Depth [m]
    !> @param data Time series data
    !> @param dt Time increment
    !> @param network Network name
    !> @param station Station name
    !> @param locatoin Location
    !> @param channel Channel name
    !> @param event Event name
    !> @param override If true existing files are overridden
    !> @param date_time Date and time [yyy, mm, dd, dm_GMT, h, m, s, ms]
    !> @param prefix
    !> @return A SAC trace object
    FUNCTION construct_trace(lat, lon, depth, data, dt, &
                             network, station, location, channel, &
                             event, override, date_time) RESULT(tr)
        TYPE(sac_cls) :: tr
        REAL(rk), INTENT(IN) :: lat, lon, depth, data(:), dt
        CHARACTER(LEN=*), INTENT(IN) :: network, station, location, channel, event
        LOGICAL, INTENT(IN) :: override
        INTEGER, INTENT(IN) :: date_time(8)
        REAL(rk) :: T
        INTEGER :: iday_

        T = dt * REAL(SIZE(data), rk)
        iday_ = iday( date_time(1), date_time(2), date_time(3) )

        CALL tr%set_header(key='STLA', value=lat) ! Station latitude (degrees, north positive)
        CALL tr%set_header(key='STLO', value=lon) ! Station longitude (degrees, east positive)
        CALL tr%set_header(key='STDP', value=depth) ! Station depth below surface (meters)

        CALL tr%set_header(key='KEVNM' , value=event) ! Event name
        CALL tr%set_header(key='KNETWK', value=network) ! Name of seismic network
        CALL tr%set_header(key='KSTNM', value=station) ! Station name
        CALL tr%set_header(key='KHOLE', value=location) ! location identifier
        CALL tr%set_header(key='KCMPNM', value=channel) ! Channel name

        CALL tr%set_header(key='LOVROK', value=override) ! TRUE if it is okay to overwrite this file on disk

        CALL tr%set_header(key='NZYEAR', value=date_time(1)) ! year
        CALL tr%set_header(key='NZJDAY', value=iday_ )
        CALL tr%set_header(key='NZHOUR', value=date_time(5)) ! hour
        CALL tr%set_header(key='NZMIN', value=date_time(6)) ! minutes
        CALL tr%set_header(key='NZSEC', value=date_time(7)) ! seconds
        CALL tr%set_header(key='NZMSEC', value=date_time(8)) ! milliseconds

        CALL tr%set_header(key='DELTA', value=dt) ! Increment between evenly spaced samples
        CALL tr%set_header(key='NPTS', value=SIZE(data)) ! Number of points per data component
        CALL tr%set_data1(data)

        CALL tr%set_header( key='B'     , value=0.0_rk ) ! Beginning value of the independent variable
        CALL tr%set_header( key='E'     , value=T ) ! Ending value of the independent variable

        CALL tr%set_header( key='IFTYPE', value=1 ) ! Type of file
        CALL tr%set_header( key='ISYNTH', value=2 ) ! Synthetic data flag
        CALL tr%set_header( key='SCALE' , value=1.0_rk ) ! Multiplying scale factor for dependent variable

    END FUNCTION construct_trace


    !> Structure constructor for reading SAC files
    !! @warning The LOGICAL-KIND-TYPE is not well defined
    !! @param file_name Path and filename as string
    !! @return Reads a SAC-file and returns a sac_cls object
    FUNCTION read_sac(file_name) RESULT(self)
        USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY : REAL32, INT32
        TYPE(sac_cls) :: self
        CHARACTER(*), INTENT(IN) :: file_name
        INTEGER :: lun
        REAL(REAL32) :: floats_(1:70)
        INTEGER(INT32) :: ints_(71:105), npts_, iftype_
        LOGICAL(4) :: bools_(106:110), leven_
        CHARACTER(LEN=4,KIND=ascii) :: chars_(111:158)
        REAL(REAL32), ALLOCATABLE :: data1_(:), data2_(:)

        ! Check boolean storage size
        IF ( STORAGE_SIZE(bools_) /= 32 ) &
            STOP 'Logical sotrage size /= 32'

        OPEN(NEWUNIT=lun, FILE=file_name, ACCESS='stream', ACTION='read', &
             FORM='unformatted', STATUS='old' )

        ! Read file header, the first 158 words
        ! REAL32, INTEGER32, LOGICAL XXX, ASCII
        READ(UNIT=lun) floats_
        READ(UNIT=lun) ints_
        READ(UNIT=lun) bools_
        READ(UNIT=lun) chars_

        ! Number of sampling points
        npts_   = ints_( header('NPTS') )
        ! See if data are evenly spaced
        leven_  = bools_( header('LEVEN') )
        ! See which type it is e.g. spectral data
        iftype_ = ints_( header('IFTYPE') )

        ! Allocate memory to for the dependent variable
        ALLOCATE( data1_(npts_) )
        ! Read First Data Section
        ! REAL32
        READ( UNIT=lun ) data1_

        ! Handle independent variable
        IF ( .NOT. leven_  &
             .OR.  iftype_==enum('IRLIM') &      ! IRLIM {Spectral file---real and imaginary}
             .OR.  iftype_==enum('IAMPH') ) THEN ! IAMPH {Spectral file---amplitude and phase}
            ALLOCATE( data2_(npts_) )
            ! Read Second Data Section
            ! REAL32
            READ( UNIT=lun ) data2_
            !self%data2(:)  = data2_
        ELSE
            ! Not needed; just allocate an empty array
            ALLOCATE( data2_(0) )
        END IF

        ! All don; close unit identifier
        CLOSE( UNIT=lun )

        ! Invoke the default structure constructor
        self = sac_cls(floats=floats_, ints=ints_, bools=bools_, &
            chars=chars_, data1=data1_, data2=data2_)

    END FUNCTION read_sac



    !> Structure constructor demanding for all definitely required field
    !! values not including data sections
    ELEMENTAL FUNCTION construct_w_req( NPTS, B, E, IFTYPE, LEVEN, DELTA ) &
            RESULT( new_sac_typ )
        INTEGER, INTENT(IN) :: NPTS, IFTYPE
        REAL(rk), INTENT(IN) :: B, E, DELTA
        LOGICAL, INTENT(IN) :: LEVEN
        TYPE( sac_cls ) :: new_sac_typ

        CALL new_sac_typ%set_header( key='NPTS'  , value=NPTS )
        CALL new_sac_typ%set_header( key='B'     , value=B )
        CALL new_sac_typ%set_header( key='E'     , value=E )
        CALL new_sac_typ%set_header( key='IFTYPE', value=IFTYPE )
        CALL new_sac_typ%set_header( key='LEVEN' , value=LEVEN )
        CALL new_sac_typ%set_header( key='DELTA' , value=DELTA )

    END FUNCTION construct_w_req


!=====================
! Member procedures
!=====================

    !> Returns true if the SAC trace has an independent data section
    !!
    !! @param self A sac_cls object (usually the passed-object dummy)
    !! @return false if there is NO independent data data section
    ELEMENTAL LOGICAL FUNCTION has_independent_data(self)
        CLASS(sac_cls), INTENT(IN) :: self

        LOGICAL :: leven
        INTEGER :: iftype

        iftype = self%ints( header( 'IFTYPE' ) )
        leven = self%bools( header( 'LEVEN' ) )

        IF ( .NOT. leven  & ! Evenly space
             .OR.  iftype==enum('IRLIM') &       ! IRLIM {Spectral file---real and imaginary}
             .OR.  iftype==enum('IAMPH') ) THEN  ! IAMPH {Spectral file---amplitude and phase}
            has_independent_data = .TRUE.
        ELSE
            has_independent_data = .FALSE.
        END IF

    END FUNCTION


    PURE FUNCTION get_header_key(self, key) RESULT(val)
        CLASS(sac_cls), INTENT(IN) :: self
        CHARACTER(LEN=*), INTENT(IN) :: key
        TYPE(header_cls) :: val

        val = get_header_word(self, header(key))

    END FUNCTION

    ELEMENTAL FUNCTION get_header_word(self, word) RESULT(val)
        CLASS(sac_cls), INTENT(IN) :: self
        INTEGER, INTENT(IN) :: word
        TYPE(header_cls) :: val

        SELECT CASE (word)
            CASE(1:70)
                val = header_cls(self%floats(word))
            CASE(71:105)
                val = header_cls(self%ints(word))
            CASE(106:110)
                val = header_cls(self%bools(word))
            CASE(111,117,119,121,123,125,127,129,131,133,135,137,139,141,143,&
                 145,147,149,151,153,155,157)
                val = header_cls(self%chars(word) // self%chars(word+1))
            CASE(113)
                val = header_cls(self%chars(word) // self%chars(word+1) // &
                                 self%chars(word+2) // self%chars(word+3) )
            !CASE DEFAULT
            !    val = header(0)
        END SELECT

    END FUNCTION


    !> Allocation data sections
    ELEMENTAL SUBROUTINE allocate_data_arrays( self )
        CLASS(sac_cls), INTENT(INOUT) :: self
        INTEGER :: npts

        npts = self%ints( header( 'NPTS' ) )

        ! First Data Section
        IF ( ALLOCATED( self%data1 ) ) &
            DEALLOCATE( self%data1 )
        ALLOCATE( self%data1( npts ), SOURCE=-1234.5_rk )

        ! Second Data Section
        IF ( ALLOCATED( self%data2 ) ) &
            DEALLOCATE( self%data2 )
        IF ( has_independent_data(self) ) THEN
            ALLOCATE( self%data2( npts ), SOURCE=-1234.5_rk )
        ELSE
            ALLOCATE( self%data2(0) )
        END IF

    END SUBROUTINE allocate_data_arrays


    !> Returns a string which is usually used auto file name
    !!
    !! Naming convention is [network].[station].[location].[channel].sac
    !!
    !! @todo To be consistent with e.g. ObsPy that procedure should better be
    !!       name 'id'
    !!
    !! @param self A sac_cls object (usually the passed-object dummy)
    !! @return a string of length auto_file_name_lenght() which is set to
    !!         [network].[station].[location].[channel].sac
    PURE FUNCTION auto_file_name(self) RESULT(fn)
        CLASS(sac_cls), INTENT(IN) :: self
        INTEGER, PARAMETER :: id(8) = [153,154,117,118,111,112,151,152]
        CHARACTER(*), PARAMETER :: fmt_ = "(4(A,A,'.'),'sac')"
        CHARACTER(LEN=auto_file_name_length(self)) :: fn
        INTEGER :: i

        WRITE(fn,fmt_) ( TRIM(self%chars(id(i))), i=1, SIZE(id) )

    END FUNCTION auto_file_name

    !> Returns the length of the auto_file_name
    !!
    !! @param A sac_cls object (usually the passed-object dummy)
    !! @return The length of the auto_file_name string
    ELEMENTAL INTEGER FUNCTION auto_file_name_length(self) RESULT(len)
        TYPE(sac_cls), INTENT(IN) :: self
        INTEGER, PARAMETER :: id(8) = [153,154,117,118,111,112,151,152]

        len = SUM(LEN_TRIM(self%chars(id))) + 9

    END FUNCTION


    !> Writes the SAC-trace into a file. The filename is obtained by using the
    !! default naming convention (see sac_io_mod::auto_file_name()).
    !!
    !! @param A sac_cls object (usually the passed-object dummy)
    !! @param I/O status (optional)
    !! @param I/O message (optional)
    SUBROUTINE write_sac_auto(self, iostat, iomsg)
        CLASS(sac_cls), INTENT(IN) :: self
        CHARACTER(*), INTENT(OUT), OPTIONAL :: iomsg
        INTEGER, INTENT(OUT), OPTIONAL :: iostat

        CALL write_sac(self=self, file='./'//auto_file_name(self), &
                       iostat=iostat, iomsg=iomsg)

    END SUBROUTINE

    !> Writes the SAC-trace to file.
    !!
    !! @param A sac_cls object (usually the passed-object dummy)
    !! @param file Filename to write the SAC-trace into
    !! @param I/O status (optional)
    !! @param I/O message (optional)
    SUBROUTINE write_sac(self, file, iostat, iomsg)
        USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY : REAL32, INT32
        CLASS(sac_cls), INTENT(IN) :: self
        CHARACTER(*), INTENT(IN) :: file
        CHARACTER(*), INTENT(OUT), OPTIONAL :: iomsg
        INTEGER, INTENT(OUT), OPTIONAL :: iostat
        CHARACTER(7) :: status
        INTEGER :: lun

        ! Determine file status at open
        IF ( self%bools(header('LOVROK')) ) THEN
            status = 'replace' ! Override existing files
        ELSE
            status = 'new' ! DO NOT override
        END IF

        ! Open file depending on optional arguments iostat and iomsg
        IF ( PRESENT(iostat) .AND. PRESENT(iomsg) ) THEN
            OPEN( NEWUNIT=lun, FILE=file, ACCESS='stream', FORM='unformatted', &
                  ACTION='write', STATUS=status, IOSTAT=iostat, IOMSG=iomsg )
            IF ( iostat /= 0 ) RETURN
        ELSE IF( PRESENT(iostat) ) THEN
            OPEN( NEWUNIT=lun, FILE=file, ACCESS='stream', FORM='unformatted', &
                  ACTION='write', STATUS=status, IOSTAT=iostat )
            IF ( iostat /= 0 ) RETURN
        ELSE
            OPEN( NEWUNIT=lun, FILE=file, ACCESS='stream', FORM='unformatted', &
                  ACTION='write', STATUS=status )
        END IF

        ! Write header
        WRITE( UNIT=lun ) REAL(    self%floats, KIND=REAL32 )
        WRITE( UNIT=lun ) INT(     self%ints  , KIND=INT32 )
        WRITE( UNIT=lun ) LOGICAL( self%bools , KIND=logical_kind )
        WRITE( UNIT=lun ) self%chars
        ! Write data
        IF ( ALLOCATED( self%data1 ) ) &
            WRITE( UNIT=lun ) REAL( self%data1, KIND=REAL32 )
        IF ( ALLOCATED( self%data2 ) ) &
            WRITE( UNIT=lun ) REAL( self%data2, KIND=REAL32 )

        ! Close file
        CLOSE( lun )

        IF( PRESENT(iostat) ) &
            iostat = 0
        IF( PRESENT(iomsg) ) &
            iomsg = "File '" // TRIM(file) // "' successfully written"

    END SUBROUTINE write_sac


    !> Writes the SAC trace header to stdout
    !!
    !! @param A sac_cls object (usually the passed-object dummy)
    SUBROUTINE print_sac_header(self)
        USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY : OUTPUT_UNIT
        CLASS(sac_cls), INTENT(IN) :: self

        CALL self%print_sac_header_unit(unit=OUTPUT_UNIT)

    END SUBROUTINE

    !> Writes the SAC trace header into a unit specifier
    !! @param self A SAC object
    !! @param unit A unit specifier to write into
    SUBROUTINE print_sac_header_unit(self, unit)
        CLASS(sac_cls), INTENT(IN) :: self
        INTEGER, INTENT(IN) :: unit

        WRITE(UNIT=unit, FMT='(5G15.7)') self%floats

        WRITE(UNIT=unit, FMT='(5I10)') self%ints

        WRITE(UNIT=unit, FMT='(5L10)') self%bools

        WRITE(UNIT=unit, FMT='( 2(A4), " ", 4(A4) )') self%chars(111:116)

        WRITE(UNIT=unit, FMT='( 3(2A4, " ") )') self%chars(117:)

    END SUBROUTINE


    !> Sets the float values of the corresponding character keywords to the SAC header
    ELEMENTAL SUBROUTINE set_hdr_char_flt( self, key, value )
        CLASS( sac_cls ), INTENT(INOUT) :: self
        CHARACTER(*), INTENT(IN) :: key
        REAL(rk), INTENT(IN) :: value

        CALL self%set_hdr_int_flt( header(key), value )

    END SUBROUTINE set_hdr_char_flt


    !> Sets the float values of the corresponding integer words to the SAC header
    ELEMENTAL SUBROUTINE set_hdr_int_flt( self, word, value )
        CLASS( sac_cls ), INTENT(INOUT) :: self
        INTEGER, INTENT(IN) :: word
        REAL(rk), INTENT(IN) :: value

        self%floats(word) = value

    END SUBROUTINE set_hdr_int_flt


    !> Sets the integer values of the corresponding character keywords to the SAC header
    ELEMENTAL SUBROUTINE set_hdr_char_int( self, key, value )
        CLASS( sac_cls ), INTENT(INOUT) :: self
        CHARACTER(*), INTENT(IN) :: key
        INTEGER, INTENT(IN) :: value

        CALL self%set_hdr_int_int( header(key), value )

    END SUBROUTINE set_hdr_char_int


    !> Sets the integer values of the corresponding integer words to the SAC header
    ELEMENTAL SUBROUTINE set_hdr_int_int( self, word, value )
        CLASS( sac_cls ), INTENT(INOUT) :: self
        INTEGER, INTENT(IN) :: word
        INTEGER, INTENT(IN) :: value

        self%ints(word) = value

    END SUBROUTINE set_hdr_int_int


    !> Sets the logical values of the corresponding character keywords to the SAC header
    ELEMENTAL SUBROUTINE set_hdr_char_bool( self, key, value )
        CLASS( sac_cls ), INTENT(INOUT) :: self
        CHARACTER(*), INTENT(IN) :: key
        LOGICAL, INTENT(IN) :: value

        CALL self%set_hdr_int_bool( header(key), value )

    END SUBROUTINE set_hdr_char_bool


    !> Sets the logical values of the corresponding integer words to the SAC header
    ELEMENTAL SUBROUTINE set_hdr_int_bool( self, word, value )
        CLASS( sac_cls ), INTENT(INOUT) :: self
        INTEGER, INTENT(IN) :: word
        LOGICAL, INTENT(IN) :: value

        self%bools(word) = value

    END SUBROUTINE set_hdr_int_bool


    !> Sets the character values of the corresponding character keywords to the SAC header
    ELEMENTAL SUBROUTINE set_hdr_char_char( self, key, value )
        CLASS( sac_cls ), INTENT(INOUT) :: self
        CHARACTER(*), INTENT(IN) :: key, value

        CALL self%set_hdr_int_char( header(key), value )

    END SUBROUTINE set_hdr_char_char


    !> Sets the character values of the corresponding integer words to the SAC header
    ELEMENTAL SUBROUTINE set_hdr_int_char( self, word, value )
        CLASS( sac_cls ), INTENT(INOUT) :: self
        INTEGER, INTENT(IN) :: word
        CHARACTER(*), INTENT(IN) :: value
        CHARACTER(16) :: val_

        val_ = '                ' ! Ensures val_ to be of white spaces
        val_ = value ! Everything longer than 16 chars gets truncated

        SELECT CASE(word)
            CASE DEFAULT
                self%chars(word)   = val_(1:4)
                self%chars(word+1) = val_(5:8)
            CASE(113)
                self%chars(word)   = val_(1:4)
                self%chars(word+1) = val_(5:8)
                self%chars(word+2) = val_(9:12)
                self%chars(word+3) = val_(13:16)
        END SELECT

    END SUBROUTINE set_hdr_int_char



    !> Set dependent data (auto re-allocation)
    !! @param self A sac_cls object (usually the passed-object dummy)
    !! @param data1 Dependent data section (array of floats)
    PURE SUBROUTINE set_data1_(self, data)
        CLASS(sac_cls), INTENT(INOUT) :: self
        REAL(rk), INTENT(IN) :: data(:)

        self%data1 = data

    END SUBROUTINE

    !> Set dependent data-point at a certain position
    !! @param self A sac_cls object (usually the passed-object dummy)
    !! @param it Position where to write data (commonly this is the time-step)
    !! @param data Dependent data (single float)
    ELEMENTAL SUBROUTINE set_data1_at(self, it, data)
        CLASS(sac_cls), INTENT(INOUT) :: self
        INTEGER, INTENT(IN) :: it
        REAL(rk), INTENT(IN) :: data

        IF ( .NOT. ALLOCATED(self%data1) ) &
            CALL self%allocate()
        self%data1(it) = data

    END SUBROUTINE



    !> Set independent data (auto re-allocation)
    !! @param self A sac_cls object (usually the passed-object dummy)
    !! @param data2 Independent data section (array of floats)
    PURE SUBROUTINE set_data2_(self, data)
        CLASS(sac_cls), INTENT(INOUT) :: self
        REAL(rk), INTENT(IN) :: data(:)

        self%data2 = data

    END SUBROUTINE

    !> Set independent data-point at a certain position
    !! @param self A sac_cls object (usually the passed-object dummy)
    !! @param it Position where to write data (commonly this is the time-step)
    !! @param data Independent data (single float)
    ELEMENTAL SUBROUTINE set_data2_at(self, it, data)
        CLASS(sac_cls), INTENT(INOUT) :: self
        INTEGER, INTENT(IN) :: it
        REAL(rk), INTENT(IN) :: data

        IF ( .NOT. ALLOCATED(self%data2) ) &
            CALL self%allocate()
        self%data2(it) = data

    END SUBROUTINE



    !> Set dependent and independent data (auto re-allocation)
    !! @param self A sac_cls object (usually the passed-object dummy)
    !! @param data1 Dependent data section (array of floats)
    !! @param data2 Independent data section (array of floats)
    PURE SUBROUTINE set_data1_data2_(self, data1, data2)
        CLASS(sac_cls), INTENT(INOUT) :: self
        REAL(rk), INTENT(IN) :: data1(:), data2(:)

        CALL set_data1_(self=self, data=data1)
        CALL set_data2_(self=self, data=data2)

    END SUBROUTINE

    !> Set dependent and independent data points at a certain position
    !! @param self A sac_cls object (usually the passed-object dummy)
    !! @param it Position where to write data (commonly this is the time-step)
    !! @param data1 Dependent data (single float)
    !! @param data2 Independent data (single float)
    ELEMENTAL SUBROUTINE set_data1_data2_at(self, it, data1, data2)
        CLASS(sac_cls), INTENT(INOUT) :: self
        INTEGER, INTENT(IN) :: it
        REAL(rk), INTENT(IN) :: data1, data2

        CALL set_data1_at(self=self, it=it, data=data1)
        CALL set_data2_at(self=self, it=it, data=data2)

    END SUBROUTINE



    !> Returns the whole dependent data section
    !! @param self A sac_cls object (usually the passed-object dummy)
    !! @return Array of NPTS floats [data1]
    PURE FUNCTION get_data1_(self)
        CLASS(sac_cls), INTENT(IN) :: self
        REAL(rk) :: get_data1_(SIZE(self%data1))

        get_data1_(:) = self%data1

    END FUNCTION

    !> Returns a point of dependent data at a certain time-step
    !! @param self A sac_cls object (usually the passed-object dummy)
    !! @param it Position where to write data (commonly this is the time-step)
    !! @return Float data1 at time-step it
    REAL(rk) ELEMENTAL FUNCTION get_data1_at(self, it)
        CLASS(sac_cls), INTENT(IN) :: self
        INTEGER, INTENT(IN) :: it

        get_data1_at = self%data1(it)

    END FUNCTION



    !> Returns the whole independent data section
    !! @param self A sac_cls object (usually the passed-object dummy)
    !! @return Array of NPTS floats [data2]
    PURE FUNCTION get_data2_( self )
        CLASS(sac_cls), INTENT(IN) :: self
        REAL(rk) :: get_data2_(SIZE(self%data2))

        get_data2_(:) = self%data2

    END FUNCTION

    !> Returns a point of dependent data at a certain time-step
    !! @param self A sac_cls object (usually the passed-object dummy)
    !! @param it Position where to write data (commonly this is the time-step)
    !! @return Float data2 at time-step it
    REAL(rk) ELEMENTAL FUNCTION  get_data2_at(self, it)
        CLASS(sac_cls), INTENT(IN) :: self
        INTEGER, INTENT(IN) :: it

        get_data2_at = self%data2(it)

    END FUNCTION



    !> Returns both, the whole dependent and independent data sections
    !! @param self A sac_cls object (usually the passed-object dummy)
    !! @return Array of 2xNPTS floats [1,data1], [2,data2]
    PURE FUNCTION get_data1_data2_( self )
        CLASS(sac_cls), INTENT(IN) :: self
        REAL(rk) :: get_data1_data2_(2,SIZE(self%data2))

        get_data1_data2_(1,:) = self%data1
        get_data1_data2_(2,:) = self%data2

    END FUNCTION

    !> Returns a point of dependent and independent data at a certain time-step
    !! @param self A sac_cls object (usually the passed-object dummy)
    !! @param it Position where to write data (commonly this is the time-step)
    !! @return Array of two floats [data1,data2] at time-step it
    PURE FUNCTION get_data1_data2_at(self, it)
        CLASS(sac_cls), INTENT(IN) :: self
        INTEGER, INTENT(IN) :: it
        REAL(rk) :: get_data1_data2_at(2)

        get_data1_data2_at(1) = get_data2_at(self=self, it=it)
        get_data1_data2_at(2) = get_data2_at(self=self, it=it)

    END FUNCTION


END MODULE sac_io_mod
