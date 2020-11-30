! Ses3d-NT - simulation of elastic wave propagation in spherical sections

! (c) by Stefan Mauerberger
!    and Maksym Melnyk
!
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


!> Module to handle the header of the SAC data format
!!
!! A SAC (Seismic Analysis Code) file consists of 158 words of header
!! accompanied by up to two data sections. This module contains the SAC
!! header definitions and procedures to handle SAC headers.
!!
!! For details see: http://www.iris.edu/software/sac/manual/file_format.html
!!
!! $Date: 2013-10-27 16:43:48 +0100 (Sun, 27 Oct 2013) $
!! $Author: mauerberger $
!! $Revision: 729 $
!! @copyright GNU General Public License version 3 or later
MODULE sac_header_mod
    USE parameters_mod, ONLY : rk => real_kind

    IMPLICIT NONE

    PRIVATE

    !> Character kind type
    INTEGER, PARAMETER :: ascii = SELECTED_CHAR_KIND("ascii")
    INTEGER, PARAMETER :: logical_kind = 4

    !> The SAC enumerated header field values (all lower case)
    !! @bug The following values are missing:
    !!  * IXYZ {General XYZ (3-D) file}
    !!  * IO (Other source of known origin)
    CHARACTER(LEN=8, KIND=ascii), PARAMETER, DIMENSION(*) :: &
        enumerated_header_variables = [ &
            'itime   ', & ! Time series file
            'irlim   ', & ! Spectral file---real and imaginary
            'iamph   ', & ! Spectral file---amplitude and phase
            'ixy     ', & ! General x versus y data
            'iunkn   ', & ! Unknown
            'idisp   ', & ! Displacement in nm
            'ivel    ', & ! Velocity in nm/sec
            'iacc    ', & ! Acceleration in nm/sec/sec
            'ib      ', & ! Begin time
            'iday    ', & ! Midnight of reference GMT day
            'io      ', & ! Event origin time
            'ia      ', & ! First arrival time
            'it0     ', & ! User defined time pick n, n=0,9
            'it1     ', &
            'it2     ', &
            'it3     ', &
            'it4     ', &
            'it5     ', &
            'it6     ', &
            'it7     ', &
            'it8     ', &
            'it9     ', &
            'iradnv  ', &
            'itannv  ', &
            'iradev  ', &
            'itanev  ', &
            'inorth  ', &
            'ieast   ', &
            'ihorza  ', &
            'idown   ', &
            'iup     ', &
            'illlbb  ', &
            'iwwsn1  ', &
            'iwwsn2  ', &
            'ihglp   ', &
            'isro    ', &
            'inucl   ', & ! Nuclear event
            'ipren   ', & ! Nuclear pre-shot event
            'ipostn  ', & ! Nuclear post-shot event
            'iquake  ', & ! Earthquake
            'ipreq   ', & ! Foreshock
            'ipostq  ', & ! Aftershock
            'ichem   ', & ! Chemical explosion
            'iother  ', & ! Other
            'igood   ', & ! Good data
            'iglch   ', & ! Glitches
            'idrop   ', & ! Dropouts
            'ilowsn  ', & ! Low signal to noise ratio
            'irldta  ', & ! Real data
            'ivolts  ', & ! Velocity in volts
            'imb     ', & ! Bodywave Magnitude
            'ims     ', & ! Surfacewave Magnitude
            'iml     ', & ! Local Magnitude
            'imw     ', & ! Moment Magnitude
            'imd     ', & ! Duration Magnitude
            'imx     ', & ! User Defined Magnitude
            'ineic   ', & ! National Earthquake Information Center
            'ipdeq   ', &
            'ipdew   ', &
            'ipde    ', & ! Preliminary Determination of Epicenter
            'iisc    ', & ! Internation Seismological Centre
            'ireb    ', & ! Reviewed Event Bulletin
            'iusgs   ', & ! US Geological Survey
            'ibrk    ', & ! UC Berkeley
            'icaltech', & ! California Institute of Technology
            'illnl   ', & ! Lawrence Livermore National Laboratory
            'ievloc  ', & ! Event Location (computer program
            'ijsop   ', & ! Joint Seismic Observation Program
            'iuser   ', & ! The individual using SAC2000
            'iunknown', & ! unknown
            'iqb     ', & ! Quarry or mine blast confirmed by quarry
            'iqb1    ', & ! Quarry/mine blast with designed shot info-ripple fired
            'iqb2    ', & ! Quarry/mine blast with observed shot info-ripple fired
            'iqbx    ', & ! Quarry or mine blast - single shot
            'iqmt    ', & ! Quarry/mining-induced events: tremors and rockbursts
            'ieq     ', & ! Earthquake
            'ieq1    ', & ! Earthquakes in a swarm or aftershock sequence
            'ieq2    ', & ! Felt earthquake
            'ime     ', & ! Marine explosion
            'iex     ', & ! Other explosion
            'inu     ', & ! Nuclear explosion
            'inc     ', & ! Nuclear cavity collapse
            'io_     ', &
            'il      ', & ! Local event of unknown origin
            'ir      ', & ! Regional event of unknown origin
            'it      ', & ! Teleseismic event of unknown origin
            'iu      ', & ! Undetermined or conflicting information
            'ieq3    ', &
            'ieq0    ', &
            'iex0    ', &
            'iqc     ', &
            'iqb0    ', &
            'igey    ', &
            'ilit    ', &
            'imet    ', &
            'iodor   ', &
            'ios     ' ]

    !> The 158 words SAC header field names (all upper case)
    CHARACTER(LEN=8, KIND=ascii), PARAMETER, DIMENSION(*) :: header_variables = &
        ! SAC-header floats - Words 1 to 70
        ['DELTA   ', 'DEPMIN  ', 'DEPMAX  ', 'SCALE   ', 'ODELTA  ', &
         'B       ', 'E       ', 'O       ', 'A       ', 'INTERNAL', &
         'T0      ', 'T1      ', 'T2      ', 'T3      ', 'T4      ', &
         'T5      ', 'T6      ', 'T7      ', 'T8      ', 'T9      ', &
         'F       ', 'RESP0   ', 'RESP1   ', 'RESP2   ', 'RESP3   ', &
         'RESP4   ', 'RESP5   ', 'RESP6   ', 'RESP7   ', 'RESP8   ', &
         'RESP9   ', 'STLA    ', 'STLO    ', 'STEL    ', 'STDP    ', &
         'EVLA    ', 'EVLO    ', 'EVEL    ', 'EVDP    ', 'MAG     ', &
         'USER0   ', 'USER1   ', 'USER2   ', 'USER3   ', 'USER4   ', &
         'USER5   ', 'USER6   ', 'USER7   ', 'USER8   ', 'USER9   ', &
         'DIST    ', 'AZ      ', 'BAZ     ', 'GCARC   ', 'INTERNAL', &
         'INTERNAL', 'DEPMEN  ', 'CMPAZ   ', 'CMPINC  ', 'XMINIMUM', &
         'XMAXIMUM', 'YMINIMUM', 'YMAXIMUM', 'UNUSED  ', 'UNUSED  ', &
         'UNUSED  ', 'UNUSED  ', 'UNUSED  ', 'UNUSED  ', 'UNUSED  ', &
        ! SAC-header integers - Words 71 to 105
         'NZYEAR  ', 'NZJDAY  ', 'NZHOUR  ', 'NZMIN   ', 'NZSEC   ', &
         'NZMSEC  ', 'NVHDR   ', 'NORID   ', 'NEVID   ', 'NPTS    ', &
         'INTERNAL', 'NWFID   ', 'NXSIZE  ', 'NYSIZE  ', 'UNUSED  ', &
         'IFTYPE  ', 'IDEP    ', 'IZTYPE  ', 'UNUSED  ', 'IINST   ', &
         'ISTREG  ', 'IEVREG  ', 'IEVTYP  ', 'IQUAL   ', 'ISYNTH  ', &
         'IMAGTYP ', 'IMAGSRC ', 'UNUSED  ', 'UNUSED  ', 'UNUSED  ', &
         'UNUSED  ', 'UNUSED  ', 'UNUSED  ', 'UNUSED  ', 'UNUSED  ', &
        ! SAC-header booleans - Words 106 to 110
         'LEVEN   ', 'LPSPOL  ', 'LOVROK  ', 'LCALDA  ', 'UNUSED  ', &
        ! SAC-header characters - Words 111 to 158
         'KSTNM   ','KSTNM_  ','KEVNM   ','KEVNM_  ','KEVNM__ ','KEVNM___',&
         'KHOLE   ','KHOLE_  ','KO      ','KO_     ','KA      ','KA_     ',&
         'KT0     ','KT0_    ','KT1     ','KT1_    ','KT2     ','KT2_    ',&
         'KT3     ','KT3_    ','KT4     ','KT4_    ','KT5     ','KT5_    ',&
         'KT6     ','KT6_    ','KT7     ','KT7_    ','KT8     ','KT8_    ',&
         'KT9     ','KT9_    ','KF      ','KF_     ','KUSER0  ','KUSER0_ ',&
         'KUSER1  ','KUSER1_ ','KUSER2  ','KUSER2_ ','KCMPNM  ','KCMPNM_ ',&
         'KNETWK  ','KNETWK_ ','KDATRD  ','KDATRD_ ','KINST   ','KINST_  ' ]

    !> Header class object
    !! @todo describe what it does
    TYPE :: header_cls
    PRIVATE
        REAL(rk) :: r
        INTEGER :: i
        LOGICAL :: l
        CHARACTER(LEN=16) :: c
    CONTAINS
    PRIVATE
        PROCEDURE, PASS(from) :: header_cls_to_real
        PROCEDURE, PASS(from) :: header_cls_to_integer
        PROCEDURE, PASS(from) :: header_cls_to_logical
        PROCEDURE, PASS(from) :: header_cls_to_character
        GENERIC, PUBLIC :: ASSIGNMENT(=) => header_cls_to_real, &
                                            header_cls_to_integer, &
                                            header_cls_to_logical, &
                                            header_cls_to_character
    END TYPE

    ! header_cls structure constructor
    INTERFACE header_cls
        MODULE PROCEDURE construct_from_real
        MODULE PROCEDURE construct_from_integer
        MODULE PROCEDURE construct_from_logical
        MODULE PROCEDURE construct_from_character
    END INTERFACE

    !> Either, for a passed string it returns the corresponding position/word in
    !! the SAC header otherwise it returns the SAC Header Variable name as
    !! string when an integer is passed
    INTERFACE header
        MODULE PROCEDURE header_int
        MODULE PROCEDURE header_char
    END INTERFACE

    !> Either, for a passed string it returns the corresponding integer value or
    !! it returns the SAC Enumerated Header Variable name as string when an
    !! integer is passed
    INTERFACE enum
        MODULE PROCEDURE enum_int
        MODULE PROCEDURE enum_char
    END INTERFACE

    PUBLIC :: ascii, rk, logical_kind, &
              header, header_variables, header_cls, &
              enum, enumerated_header_variables, &
              iday

CONTAINS

    !> Returns the position/word in the SAC header corresponding to the
    !! passed header keyword
    !!
    !! @exception If key is not a SAC header keyword -1 is returned
    !! @exception If key is UNUSED -2 is returned
    !! @exception If key is INTERNAL -3 is returned
    !!
    !> @param key SAC header keyword
    !> @return Position/word in the SAC header
    ELEMENTAL INTEGER FUNCTION header_int(key) RESULT(word)
        USE string_utilities_mod, ONLY : raise_case
        CHARACTER(LEN=*), INTENT(IN) :: key
        CHARACTER(LEN=LEN(key)) :: key_

        key_ = raise_case(ADJUSTL(key))

        IF ( key_ == 'UNUSED' ) THEN
            word = -2
            RETURN
        END IF

        IF ( key_ == 'INTERNAL' ) THEN
            word = -3
            RETURN
        END IF

        DO word=1, SIZE(header_variables)
            IF ( header_variables(word) == key_ ) &
                RETURN
        END DO

        word = -1

    END FUNCTION header_int

    !> Returns the SAC Header Variable keyword corresponding to passed word/
    !! position
    !!
    !! @exception If word is not in range [1-158] 'OUT OF BOUNDS' is returned
    !!
    !> @param word Position/word in the SAC header
    !> @return Corresponding SAC Header Variable keyword
    PURE FUNCTION header_char(word) RESULT(key)
        INTEGER, INTENT(IN) :: word
        CHARACTER(LEN=str_len_header(word)) :: key

        IF ( 1 <= word  .AND.  word <= SIZE(header_variables) ) THEN
            key = TRIM(header_variables(word))
        ELSE
            key = 'OUT OF BOUNDS'
        END IF

    END FUNCTION header_char

    !> Returns the length of a SAC Header Keyword at a position/word
    !! Trailing whitespaces are removed.
    !!
    !> @param word Position/word in the SAC header
    !> @return length of a SAC Header Keyword
    ELEMENTAL INTEGER FUNCTION str_len_header( word )
        INTEGER, INTENT(IN) :: word

        IF ( 1 <= word  .AND.  word <= SIZE( header_variables ) ) THEN
            str_len_header = LEN( TRIM( header_variables(word) ) )
        ELSE
            str_len_header = LEN( 'OUT OF BOUNDS' )
        END IF

    END FUNCTION str_len_header



    !> Returns the index/position of the passed SAC Enumerated Header Keyword
    !! (see sac_header_mod::enumerated_header_variables)
    !!
    !! @exception If keyword is not a SAC Enumerated Variables -1 is returned
    !!
    !> @param key SAC enumerator keyword
    !> @return Position/index in the SAC Enumerated Header Variables
    ELEMENTAL INTEGER FUNCTION enum_int(key) RESULT(word)
        USE string_utilities_mod, ONLY : lower_case
        CHARACTER(*), INTENT(IN) :: key
        CHARACTER(LEN=LEN(key)) :: key_

        key_ = lower_case(ADJUSTL(key))

        DO word=1, SIZE( enumerated_header_variables )
            IF ( enumerated_header_variables(word) == key_ ) &
                RETURN
        END DO

        word = -1

    END FUNCTION enum_int

    !> Returns the Enumerated Header Keyword at a certain position/index
    !!
    !! @exception If index/word is not in range 'OUT OF BOUNDS' is returned
    !!
    !! @param word Enumerated Header Index
    !! @return Enumerated Header Keyword with trailing whitespace removed
    PURE FUNCTION enum_char(word) RESULT(key)
        INTEGER, INTENT(IN) :: word
        CHARACTER(LEN=str_len_enum(word)) :: key

        IF ( 1 <= word  .AND.  word <= SIZE(enumerated_header_variables) ) THEN
            key = TRIM(enumerated_header_variables(word))
        ELSE
            key = 'OUT OF BOUNDS'
        END IF

    END FUNCTION enum_char

    !> Returns the length of the enumerate header keyword
    !! with trailing whitespace removed
    !!
    !! @exception If word is not in range 13 is returned
    !!
    !> @param word Passed position/index
    !> @return Keyword-length
    ELEMENTAL INTEGER FUNCTION str_len_enum(word)
        INTEGER, INTENT(IN) :: word

        IF ( 1 <= word .and. word <= SIZE( enumerated_header_variables ) ) THEN
            str_len_enum = LEN_TRIM( enumerated_header_variables(word) )
        ELSE
            str_len_enum = LEN( 'OUT OF BOUNDS' )
        END IF

    END FUNCTION str_len_enum



    !> Returns the day of the year number
    !!
    !! Examples:
    !!
    !!     iday( 2012,  1,  1) =>   1
    !!     iday( 2012, 10, 11) => 285
    !!
    !! @todo Check if passed values are valid
    !! @param year Year (e.g. 2000)
    !! @param month Month [1-12]
    !! @param day Day [1-31]
    !! @return Day of the year number [1-366]
    ELEMENTAL INTEGER FUNCTION iday( year, month, day )
        INTEGER, INTENT(IN) :: year, month, day

        iday = 3055 * ( month + 2 ) / 100 &
               - ( month + 10 ) / 13 * 2 - 91 + &
               ( 1 - ( MODULO( year, 4 ) + 3 ) / 4 &
               + ( MODULO( year, 100 ) + 99 ) / 100 - &
               ( MODULO( year, 400 ) +399 ) /400 ) * ( month + 10 ) / 13 + day

    END FUNCTION iday



    !> For the assignment interface it converts header_cls to real(rk)
    !!
    !! @param from A header_cls object
    !! @param to The corresponding real
    ELEMENTAL SUBROUTINE header_cls_to_real(to,from)
        REAL(rk), INTENT(OUT) :: to
        CLASS(header_cls), INTENT(IN) :: from
        to = from%r
    END SUBROUTINE

    TYPE(header_cls) ELEMENTAL FUNCTION construct_from_real(from) RESULT(new)
        REAL(rk), INTENT(IN) :: from
        INTEGER :: n
        n = STORAGE_SIZE(from)
        new = header_cls(r=from, &
                         i=TRANSFER(from,1),&
                         l=TRANSFER(from,.TRUE.), &
                         c=TRANSFER(from,REPEAT('c',n)) )
    END FUNCTION


    !> For the assignment interface it converts header_cls to integer
    !!
    !! @param from A header_cls object
    !! @param to The corresponding integer
    ELEMENTAL SUBROUTINE header_cls_to_integer(to,from)
        INTEGER, INTENT(out) :: to
        CLASS(header_cls), INTENT(in) :: from
        to = from%i
    END SUBROUTINE

    TYPE(header_cls) ELEMENTAL FUNCTION construct_from_integer(from) RESULT(new)
        INTEGER, INTENT(in) :: from
        INTEGER :: n
        n = STORAGE_SIZE(from)
        new = header_cls(r=TRANSFER(from,1.0_rk), &
                         i=from, &
                         l=TRANSFER(from,.TRUE.), &
                         c=TRANSFER(from,REPEAT('c',n)) )
    END FUNCTION


    !> For the assignment interface it converts header_cls to logical
    !!
    !! @param from A header_cls object
    !! @param to The corresponding logical
    ELEMENTAL SUBROUTINE header_cls_to_logical(to,from)
        LOGICAL, INTENT(out) :: to
        CLASS(header_cls), INTENT(in) :: from
        to = from%l
    END SUBROUTINE

    TYPE(header_cls) ELEMENTAL FUNCTION construct_from_logical(from) RESULT(new)
        LOGICAL, INTENT(in) :: from
        INTEGER :: n
        n = STORAGE_SIZE(from)
        new = header_cls(r=TRANSFER(from,1.0_rk), &
                         i=TRANSFER(from,1),&
                         l=from, &
                         c=TRANSFER(from,REPEAT('c',n)) )
    END FUNCTION


    !> For the assignment interface it converts header_cls to character(*)
    !!
    !! @param from A header_cls object
    !! @param to The corresponding character(*)
    ELEMENTAL SUBROUTINE header_cls_to_character(to,from)
        CHARACTER(LEN=*), INTENT(out) :: to
        CLASS(header_cls), INTENT(in) :: from
        to = from%c
    END SUBROUTINE

    TYPE(header_cls) ELEMENTAL FUNCTION construct_from_character(from) RESULT(new)
        CHARACTER(LEN=*), INTENT(in) :: from
        new = header_cls(r=TRANSFER(from,1.0_rk), &
                         i=TRANSFER(from,1),&
                         l=TRANSFER(from,.TRUE.), &
                         c=from )
    END FUNCTION


END MODULE sac_header_mod
