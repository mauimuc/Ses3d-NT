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
!! $Date: 2013-12-03 19:32:32 +0100 (Tue, 03 Dec 2013) $
!! $Author: mauerberger $
!! $Revision: 837 $
!! @copyright GNU General Public License version 3 or later


!> Parser module
!!
MODULE parser_mod
    USE error_mod, ONLY : abort
    USE parameters_mod, ONLY : max_attributes, fsl, fnl, real_kind

    IMPLICIT NONE

    PRIVATE

    ! In case a NAMELIST component is not given assume default values
    INTEGER, PARAMETER :: default_integer = HUGE(1)
    REAL(real_kind), PARAMETER :: default_real = HUGE(0.0_real_kind)
    CHARACTER(LEN=*), PARAMETER :: default_character = 'N/A'

    ! Abstract namelist class
    TYPE, ABSTRACT :: namelist_cls
    END TYPE

    ! Model
    TYPE, EXTENDS(namelist_cls) :: model_typ
        REAL(real_kind) :: lat_min, lat_max
        REAL(real_kind) :: lon_min, lon_max
        REAL(real_kind) :: rad_min, rad_max
        INTEGER :: model_type
        CHARACTER(LEN=fnl) :: rhoinv, lambda, mu, a, b, c, q
        LOGICAL :: override
    END TYPE

    ! General
    TYPE, EXTENDS(namelist_cls) :: general_typ
        CHARACTER(LEN=16) :: event_name, workflow
        CHARACTER(LEN=fnl) :: log_file_dir
    END TYPE

    ! Grid
    TYPE, EXTENDS(namelist_cls) :: grid_typ
        INTEGER :: nx, ny, nz, lpd, taper_width
        REAL(real_kind) :: taper_slope
    END TYPE

    ! Time
    TYPE, EXTENDS(namelist_cls) :: time_typ
        INTEGER :: nt, date_time(8)
        REAL(real_kind) :: dt
    END TYPE

    ! Receivers
    TYPE, EXTENDS(namelist_cls) :: receiver_typ
        REAL(real_kind) :: lat, lon, depth
        CHARACTER(LEN=8) :: network, station, location, &
                            attributes(max_attributes)
        LOGICAL :: override
        CHARACTER(LEN=fnl) :: prefix
    END TYPE

    ! Abstract Source Class
    TYPE, ABSTRACT, EXTENDS(namelist_cls) :: source_cls
        REAL(real_kind) :: lat, lon, depth, onset, width
        CHARACTER(LEN=fnl) :: wavelet
    END TYPE

    ! Monopole Sources
    TYPE, EXTENDS(source_cls) :: source_sf_typ
        REAL(real_kind) :: direction(3)
    END TYPE

    ! Dipole Source
    TYPE, EXTENDS(source_cls) :: source_mt_typ
        REAL(real_kind) :: moment_tensor(6)
    END TYPE

    ! SAC monopole source
    TYPE, EXTENDS(namelist_cls) :: source_sac_typ
        CHARACTER(LEN=fnl) :: file_name
    END TYPE

    ! Abstract Volume Class
    TYPE, ABSTRACT, EXTENDS(namelist_cls) :: volume_cls
        INTEGER :: timestep_start, timestep_increment, timestep_end
        CHARACTER(LEN=8) :: attributes(max_attributes)
    END TYPE

    ! Gradient
    TYPE, EXTENDS(volume_cls) :: grad_typ
    END TYPE

    !> Abstract derived type serving as base-type for parsing output related
    !! NAMELISTs. There is a subroutine bcast_mod::bcast_output_cls() which
    !! distributes values amongst all MPI ranks.
    TYPE, ABSTRACT, EXTENDS(volume_cls) :: output_cls
        LOGICAL :: override
        CHARACTER(LEN=fnl) :: prefix
    END TYPE

    !> Derived type used for parsing &input-NAMELISTs. The corresponding
    !! subroutine which actually parses the NAMELISTs is input_raw_parser().
    !!
    !! @todo Actually override is not needed
    TYPE, EXTENDS(output_cls) :: input_raw_typ
    END TYPE

    !> Derived type used for parsing &output_netcdf-NAMELISTs. The
    !! corresponding parser subroutine is output_netcdf_parser(). The
    !! subroutine for distributing values is bcast_mod::bcast_output_cls()
    TYPE, EXTENDS(output_cls) :: output_netcdf_typ
    END TYPE

    !> Derived type used for parsing &output_raw-NAMELISTs. The
    !! corresponding parser subroutine is output_netcdf_parser(). The
    !! subroutine for distributing values is bcast_mod::bcast_output_cls()
    TYPE, EXTENDS(output_cls) :: output_raw_typ
    END TYPE

    PUBLIC :: general_typ, model_typ, grad_typ, grid_typ, time_typ, &
              source_cls, source_sf_typ, source_mt_typ, source_sac_typ, &
              receiver_typ, &
              output_cls, output_raw_typ, output_netcdf_typ, &
              namelist_cls, count_namelists, parse_namelist, &
              default_integer, default_real, default_character

CONTAINS

    SUBROUTINE model_parser(unit, typ, iostat, iomsg)
        ! Passed dummy arguments
        INTEGER, INTENT(IN) :: unit
        TYPE(model_typ), INTENT(OUT) :: typ
        INTEGER, INTENT(OUT), OPTIONAL :: iostat
        CHARACTER(LEN=*), INTENT(OUT), OPTIONAL :: iomsg
        ! Local declarations initialized to default values
        INTEGER :: iostat_ = 0
        CHARACTER(LEN=fsl) :: iomsg_ = ''
        CHARACTER(LEN=*), PARAMETER :: err_suffix = &
            ' while parsing &model-group!'
        ! NAMELIST variable list
        REAL(real_kind) :: lat_max = 0.0, lat_min = 0.0, &
                           lon_max = 0.0, lon_min = 0.0, &
                           rad_max = 0.0, rad_min = 0.0
        INTEGER :: model_type = -1
        CHARACTER(LEN=fnl) :: rhoinv = 'N/A', lambda = 'N/A', mu = 'N/A', &
                              a = 'N/A', b = 'N/A', c = 'N/A', q = 'N/A'
        LOGICAL :: override = .FALSE.
        ! NAMELIST definition
        NAMELIST /model/ lat_max, lat_min, lon_max, lon_min, rad_max, rad_min,&
            model_type, rhoinv, lambda, mu, a, b, c, q, override

        IF ( PRESENT(iostat) .AND. PRESENT(iomsg) ) THEN
            READ(UNIT=unit, NML=model, IOSTAT=iostat_, IOMSG=iomsg_)
            iostat = iostat_; iomsg=TRIM(iomsg_)//err_suffix
        ELSE IF ( PRESENT(iostat) ) THEN
            READ(UNIT=unit, NML=model, IOSTAT=iostat_)
            iostat = iostat_
        ELSE
            READ(UNIT=unit, NML=model)
        END IF

        typ = model_typ(lat_min=lat_min, lat_max=lat_max, &
                        lon_min=lon_min, lon_max=lon_max, &
                        rad_min=rad_min, rad_max=rad_max, &
                        model_type=model_type, &
                        rhoinv=rhoinv, lambda=lambda, mu=mu, &
                        a=a, b=b, c=c, q=q, &
                        override=override )

    END SUBROUTINE model_parser

    SUBROUTINE grid_parser(unit, typ, iostat, iomsg)
        ! Passed dummy arguments
        INTEGER, INTENT(IN) :: unit
        TYPE(grid_typ), INTENT(OUT) :: typ
        INTEGER, INTENT(OUT), OPTIONAL :: iostat
        CHARACTER(LEN=*), INTENT(OUT), OPTIONAL :: iomsg
        ! Local declarations initialized to default values
        INTEGER :: iostat_ = 0
        CHARACTER(LEN=fsl) :: iomsg_ = ''
        CHARACTER(LEN=*), PARAMETER :: err_suffix = &
            ' while parsing &grid-group!'
        ! NAMELIST variable list
        INTEGER :: nx, ny, nz, lpd, taper_width
        REAL(real_kind) :: taper_slope
        ! NAMELIST definition
        NAMELIST /grid/ nx, ny, nz, lpd, taper_width, taper_slope
        ! Initialize local variables to default values
        nx = default_integer
        ny = default_integer
        nz = default_integer
        lpd = default_integer
        taper_width = default_integer
        taper_slope = default_integer

        IF ( PRESENT(iostat) .AND. PRESENT(iomsg) ) THEN
            READ(UNIT=unit, NML=grid, IOSTAT=iostat_, IOMSG=iomsg_)
            iostat = iostat_; iomsg=TRIM(iomsg_)//err_suffix
        ELSE IF ( PRESENT(iostat) ) THEN
            READ(UNIT=unit, NML=grid, IOSTAT=iostat_)
            iostat = iostat_
        ELSE
            READ(UNIT=unit, NML=grid)
        END IF

        typ = grid_typ(nx=nx, ny=ny, nz=nz, lpd=lpd, &
                       taper_width=taper_width, taper_slope=taper_slope)

    END SUBROUTINE grid_parser


    SUBROUTINE time_parser(unit, typ, iostat, iomsg)
        ! Passed dummy arguments
        INTEGER, INTENT(IN) :: unit
        TYPE(time_typ), INTENT(OUT) :: typ
        INTEGER, INTENT(OUT), OPTIONAL :: iostat
        CHARACTER(LEN=*), INTENT(OUT), OPTIONAL :: iomsg
        ! Local declarations initialized to default values
        INTEGER :: iostat_
        CHARACTER(LEN=fsl) :: iomsg_
        CHARACTER(LEN=*), PARAMETER :: err_suffix = &
            ' while parsing &time-group!'
        ! NAMELIST variable list
        INTEGER :: nt = -1, date_time(8) = -1
        REAL(real_kind) :: dt = -1.0
        ! NAMELIST definition
        NAMELIST /time/ dt, nt, date_time

        IF ( PRESENT(iostat) .AND. PRESENT(iomsg) ) THEN
            READ(UNIT=unit, NML=time, IOSTAT=iostat_, IOMSG=iomsg_)
            iostat = iostat_; iomsg=TRIM(iomsg_)//err_suffix
        ELSE IF ( PRESENT(iostat) ) THEN
            READ(UNIT=unit, NML=time, IOSTAT=iostat_)
            iostat = iostat_
        ELSE
            READ(UNIT=unit, NML=time)
        END IF

        typ = time_typ(dt=dt, nt=nt, date_time=date_time)

    END SUBROUTINE time_parser


    SUBROUTINE general_parser(unit, typ, iostat, iomsg)
        ! Passed dummy arguments
        INTEGER, INTENT(IN) :: unit
        TYPE(general_typ), INTENT(OUT) :: typ
        INTEGER, INTENT(OUT), OPTIONAL :: iostat
        CHARACTER(LEN=*), INTENT(OUT), OPTIONAL :: iomsg
        ! Local declarations initialized to default values
        INTEGER :: iostat_ = 0
        CHARACTER(LEN=fsl) :: iomsg_ = ''
        CHARACTER(LEN=*), PARAMETER :: err_suffix = &
            ' while parsing &general-group!'
        ! NAMELIST variable list
        CHARACTER(LEN=16) :: event_name = 'N/A', workflow = ''
        CHARACTER(LEN=fnl) :: log_file_dir = 'N/A'
        ! NAMELIST definition
        NAMELIST /general/ log_file_dir, event_name, workflow

        IF ( PRESENT(iostat) .AND. PRESENT(iomsg) ) THEN
            READ(UNIT=unit, NML=general, IOSTAT=iostat_, IOMSG=iomsg_)
            iostat = iostat_; iomsg=TRIM(iomsg_)//err_suffix
        ELSE IF ( PRESENT(iostat) ) THEN
            READ(UNIT=unit, NML=general, IOSTAT=iostat_)
            iostat = iostat_
        ELSE
            READ(UNIT=unit, NML=general)
        END IF

        typ = general_typ(event_name=event_name, &
                          workflow=workflow, &
                          log_file_dir=log_file_dir)

    END SUBROUTINE general_parser


    SUBROUTINE source_mt_parser(unit, typ, iostat, iomsg)
        ! Passed dummy arguments
        INTEGER, INTENT(IN) :: unit
        TYPE(source_mt_typ), INTENT(OUT) :: typ
        INTEGER, INTENT(OUT), OPTIONAL :: iostat
        CHARACTER(LEN=*), INTENT(OUT), OPTIONAL :: iomsg
        ! Local declarations initialized to default values
        INTEGER :: iostat_ = 0
        CHARACTER(LEN=fsl) :: iomsg_ = ''
        CHARACTER(LEN=*), PARAMETER :: err_suffix = &
            ' while parsing &source_mt-group!'
        ! NAMELIST variable list
        REAL(real_kind) :: lon, lat, depth, onset, width
        REAL(real_kind) :: moment_tensor(6)
        CHARACTER(LEN=fnl) :: wavelet
        ! NAMELIST definition
        NAMELIST /source_mt/ lat, lon, depth, moment_tensor, &
            wavelet, onset, width
        ! Initialize local variables to default values
        lon = default_real
        lat = default_real
        depth = default_real
        onset = default_real
        width = default_real
        moment_tensor(6) = default_real
        wavelet = default_character

        IF ( PRESENT(iostat) .AND. PRESENT(iomsg) ) THEN
            READ(UNIT=unit, NML=source_mt, IOSTAT=iostat_, IOMSG=iomsg_)
            iostat = iostat_; iomsg=TRIM(iomsg_)//err_suffix
        ELSE IF ( PRESENT(iostat) ) THEN
            READ(UNIT=unit, NML=source_mt, IOSTAT=iostat_)
            iostat = iostat_
        ELSE
            READ(UNIT=unit, NML=source_mt)
        END IF

        typ = source_mt_typ(lat=lat, lon=lon, depth=depth, onset=onset, &
                            width=width, moment_tensor=moment_tensor(:), &
                            wavelet=wavelet)

    END SUBROUTINE source_mt_parser


    SUBROUTINE source_sf_parser(unit, typ, iostat, iomsg)
        ! Passed dummy arguments
        INTEGER, INTENT(IN) :: unit
        TYPE(source_sf_typ), INTENT(OUT) :: typ
        INTEGER, INTENT(OUT), OPTIONAL :: iostat
        CHARACTER(LEN=*), INTENT(OUT), OPTIONAL :: iomsg
        ! Local declarations initialized to default values
        INTEGER :: iostat_ = 0
        CHARACTER(LEN=fsl) :: iomsg_ = ''
        CHARACTER(LEN=*), PARAMETER :: err_suffix = &
            ' while parsing &source_sf-group!'
        ! NAMELIST variable list
        REAL(real_kind) :: lon, lat, depth, onset, width
        REAL(real_kind) :: direction(3)
        CHARACTER(LEN=fnl) :: wavelet
        ! NAMELIST definition
        NAMELIST /source_sf/ lat, lon, depth, direction, wavelet, onset, width
        ! Initialize local variables to default values
        lon = default_real
        lat = default_real
        depth = default_real
        onset = default_real
        width = default_real
        direction(3) = default_real
        wavelet = default_character

        IF ( PRESENT(iostat) .AND. PRESENT(iomsg) ) THEN
            READ(UNIT=unit, NML=source_sf, IOSTAT=iostat_, IOMSG=iomsg_)
            iostat = iostat_; iomsg=TRIM(iomsg_)//err_suffix
        ELSE IF ( PRESENT(iostat) ) THEN
            READ(UNIT=unit, NML=source_sf, IOSTAT=iostat_)
            iostat = iostat_
        ELSE
            READ(UNIT=unit, NML=source_sf)
        END IF

        typ = source_sf_typ(lat=lat, lon=lon, depth=depth, onset=onset, &
                            width=width, direction=direction(:), &
                            wavelet=wavelet)

    END SUBROUTINE source_sf_parser

    SUBROUTINE source_sac_parser(unit, typ, iostat, iomsg)
        ! Passed dummy arguments
        INTEGER, INTENT(IN) :: unit
        TYPE(source_sac_typ), INTENT(OUT) :: typ
        INTEGER, INTENT(OUT), OPTIONAL :: iostat
        CHARACTER(LEN=*), INTENT(OUT), OPTIONAL :: iomsg
        ! Local declarations initialized to default values
        INTEGER :: iostat_ = 0
        CHARACTER(LEN=fsl) :: iomsg_ = ''
        CHARACTER(LEN=*), PARAMETER :: err_suffix = &
            ' while parsing &source_sac-group!'
        ! NAMELIST variable list
        CHARACTER(LEN=fnl) :: file_name
        ! NAMELIST definition
        NAMELIST /source_sac/ file_name
        ! Initialize local variables to default values
        file_name = default_character

        IF ( PRESENT(iostat) .AND. PRESENT(iomsg) ) THEN
            READ(UNIT=unit, NML=source_sac, IOSTAT=iostat_, IOMSG=iomsg_)
            iostat = iostat_; iomsg=TRIM(iomsg_)//err_suffix
        ELSE IF ( PRESENT(iostat) ) THEN
            READ(UNIT=unit, NML=source_sac, IOSTAT=iostat_)
            iostat = iostat_
        ELSE
            READ(UNIT=unit, NML=source_sac)
        END IF

        typ = source_sac_typ(file_name=file_name)

    END SUBROUTINE source_sac_parser

    SUBROUTINE receiver_parser(unit, typ, iostat, iomsg)
        ! Passed dummy arguments
        INTEGER, INTENT(IN) :: unit
        TYPE(receiver_typ), INTENT(OUT) :: typ
        INTEGER, INTENT(OUT), OPTIONAL :: iostat
        CHARACTER(LEN=*), INTENT(OUT), OPTIONAL :: iomsg
        ! Local declarations initialized to default values
        INTEGER :: iostat_ = 0
        CHARACTER(LEN=fsl) :: iomsg_ = ''
        CHARACTER(LEN=*), PARAMETER :: err_suffix = &
            ' while parsing &receiver-group!'
        ! NAMELIST variable list
        REAL(real_kind) :: lon, lat, depth
        CHARACTER(LEN=8) :: network, station, location
        CHARACTER(LEN=8) :: attributes(max_attributes)
        LOGICAL :: override
        CHARACTER(LEN=fnl) :: prefix
        ! NAMELIST definition
        NAMELIST /receiver/ lat, lon, depth, network, station, location, &
                            attributes, override, prefix
        ! Initialize local variables to default values
        lon = default_real
        lat = default_real
        depth = default_real
        network = 'Ses3d-NT'
        station = ''
        location = ''
        attributes(:) = ''
        override = .FALSE.
        prefix = './'

        IF ( PRESENT(iostat) .AND. PRESENT(iomsg) ) THEN
            READ(UNIT=unit, NML=receiver, IOSTAT=iostat_, IOMSG=iomsg_)
            iostat = iostat_; iomsg=TRIM(iomsg_)//err_suffix
        ELSE IF ( PRESENT(iostat) ) THEN
            READ(UNIT=unit, NML=receiver, IOSTAT=iostat_)
            iostat = iostat_
        ELSE
            READ(UNIT=unit, NML=receiver)
        END IF

        typ = receiver_typ(lat=lat, lon=lon, depth=depth, &
                           network=network, location=location, &
                           station=station, override=override, &
                           attributes=attributes(:), prefix=prefix)

    END SUBROUTINE receiver_parser


    SUBROUTINE grad_parser(unit, typ, iostat, iomsg)
        ! Passed dummy arguments
        INTEGER, INTENT(IN) :: unit
        TYPE(grad_typ), INTENT(OUT) :: typ
        INTEGER, INTENT(OUT), OPTIONAL :: iostat
        CHARACTER(LEN=*), INTENT(OUT), OPTIONAL :: iomsg
        ! Local declarations initialized to default values
        INTEGER :: iostat_ = 0
        CHARACTER(LEN=fsl) :: iomsg_ = ''
        CHARACTER(LEN=*), PARAMETER :: err_suffix = &
            ' while parsing &grad-group!'
        ! NAMELIST variable list
        INTEGER :: timestep_start = -1, timestep_increment = -1, timestep_end = -1
        CHARACTER(LEN=8) :: attributes(max_attributes) = ''
        ! NAMELIST definition
        NAMELIST /grad/ timestep_start, timestep_increment, timestep_end, &
            attributes
        ! XXX weird attributes is not initialized to ''
        attributes(:) = SPREAD('', 1, max_attributes)

        IF ( PRESENT(iostat) .AND. PRESENT(iomsg) ) THEN
            READ(UNIT=unit, NML=grad, IOSTAT=iostat_, IOMSG=iomsg_)
            iostat = iostat_; iomsg=TRIM(iomsg_)//err_suffix
        ELSE IF ( PRESENT(iostat) ) THEN
            READ(UNIT=unit, NML=grad, IOSTAT=iostat_)
            iostat = iostat_
        ELSE
            READ(UNIT=unit, NML=grad)
        END IF

        typ = grad_typ(timestep_start=timestep_start, &
                       timestep_increment=timestep_increment, &
                       timestep_end=timestep_end, &
                       attributes=attributes(:))

    END SUBROUTINE grad_parser


    SUBROUTINE output_netcdf_parser(unit, typ, iostat, iomsg)
        ! Passed dummy arguments
        INTEGER, INTENT(IN) :: unit
        TYPE(output_netcdf_typ), INTENT(OUT) :: typ
        INTEGER, INTENT(OUT), OPTIONAL :: iostat
        CHARACTER(LEN=*), INTENT(OUT), OPTIONAL :: iomsg
        ! Local declarations
        INTEGER :: iostat_ = 0
        CHARACTER(LEN=fsl) :: iomsg_ = ''
        CHARACTER(LEN=*), PARAMETER :: err_suffix = &
            ' while parsing &output_netcdf-group!'
        ! NAMELIST variable list
        INTEGER :: timestep_start, timestep_increment, timestep_end
        LOGICAL :: override
        CHARACTER(LEN=fnl) :: prefix
        CHARACTER(LEN=8) :: attributes(max_attributes) = ''
        ! NAMELIST definition
        NAMELIST /output_netcdf/ timestep_start, timestep_increment, &
            timestep_end, override, prefix, attributes
        ! Initialize local variables to default values
        timestep_start = default_integer
        timestep_increment = default_integer
        timestep_end = default_integer
        override = .FALSE.
        prefix = default_character
        attributes(:) = ''

        IF ( PRESENT(iostat) .AND. PRESENT(iomsg) ) THEN
            READ(UNIT=unit, NML=output_netcdf, IOSTAT=iostat_, IOMSG=iomsg_)
            iostat = iostat_; iomsg=TRIM(iomsg_)//err_suffix
        ELSE IF ( PRESENT(iostat) ) THEN
            READ(UNIT=unit, NML=output_netcdf, IOSTAT=iostat_)
            iostat = iostat_
        ELSE
            READ(UNIT=unit, NML=output_netcdf)
        END IF

        typ = output_netcdf_typ(timestep_start=timestep_start, &
                                timestep_increment=timestep_increment, &
                                timestep_end=timestep_end, &
                                override=override, &
                                prefix=prefix, &
                                attributes=attributes(:))

    END SUBROUTINE output_netcdf_parser


    SUBROUTINE output_raw_parser(unit, typ, iostat, iomsg)
        ! Passed dummy arguments
        INTEGER, INTENT(IN) :: unit
        TYPE(output_raw_typ), INTENT(OUT) :: typ
        INTEGER, INTENT(OUT), OPTIONAL :: iostat
        CHARACTER(LEN=*), INTENT(OUT), OPTIONAL :: iomsg
        ! Local declarations initialized to default values
        INTEGER :: iostat_ = 0
        CHARACTER(LEN=fsl) :: iomsg_ = ''
        CHARACTER(LEN=*), PARAMETER :: err_suffix = &
            ' while parsing &output_raw-group!'
        ! NAMELIST variable list
        INTEGER :: timestep_start, timestep_increment, timestep_end
        LOGICAL :: override
        CHARACTER(LEN=fnl) :: prefix
        CHARACTER(LEN=8) :: attributes(max_attributes)
        ! NAMELIST definition
        NAMELIST /output_raw/ timestep_start, timestep_increment, timestep_end, &
            attributes, override, prefix
        ! Initialize local variables to default values
        timestep_start = default_integer
        timestep_increment = default_integer
        timestep_end = default_integer
        override = .FALSE.
        prefix = default_character
        attributes(:) = ''

        IF ( PRESENT(iostat) .AND. PRESENT(iomsg) ) THEN
            READ(UNIT=unit, NML=output_raw, IOSTAT=iostat_, IOMSG=iomsg_)
            iostat = iostat_; iomsg=TRIM(iomsg_)//err_suffix
        ELSE IF ( PRESENT(iostat) ) THEN
            READ(UNIT=unit, NML=output_raw, IOSTAT=iostat_)
            iostat = iostat_
        ELSE
            READ(UNIT=unit, NML=output_raw)
        END IF

        typ = output_raw_typ(timestep_start=timestep_start, &
                             timestep_increment=timestep_increment, &
                             timestep_end=timestep_end, &
                             override=override, &
                             prefix=prefix, &
                             attributes=attributes(:))

    END SUBROUTINE output_raw_parser


    SUBROUTINE input_raw_parser(unit, typ, iostat, iomsg)
        ! Passed dummy arguments
        INTEGER, INTENT(IN) :: unit
        TYPE(input_raw_typ), INTENT(OUT) :: typ
        INTEGER, INTENT(OUT), OPTIONAL :: iostat
        CHARACTER(LEN=*), INTENT(OUT), OPTIONAL :: iomsg
        ! Local declarations initialized to default values
        INTEGER :: iostat_ = 0
        CHARACTER(LEN=fsl) :: iomsg_ = ''
        CHARACTER(LEN=*), PARAMETER :: err_suffix = &
            ' while parsing &input_raw-group!'
        ! NAMELIST variable list
        INTEGER :: timestep_start, timestep_increment, timestep_end
        CHARACTER(LEN=fnl) :: prefix
        CHARACTER(LEN=8) :: attributes(max_attributes)
        ! NAMELIST definition
        NAMELIST /input_raw/ timestep_start, timestep_increment, timestep_end, &
            attributes, prefix
        ! Initialize local variables to default values
        timestep_start = default_integer
        timestep_increment = default_integer
        timestep_end = default_integer
        prefix = default_character
        attributes(:) = ''

        IF ( PRESENT(iostat) .AND. PRESENT(iomsg) ) THEN
            READ(UNIT=unit, NML=input_raw, IOSTAT=iostat_, IOMSG=iomsg_)
            iostat = iostat_; iomsg=TRIM(iomsg_)//err_suffix
        ELSE IF ( PRESENT(iostat) ) THEN
            READ(UNIT=unit, NML=input_raw, IOSTAT=iostat_)
            iostat = iostat_
        ELSE
            READ(UNIT=unit, NML=input_raw)
        END IF

        typ = input_raw_typ(timestep_start=timestep_start, &
                            timestep_increment=timestep_increment, &
                            timestep_end=timestep_end, &
                            override=.TRUE., &
                            prefix=prefix, &
                            attributes=attributes(:))

    END SUBROUTINE input_raw_parser



    INTEGER FUNCTION count_namelists( unit, nml_cls ) RESULT(n)
        USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY : IOSTAT_END
        INTEGER, INTENT(IN) :: unit
        CLASS(namelist_cls), INTENT(OUT) :: nml_cls
        INTEGER :: iostat_ = 0
        CHARACTER(fsl) :: iomsg_ = ''

        n = 0
        REWIND( UNIT=unit )

        DO
            CALL parse_namelist(unit, nml_cls, iostat=iostat_, iomsg=iomsg_)
            SELECT CASE (iostat_)
                CASE (IOSTAT_END)
                    EXIT
                CASE (0)
                    n = n + 1
                CASE DEFAULT
                    CALL abort( "Soemthing went in function 'count_namelists' &
                               &- " // TRIM(iomsg_) )

            END SELECT

        END DO

    END FUNCTION count_namelists


    SUBROUTINE parse_namelist( unit, nml_cls, iostat, iomsg )
        INTEGER, INTENT(IN) :: unit
        CLASS(namelist_cls), INTENT(OUT) :: nml_cls
        INTEGER, INTENT(OUT), OPTIONAL :: iostat
        CHARACTER(LEN=*), INTENT(OUT), OPTIONAL :: iomsg

        SELECT TYPE (nml_cls)
            TYPE IS (model_typ)
                CALL model_parser(unit=unit, typ=nml_cls, &
                    iostat=iostat, iomsg=iomsg)
            TYPE IS (grid_typ)
                CALL grid_parser(unit=unit, typ=nml_cls, &
                    iostat=iostat, iomsg=iomsg)
            TYPE IS (time_typ)
                CALL time_parser(unit=unit, typ=nml_cls, &
                    iostat=iostat, iomsg=iomsg)
            TYPE IS (general_typ)
                CALL general_parser(unit=unit, typ=nml_cls, &
                    iostat=iostat, iomsg=iomsg)

            TYPE IS (source_sf_typ)
                CALL source_sf_parser(unit=unit, typ=nml_cls, &
                    iostat=iostat, iomsg=iomsg)
            TYPE IS (source_mt_typ)
                CALL source_mt_parser(unit=unit, typ=nml_cls, &
                    iostat=iostat, iomsg=iomsg)
            TYPE IS (source_sac_typ)
                CALL source_sac_parser(unit=unit, typ=nml_cls, &
                    iostat=iostat, iomsg=iomsg)

            TYPE IS (receiver_typ)
                CALL receiver_parser(unit=unit, typ=nml_cls, &
                    iostat=iostat, iomsg=iomsg)

            TYPE IS (grad_typ)
                CALL grad_parser(unit=unit, typ=nml_cls, &
                    iostat=iostat, iomsg=iomsg)

            TYPE IS (input_raw_typ)
                CALL input_raw_parser(unit=unit, typ=nml_cls, &
                    iostat=iostat, iomsg=iomsg)

            TYPE IS (output_netcdf_typ)
                CALL output_netcdf_parser(unit=unit, typ=nml_cls, &
                    iostat=iostat, iomsg=iomsg)
            TYPE IS (output_raw_typ)
                CALL output_raw_parser(unit=unit, typ=nml_cls, &
                    iostat=iostat, iomsg=iomsg)

            CLASS DEFAULT
                CALL abort( "A unknown class was passed to subroutine &
                    &'parse_namelist'" )

        END SELECT

    END SUBROUTINE parse_namelist


END MODULE parser_mod

