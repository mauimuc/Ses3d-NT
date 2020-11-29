!> @file
!! Ses3d-NT - simulation of elastic wave propagation in spherical sections
!!
!! (c) by Michael Haas <mhaas@geophysik.uni-muenchen.de>
!!        Stefan Mauerberger <mauerberger@geophysik.uni-muenchen.de>
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
!! $Date: 2013-11-25 20:20:46 +0100 (Mon, 25 Nov 2013) $
!! $Author: mauerberger $
!! $Revision: 814 $
!! @copyright GNU General Public License version 3 or later



!> Module to write output of a Ses3d-NT simulation (6D array)
!! into a 3d array in netCDF format
MODULE netcdf_parallel_io_mod
    USE configuration_mod, ONLY : config => configuration
    USE parameters_mod
#ifdef NetCDF
    USE netcdf, ONLY : NF90_NOERR
#endif

    IMPLICIT NONE

    PRIVATE

#ifdef NetCDF
    PUBLIC :: open_netcdf_file, my_nf90_end_def, my_nf90_close, &
              my_nf90_def_var, my_nf90_put_var, add_metadata, nf90_noerr

CONTAINS


    ! -------------------------------------------------------------------------
    ! Opens a netCDF file
    ! -------------------------------------------------------------------------
    SUBROUTINE open_netCDF_file( file_name, ncid, dimid, time, &
                                 nf90_err, override )
        USE netcdf, ONLY : NF90_CLOSE, NF90_CREATE, NF90_ENDDEF, &
                           NF90_HDF5, NF90_MPIIO, NF90_NOCLOBBER, &
                           NF90_NOERR
#ifndef MPI_INCLUDE
        !! RECOMMENDED !!
        ! in case MPI_INCLUDE is NOT defined the mpi module will be used
        USE mpi, ONLY : MPI_INFO_NULL
#else
        !! DEPRECATED !!
        ! when MPI_INCLUDE is defined, mpi header-files will be included
        INCLUDE 'mpif.h'
#endif
        ! File ID, Dimension ID, Coordinate ID
        INTEGER, INTENT(OUT) :: ncid, dimid(4)
        INTEGER :: coordid(4), nf90_err_, cmode
        CHARACTER(LEN=fnl) :: file_name
        REAL, INTENT(IN) :: time
        LOGICAL, INTENT(IN), OPTIONAL :: override
        INTEGER, INTENT(OUT), OPTIONAL :: nf90_err

        cmode = NF90_HDF5+NF90_MPIIO+NF90_NOCLOBBER
        IF ( PRESENT(override) ) THEN
            IF ( override .eqv. .TRUE. ) &
            cmode = NF90_HDF5+NF90_MPIIO
        END IF

        nf90_err_ = NF90_CREATE( path=file_name, &
                                 cmode=cmode, &
                                 ncid=ncid, &
                                 comm=config%mpi_comm_cart(), &
                                 info=MPI_INFO_NULL )

        IF ( PRESENT(nf90_err) ) THEN
            IF ( nf90_err_ /= NF90_NOERR ) THEN
                nf90_err = nf90_err_
                RETURN
            ELSE
                nf90_err = nf90_err_
            END IF
        ELSE
            CALL check( nf90_err_ )
        END IF

        CALL create_netCDF_volume_dims( ncid, dimid, coordid(:), time )

    END SUBROUTINE open_netCDF_file


    ! -------------------------------------------------------------------------
    ! Creating dimensions of a netCDF file
    ! and writing current time-step as meta-data
    ! -------------------------------------------------------------------------
    SUBROUTINE create_netCDF_volume_dims( ncid, dimid, coordid, time )
        USE coordinate_utilities_mod, ONLY : rad2deg, colat2lat
        USE netcdf, ONLY: NF90_DEF_DIM, NF90_DEF_VAR, NF90_PUT_VAR, NF90_PUT_ATT
        INTEGER, INTENT(IN)  :: ncid
        INTEGER, INTENT(OUT) :: coordid(4), dimid(4)
        INTEGER :: nf90_err, date_time(8)
        REAL(real_kind) :: x(config%nx_glb()*config%lpd() + 1), &
                           y(config%ny_glb()*config%lpd() + 1), &
                           z(config%nz_glb()*config%lpd() + 1)
        CHARACTER(LEN=fsl) :: fmt
        CHARACTER(LEN=40) :: time_string
        real, intent(in) :: time

        ! Fetch coordinates
        x(:) = config%x_coord_glb_r1()
        y(:) = config%y_coord_glb_r1()
        z(:) = config%z_coord_glb_r1()
        x(:) = colat2lat( rad2deg( x ) )
        y(:) = rad2deg( y(:) )

        ! Define the dimensions. NetCDF will hand back an ID for each.
        nf90_err = NF90_DEF_DIM( ncid, "latitude",  SIZE(x), dimid(1) )
        CALL check( nf90_err )
        nf90_err = NF90_DEF_DIM( ncid, "longitude", SIZE(y), dimid(2) )
        CALL check( nf90_err )
        nf90_err = NF90_DEF_DIM( ncid, "r",         SIZE(z), dimid(3) )
        CALL check( nf90_err )
        nf90_err = NF90_DEF_DIM( ncid, "time",      1,       dimid(4) )
        CALL check( nf90_err )

        ! Define the coordinate variables.
        nf90_err = NF90_DEF_VAR( ncid, "latitude",  my_nf90_real, &
                                 dimid(1), coordid(1) )
        CALL check( nf90_err )
        nf90_err = NF90_DEF_VAR( ncid, "longitude", my_nf90_real, &
                                 dimid(2), coordid(2) )
        CALL check( nf90_err )
        nf90_err = NF90_DEF_VAR( ncid, "r",         my_nf90_real, &
                                 dimid(3), coordid(3))
        CALL check( nf90_err )
        nf90_err = NF90_DEF_VAR( ncid, "time",      my_nf90_real, &
                                 dimid(4), coordid(4))
        CALL check( nf90_err )

        ! Assign units attributes to coordinate var data.
        nf90_err = NF90_PUT_ATT( ncid, coordid(1), 'units', 'degrees_north' )
        CALL check( nf90_err )
        nf90_err = NF90_PUT_ATT( ncid, coordid(2), 'units', 'degrees_east' )
        CALL check( nf90_err )
        nf90_err = NF90_PUT_ATT( ncid, coordid(3), 'units', 'meters' )
        CALL check( nf90_err )
        ! Generate udunits date and time string
        date_time(:) = config%date_time()
        fmt = '( "seconds since ", I0,"-",I0,"-",I0, &
              & 1x , I0":",I0,":",I0,".",I0 )'
        WRITE(UNIT=time_string, FMT=fmt) date_time([1,2,3,5,6,7,8])
        ! Assign start time
        nf90_err = NF90_PUT_ATT( ncid, coordid(4), 'units', time_string )
        CALL check( nf90_err )

        ! Write the coordinate variable data
        nf90_err = NF90_PUT_VAR( ncid, coordid(1), REAL(x(:),my_nf90_real_kind) )
        CALL check( nf90_err )
        nf90_err = NF90_PUT_VAR( ncid, coordid(2), REAL(y(:),my_nf90_real_kind) )
        CALL check( nf90_err )
        nf90_err = NF90_PUT_VAR( ncid, coordid(3), REAL(z(:),my_nf90_real_kind) )
        CALL check( nf90_err )
        nf90_err = NF90_PUT_VAR( ncid, coordid(4), time )
        CALL check( nf90_err )

    END SUBROUTINE create_netCDF_volume_dims


    ! -------------------------------------------------------------------------
    !
    ! -------------------------------------------------------------------------
    SUBROUTINE add_metadata( ncid, it )
        USE netcdf, ONLY: NF90_PUT_ATT, NF90_GLOBAL
        INTEGER, INTENT(IN) :: ncid, it
        INTEGER :: ierr

        ierr = NF90_PUT_ATT(ncid, NF90_GLOBAL, 'timestep', it)
        CALL check( ierr )
        ierr = NF90_PUT_ATT(ncid, NF90_GLOBAL, 'workflow', config%workflow() )
        CALL check( ierr )
        ! May more shall follow

    END SUBROUTINE add_metadata


    ! -------------------------------------------------------------------------
    !
    ! -------------------------------------------------------------------------
    SUBROUTINE my_nf90_end_def( ncid )
        USE netcdf, ONLY : NF90_ENDDEF
        INTEGER, INTENT(IN) :: ncid
        INTEGER :: nf90_err

        nf90_err = NF90_ENDDEF( ncid )

        CALL CHECK( nf90_err, 'while taking netcdf out of define mode!')

    END SUBROUTINE my_nf90_end_def



    IMPURE ELEMENTAL SUBROUTINE my_nf90_def_var( ncid, dimidx, dimidy, dimidz, &
                                                 dimidt, attribute, varid )
        USE netcdf, ONLY : nf90_def_var
        INTEGER, INTENT(IN) :: ncid, dimidx, dimidy, dimidz, dimidt
        CHARACTER(*), INTENT(IN) :: attribute
        INTEGER, INTENT(OUT) :: varid
        INTEGER :: nf90_err

        nf90_err = nf90_def_var( ncid, attribute, my_nf90_real, &
                                 [dimidx, dimidy, dimidz, dimidt], varid )
        CALL check( nf90_err )

    END SUBROUTINE my_nf90_def_var


    ! Define variables
    SUBROUTINE my_nf90_put_var( ncid, varid, field )
        USE netcdf, ONLY : nf90_put_var
        INTEGER, INTENT(IN) :: ncid, varid
        REAL(real_kind), INTENT(IN) :: field(:,:,:,:,:,:)
        INTEGER :: nf90_err, nx, ny, nz, lpd, shape_r3(4), start(4)

        nx  = UBOUND(field,DIM=1)
        ny  = UBOUND(field,DIM=2)
        nz  = UBOUND(field,DIM=3)
        lpd = UBOUND(field,DIM=4)-1

        shape_r3(:) = [ nx*lpd+1, &
                        ny*lpd+1, &
                        nz*lpd+1, 1 ]

        start(:) = [ config%mpi_x_coord(), &
                     config%mpi_y_coord(), &
                    ! Why the hell are z-coordinates upside down?
                    ! I don't really understand what is going on here!
                     config%pz() - 1 - config%mpi_z_coord(), 1 ]

        start(:) = start*(shape_r3-1) + 1

        ! Borders are overlapping; this is not a good approach
        ! and needs some improvements
        nf90_err = nf90_put_var( ncid, varid, &
                                 pack_r6_to_r3( field ), &
                                 start=start, count=shape_r3 )
        CALL check( nf90_err )

    END SUBROUTINE my_nf90_put_var


        !ccc(:) = count+1
        !IF ( config%mpi_south() < 0 ) &
        !    count(1) = count(1) + 1
        !IF ( config%mpi_east() < 0 ) &
        !    count(2) = count(2) + 1
        !IF ( config%mpi_top() < 0 ) &
        !    count(3) = count(3) + 1
        !slice => r3_temp(1:count(1),1:count(2),:)

    ! -------------------------------------------------------------------------
    !
    ! -------------------------------------------------------------------------
    SUBROUTINE my_nf90_close( ncid )
        USE netcdf, ONLY : NF90_CLOSE
        INTEGER, INTENT(IN) :: ncid
        INTEGER :: nf90_err

        nf90_err = NF90_CLOSE( ncid )
        CALL CHECK( nf90_err, 'while closing netcdf-file!')

    END SUBROUTINE my_nf90_close


    ! -------------------------------------------------------------------------
    ! Check netCDF procedures for errors
    ! -------------------------------------------------------------------------
    SUBROUTINE check( status, msg )
        USE netcdf, ONLY : NF90_NOERR, NF90_STRERROR
        USE error_mod, ONLY : abort
        INTEGER, INTENT(IN) :: status
        CHARACTER(LEN=*), OPTIONAL :: msg

        IF(STATUS /= NF90_NOERR) THEN
            IF ( PRESENT(msg) ) THEN
                CALL abort( '"'//TRIM( NF90_STRERROR(status) )//'" '//msg )
            ELSE
                CALL abort( TRIM( NF90_STRERROR(status) ) )
            END IF
        END IF
    END SUBROUTINE check


    ! -------------------------------------------------------------------------
    ! Subroutine packing rank-6 arrays into rank-3 arrays
    ! -------------------------------------------------------------------------
    PURE FUNCTION pack_r6_to_r3(r6) RESULT(r3)
        INTEGER :: i, j, k , il, jl, kl ! Indices
        INTEGER :: i3, j3, k3           ! New r3 indices
        INTEGER :: nx, ny, nz, lpd      ! Number of elements and nodes
        REAL(real_kind), INTENT(IN) :: r6(:,:,:,:,:,:) ! 6d field
        REAL(my_nf90_real_kind), ALLOCATABLE :: r3(:,:,:) ! 3d field

        nx  = UBOUND(r6,DIM=1)
        ny  = UBOUND(r6,DIM=2)
        nz  = UBOUND(r6,DIM=3)
        lpd = UBOUND(r6,DIM=4)

        ALLOCATE( r3( nx*(lpd-1)+1, &
                      ny*(lpd-1)+1, &
                      nz*(lpd-1)+1) )

        ! Pack r6 field to r3 field
        ! Total number of coordinates in each direction: nx*(lpd-1)+1

        DO i=1,nx
            DO j=1,ny
                DO k=1,nz
                    DO il=1,lpd
                        DO jl=1,lpd
                            DO kl=1,lpd

                                ! Coordinate transformation
                                i3 = il*i+(i-1)*(lpd-il-1)
                                j3 = jl*j+(j-1)*(lpd-jl-1)
                                k3 = kl*k+(k-1)*(lpd-kl-1)

                                ! Store values in new field
                                r3(i3,j3,k3) = &
                                    REAL(r6(i,j,k,il,jl,kl), my_nf90_real_kind)

                            ENDDO
                        ENDDO
                    ENDDO
                ENDDO
            ENDDO
        ENDDO

    END FUNCTION pack_r6_to_r3

#endif
END MODULE netcdf_parallel_io_mod
