!> @file
!! Ses3d-NT - simulation of elastic wave propagation in spherical sections
!!
!! (c) by Gilbert Brietzke
!!    and Stefan Mauerberger
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
!! $Date: 2013-11-30 10:55:27 +0100 (Sat, 30 Nov 2013) $
!! $Author: mauerberger $
!! $Revision: 828 $
!! @copyright GNU General Public License version 3 or later



!> Module for communicating/exchanging values between adjacent elements.
!!
!! It provides two subroutines. One communicates values of adjacent
!! elements locally. The other exchanges values amongst processes/ranks.
MODULE communicate_fields_mod
    USE parameters_mod, ONLY : real_kind, my_mpi_real
    USE configuration_mod, ONLY : conf => configuration
    ! A GCC <= 4.8.1 workaround
    USE grid_mod, ONLY : nx_loc, ny_loc, nz_loc, lpd

    IMPLICIT NONE

    PRIVATE

    PUBLIC :: communicate_local_boundaries, communicate_global_boundaries

CONTAINS

    !> Communicate boundaries inside local process ranks
    !!
    !! @param field Any Ses3d-NT array of rank six.
    SUBROUTINE communicate_local_boundaries(field)
        ! Passed dummy variable
        REAL(real_kind), INTENT(INOUT) :: field(0:nx_loc(conf),0:ny_loc(conf),&
                                                0:nz_loc(conf),0:lpd(conf),&
                                                0:lpd(conf),0:lpd(conf))
        ! Local variables
        INTEGER :: i, j, k, nx_, ny_, nz_, lpd_

        nx_ = conf%nx_loc()
        ny_ = conf%ny_loc()
        nz_ = conf%nz_loc()
        lpd_ = conf%lpd()

        ! Communicate theta direction
        DO i=1, nx_
           field(  i,:,:,   0,:,:) = field(i,:,:,0,:,:) + field(i-1,:,:,lpd_,:,:)
           field(i-1,:,:,lpd_,:,:) = field(i,:,:,0,:,:)
        END DO

        ! Communicate phi direction
        DO j=1, ny_
           field(:,  j,:,:,   0,:) = field(:,j,:,:,0,:) + field(:,j-1,:,:,lpd_,:)
           field(:,j-1,:,:,lpd_,:) = field(:,j,:,:,0,:)
        ENDDO

        ! Communicate r direction
        DO k=1, nz_
           field(:,:,  k,:,:,   0) = field(:,:,k,:,:,0) + field(:,:,k-1,:,:,lpd_)
           field(:,:,k-1,:,:,lpd_) = field(:,:,k,:,:,0)
        END DO

    END SUBROUTINE communicate_local_boundaries


    !> Communication between processors
    !!
    !! @note Could be more efficient using nonblocking communication.
    !!       However, this is difficult because edges and corners have to be
    !!       exchanges explicitly. This results in a total of 26 MPI calls.
    !!        * 6 faces
    !!        * 12 edges
    !!        * 8 points
    !!
    !! @param field Any Ses3d-NT array of rank six.
    SUBROUTINE communicate_global_boundaries(field)
#ifndef MPI_INCLUDE
        !! RECOMMENDED !!
        ! in case MPI_INCLUDE is NOT defined the mpi module will be used
        USE mpi
#else
        !! DEPRECATED !!
        ! when MPI_INCLUDE is defined, mpi header-files will be included
        INCLUDE 'mpif.h'
#endif
        ! Passed dummy variable
        REAL(real_kind), INTENT(INOUT) :: field(0:nx_loc(conf),0:ny_loc(conf),&
                                                0:nz_loc(conf),0:lpd(conf),&
                                                0:lpd(conf),0:lpd(conf))
        ! Local variables
        INTEGER :: ierr, statuses(MPI_Status_size,4), requests(4), nx_, ny_, &
            nz_, lpd_, bottom, top, north, south, west, east, mpi_comm_cart
        REAL(real_kind) :: &
            s_south(0:ny_loc(conf),0:nz_loc(conf),0:lpd(conf),0:lpd(conf)), &
            s_north(0:ny_loc(conf),0:nz_loc(conf),0:lpd(conf),0:lpd(conf)), &
            s_west(0:nx_loc(conf),0:nz_loc(conf),0:lpd(conf),0:lpd(conf)), &
            s_east(0:nx_loc(conf),0:nz_loc(conf),0:lpd(conf),0:lpd(conf)), &
            s_bottom(0:nx_loc(conf),0:ny_loc(conf),0:lpd(conf),0:lpd(conf)), &
            s_top(0:nx_loc(conf),0:ny_loc(conf),0:lpd(conf),0:lpd(conf)), &
            r_south(0:ny_loc(conf),0:nz_loc(conf),0:lpd(conf),0:lpd(conf)), &
            r_north(0:ny_loc(conf),0:nz_loc(conf),0:lpd(conf),0:lpd(conf)), &
            r_west(0:nx_loc(conf),0:nz_loc(conf),0:lpd(conf),0:lpd(conf)), &
            r_east(0:nx_loc(conf),0:nz_loc(conf),0:lpd(conf),0:lpd(conf)), &
            r_bottom(0:nx_loc(conf),0:ny_loc(conf),0:lpd(conf),0:lpd(conf)), &
            r_top(0:nx_loc(conf),0:ny_loc(conf),0:lpd(conf),0:lpd(conf))

        ! Initialize local variables
        requests(:) = MPI_Request_null
        nx_ = conf%nx_loc()
        ny_ = conf%ny_loc()
        nz_ = conf%nz_loc()
        lpd_ = conf%lpd()
        bottom = conf%mpi_bottom()
        top = conf%mpi_top()
        north = conf%mpi_north()
        south = conf%mpi_south()
        west = conf%mpi_west()
        east = conf%mpi_east()
        mpi_comm_cart = conf%mpi_comm_cart()

        !!! Communication in x direction
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Send to and receive from southern neighbour
        IF ( south >= 0 ) THEN
            ! Buffer; prevents overriding edges and corners while sending
            s_south(:,:,:,:) = field(nx_,:,:,lpd_,:,:)
            ! Send local norther-face to western neighbouring rank
            CALL MPI_ISend( s_south, SIZE(s_south), my_mpi_real, south, &
                            110, mpi_comm_cart, requests(1), ierr )
            ! Receive from southern neighbouring rank
            CALL MPI_IRecv( r_south, SIZE(r_south), my_mpi_real, south, &
                            120, mpi_comm_cart, requests(2), ierr )
        ELSE
            ! If there is no neighbour
            r_south(:,:,:,:) = 0.0
        END IF
        ! Send to and receive from northern neighbour
        IF ( north >= 0 ) THEN
            ! Buffer; prevents overriding edges and corners while sending
            s_north = field(0,:,:,0,:,:)
            ! Send local southern-face to western neighbouring rank
            CALL MPI_ISend( s_north, SIZE(s_north), my_mpi_real, north, &
                            120, mpi_comm_cart, requests(3), ierr )
            ! Receive from northern neighbouring rank
            CALL MPI_IRecv( r_north, SIZE(r_north), my_mpi_real, north, &
                            110, mpi_comm_cart, requests(4), ierr )
        ELSE
            ! If there is no neighbour
            r_north(:,:,:,:) = 0.0
        END IF
        ! Wait till MPI communication has finished
        CALL MPI_Waitall( 4, requests(:), statuses(:,:), ierr )
        ! Exchange boundaries
        field(nx_,:,:,lpd_,:,:) = field(nx_,:,:,lpd_,:,:) + r_south
        field(0,:,:,0,:,:) = field(0,:,:,0,:,:) + r_north


        !!! Communication in y direction
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Send to and receive from eastern neighbour
        IF ( east >= 0 ) THEN
            ! Buffer; prevents overriding edges and corners while sending
            s_east = field(:,ny_,:,:,lpd_,:)
            ! Send local western-face to western neighbouring rank
            CALL MPI_ISend( s_east, SIZE(s_east), my_mpi_real, east, &
                            130, mpi_comm_cart, requests(1), ierr )
            ! Receive from eastern neighbouring rank
            CALL MPI_IRecv( r_east, SIZE(r_east), my_mpi_real, east, &
                            140, mpi_comm_cart, requests(2), ierr )
        ELSE
            ! If there is no neighbour
            r_east(:,:,:,:) = 0.0
        END IF
        ! Send to and receive from western neighbour
        IF ( west >= 0 ) THEN
            ! Buffer; prevents overriding edges and corners while sending
            s_west = field(:,0,:,:,0,:)
            ! Send local eastern-face to western neighbouring rank
            CALL MPI_ISend( s_west, SIZE(s_west), my_mpi_real, west, &
                            140, mpi_comm_cart, requests(3), ierr )
            ! Receive from western neighbouring rank
            CALL MPI_IRecv( r_west, SIZE(r_west), my_mpi_real, west, &
                            130, mpi_comm_cart, requests(4), ierr )
        ELSE
            ! If there is no neighbour
            r_west(:,:,:,:) = 0.0
        END IF
        ! Wait till MPI communication has finished
        CALL MPI_Waitall( 4, requests(:), statuses(:,:), ierr )
        ! Exchange boundaries
        field(:,ny_,:,:,lpd_,:) = field(:,ny_,:,:,lpd_,:) + r_east
        field(:,0,:,:,0,:) = field(:,0,:,:,0,:) + r_west


        !!! Communication in z direction
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Send to and receive from top neighbour
        IF ( top >= 0  ) THEN
            ! Buffer; prevents overriding edges and corners while sending
            s_top = field(:,:,nz_,:,:,lpd_)
            ! Send local bottom-face to top neighbouring rank
            CALL MPI_ISend( s_top, SIZE(s_top), my_mpi_real, top, &
                            150, mpi_comm_cart, requests(1), ierr )
            ! Receive from top neighbouring rank
            CALL MPI_IRecv( r_top, SIZE(r_top), my_mpi_real, top, &
                            160, mpi_comm_cart, requests(2), ierr )
        ELSE
            ! If there is no neighbour
            r_top(:,:,:,:) = 0.0
        END IF
        ! Send to and receive from bottom neighbour
        IF ( bottom >= 0 ) THEN
            ! Buffer; prevents overriding edges and corners while sending
            s_bottom = field(:,:,0,:,:,0)
            ! Send local top-face to bottom neighbouring rank
            CALL MPI_ISend( s_bottom, SIZE(s_bottom), my_mpi_real, bottom, &
                            160, mpi_comm_cart, requests(3), ierr )
            ! Receive from bottom neighbouring rank
            CALL MPI_IRecv( r_bottom, SIZE(r_bottom), my_mpi_real, bottom, &
                            150, mpi_comm_cart, requests(4), ierr )
        ELSE
            ! If there is no bottom neighbour
            r_bottom(:,:,:,:) = 0.0
        END IF
        ! Wait till MPI communication has finished
        CALL MPI_Waitall( 4, requests(:), statuses(:,:), ierr )
        ! Exchange boundaries
        field(:,:,nz_,:,:,lpd_) = field(:,:,nz_,:,:,lpd_) + r_top
        field(:,:,0,:,:,0) = field(:,:,0,:,:,0) + r_bottom

    END SUBROUTINE communicate_global_boundaries

END MODULE communicate_fields_mod
