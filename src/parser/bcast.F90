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
!! $Date: 2013-08-14 15:24:51 +0200 (Wed, 14 Aug 2013) $
!! $Author: mauerberger $
!! $Revision: 576 $
!! @copyright GNU General Public License version 3 or later



!> Module to broadcast data amongst processes
!!
!! This module is dedicated to ...
MODULE bcast_mod
    USE error_mod, ONLY : abort
    USE parameters_mod, ONLY : real_kind, my_mpi_real, my_mpi_integer
    USE parser_mod
#ifndef MPI_INCLUDE
    !! RECOMMENDED !!
    ! in case MPI_INCLUDE is NOT defined the mpi module will be used
    USE mpi
#endif
    IMPLICIT NONE
#ifdef MPI_INCLUDE
    !! DEPRECATED !!
    ! when MPI_INCLUDE is defined, mpi header-files will be included
    INCLUDE 'mpif.h'
#endif

    PRIVATE

    INTERFACE bcast_namelist
        MODULE PROCEDURE bcast_grid_typ
        MODULE PROCEDURE bcast_grad_typ
        MODULE PROCEDURE bcast_general_typ
        MODULE PROCEDURE bcast_model_typ
        MODULE PROCEDURE bcast_time_typ
        MODULE PROCEDURE bcast_source_cls
        MODULE PROCEDURE bcast_source_sac_typ
        MODULE PROCEDURE bcast_receiver_typ
        MODULE PROCEDURE bcast_output_cls
    END INTERFACE

    PUBLIC :: bcast_namelist, bcast_integer, reduce_sum

CONTAINS

    SUBROUTINE bcast_model_typ(model, mr)
        TYPE(model_typ), INTENT(INOUT) :: model
        INTEGER, INTENT(IN) :: mr
        INTEGER :: mpi_err

        CALL MPI_BCAST(model%lat_max     , 1,   my_mpi_real,    &
            mr, MPI_COMM_WORLD, mpi_err)
        CALL MPI_BCAST(model%lat_min     , 1,   my_mpi_real,    &
            mr, MPI_COMM_WORLD, mpi_err)
        CALL MPI_BCAST(model%lon_max     , 1,   my_mpi_real,    &
            mr, MPI_COMM_WORLD, mpi_err)
        CALL MPI_BCAST(model%lon_min     , 1,   my_mpi_real,    &
            mr, MPI_COMM_WORLD, mpi_err)
        CALL MPI_BCAST(model%rad_max   , 1,   my_mpi_real,    &
            mr, MPI_COMM_WORLD, mpi_err)
        CALL MPI_BCAST(model%rad_min   , 1,   my_mpi_real,    &
            mr, MPI_COMM_WORLD, mpi_err)
        CALL MPI_BCAST(model%model_type, 1,   my_mpi_integer, &
            mr, MPI_COMM_WORLD, mpi_err)
        CALL MPI_BCAST(model%rhoinv      , LEN(model%rhoinv), mpi_character,  &
            mr, MPI_COMM_WORLD, mpi_err)
        CALL MPI_BCAST(model%lambda      , LEN(model%lambda), mpi_character,  &
            mr, MPI_COMM_WORLD, mpi_err)
        CALL MPI_BCAST(model%mu          , LEN(model%mu    ), mpi_character,  &
            mr, MPI_COMM_WORLD, mpi_err)
        CALL MPI_BCAST(model%a           , LEN(model%a     ), mpi_character,  &
            mr, MPI_COMM_WORLD, mpi_err)
        CALL MPI_BCAST(model%b           , LEN(model%b     ), mpi_character,  &
            mr, MPI_COMM_WORLD, mpi_err)
        CALL MPI_BCAST(model%c           , LEN(model%c     ), mpi_character,  &
            mr, MPI_COMM_WORLD, mpi_err)
        CALL MPI_BCAST(model%q           , LEN(model%q     ), mpi_character,  &
            mr, MPI_COMM_WORLD, mpi_err)
        CALL MPI_BCAST(model%override    , 1                , mpi_logical  ,  &
            mr, MPI_COMM_WORLD, mpi_err)
    END SUBROUTINE bcast_model_typ

    SUBROUTINE bcast_general_typ(general, mr)
        TYPE(general_typ), INTENT(INOUT) :: general
        INTEGER, INTENT(IN) :: mr
        INTEGER :: mpi_err

        CALL MPI_BCAST(general%workflow    , LEN(general%workflow),     &
            mpi_character, mr, MPI_COMM_WORLD, mpi_err)
        CALL MPI_BCAST(general%log_file_dir, LEN(general%log_file_dir), &
            mpi_character, mr, MPI_COMM_WORLD, mpi_err)
        CALL MPI_BCAST(general%event_name  , LEN(general%event_name),   &
            mpi_character,  mr, MPI_COMM_WORLD, mpi_err)

    END SUBROUTINE bcast_general_typ

    SUBROUTINE bcast_grid_typ(grid, mr)
        TYPE(grid_typ), INTENT(INOUT) :: grid
        INTEGER, INTENT(IN) :: mr
        INTEGER :: mpi_err

        CALL MPI_BCAST(grid%nx,  1, my_mpi_integer, mr, MPI_COMM_WORLD, mpi_err)
        CALL MPI_BCAST(grid%ny,  1, my_mpi_integer, mr, MPI_COMM_WORLD, mpi_err)
        CALL MPI_BCAST(grid%nz,  1, my_mpi_integer, mr, MPI_COMM_WORLD, mpi_err)
        CALL MPI_BCAST(grid%lpd, 1, my_mpi_integer, mr, MPI_COMM_WORLD, mpi_err)
        CALL MPI_BCAST(grid%taper_width, 1, my_mpi_integer, mr, MPI_COMM_WORLD, mpi_err)
        CALL MPI_BCAST(grid%taper_slope, 1, my_mpi_real, mr, MPI_COMM_WORLD, mpi_err)

    END SUBROUTINE bcast_grid_typ

    SUBROUTINE bcast_source_cls(src, mr)
        CLASS(source_cls), INTENT(INOUT) :: src
        INTEGER, INTENT(IN) :: mr
        INTEGER :: mpi_err

        CALL MPI_BCAST(src%lat, 1, my_mpi_real, &
            mr, MPI_COMM_WORLD, mpi_err)
        CALL MPI_BCAST(src%lon, 1, my_mpi_real, &
            mr, MPI_COMM_WORLD, mpi_err)
        CALL MPI_BCAST(src%depth, 1, my_mpi_real, &
            mr, MPI_COMM_WORLD, mpi_err)
        CALL MPI_BCAST(src%wavelet, LEN(src%wavelet), &
            mpi_character, mr, MPI_COMM_WORLD, mpi_err)
        CALL MPI_BCAST(src%onset, 1, my_mpi_real, &
            mr, MPI_COMM_WORLD, mpi_err)
        CALL MPI_BCAST(src%width, 1, my_mpi_real, &
            mr, MPI_COMM_WORLD, mpi_err)

        SELECT TYPE (src)
            TYPE IS (source_sf_typ)
                CALL MPI_BCAST(src%direction(:), 3, my_mpi_real, &
                    mr, MPI_COMM_WORLD, mpi_err)
            TYPE IS (source_mt_typ)
                CALL MPI_BCAST(src%moment_tensor(:), 6, my_mpi_real, &
                    mr, MPI_COMM_WORLD, mpi_err)
            CLASS DEFAULT
                CALL abort("source class not known in subroutine 'bcast'")
        END SELECT

    END SUBROUTINE bcast_source_cls

    SUBROUTINE bcast_source_sac_typ(src, mr)
        CLASS(source_sac_typ), INTENT(INOUT) :: src
        INTEGER, INTENT(IN) :: mr
        INTEGER :: mpi_err

        CALL MPI_BCAST(src%file_name, LEN(src%file_name), &
            mpi_character, mr, MPI_COMM_WORLD, mpi_err)

    END SUBROUTINE bcast_source_sac_typ

    SUBROUTINE bcast_time_typ(time, mr)
        TYPE(time_typ), INTENT(INOUT) :: time
        INTEGER, INTENT(IN) :: mr
        INTEGER :: mpi_err

        CALL MPI_BCAST(time%dt          , 1                   , my_mpi_real   , &
            mr, MPI_COMM_WORLD, mpi_err)
        CALL MPI_BCAST(time%nt          , 1                   , my_mpi_integer, &
            mr, MPI_COMM_WORLD, mpi_err)
        CALL MPI_BCAST(time%date_time(:), SIZE(time%date_time), my_mpi_integer, &
            mr, MPI_COMM_WORLD, mpi_err)

    END SUBROUTINE bcast_time_typ


    SUBROUTINE bcast_grad_typ(typ, mr)
        TYPE(grad_typ), INTENT(INOUT) :: typ
        INTEGER, INTENT(IN) :: mr
        INTEGER :: mpi_err
        INTEGER :: k
        CHARACTER(LEN=8) :: attribute

        CALL MPI_BCAST(typ%timestep_start, 1, my_mpi_integer, &
            mr, MPI_COMM_WORLD, mpi_err)
        CALL MPI_BCAST(typ%timestep_increment, 1, my_mpi_integer, &
            mr, MPI_COMM_WORLD, mpi_err)
        CALL MPI_BCAST(typ%timestep_end, 1, my_mpi_integer, &
            mr, MPI_COMM_WORLD, mpi_err)
        DO k=1, SIZE(typ%attributes(:))
            attribute = typ%attributes(k)
            CALL MPI_BCAST(attribute, 8, MPI_CHARACTER, &
                mr, MPI_COMM_WORLD, mpi_err)
            typ%attributes(k) = attribute
        END DO

    END SUBROUTINE bcast_grad_typ


    SUBROUTINE bcast_output_cls(typ, mr)
        CLASS(output_cls), INTENT(INOUT) :: typ
        INTEGER, INTENT(IN) :: mr
        INTEGER :: mpi_err
        INTEGER :: k
        CHARACTER(LEN=8) :: attribute

        CALL MPI_BCAST(typ%timestep_start, 1, my_mpi_integer, &
            mr, MPI_COMM_WORLD, mpi_err)
        CALL MPI_BCAST(typ%timestep_increment, 1, my_mpi_integer, &
            mr, MPI_COMM_WORLD, mpi_err)
        CALL MPI_BCAST(typ%timestep_end, 1, my_mpi_integer, &
            mr, MPI_COMM_WORLD, mpi_err)
        DO k=1, SIZE(typ%attributes(:))
            attribute = typ%attributes(k)
            CALL MPI_BCAST(attribute, 8, MPI_CHARACTER, &
                mr, MPI_COMM_WORLD, mpi_err)
            typ%attributes(k) = attribute
        END DO
        CALL MPI_BCAST(typ%override, 1, MPI_LOGICAL, &
            mr, MPI_COMM_WORLD, mpi_err)
        CALL MPI_BCAST(typ%prefix, LEN(typ%prefix), MPI_CHARACTER, &
            mr, MPI_COMM_WORLD, mpi_err)

    END SUBROUTINE bcast_output_cls


    SUBROUTINE bcast_receiver_typ(rec, mr)
        TYPE(receiver_typ), INTENT(INOUT) :: rec
        INTEGER, INTENT(IN) :: mr
        INTEGER :: mpi_err
        INTEGER :: k
        CHARACTER(LEN=8) :: attribute

        CALL MPI_BCAST(rec%lat, 1, my_mpi_real, &
            mr, MPI_COMM_WORLD, mpi_err)
        CALL MPI_BCAST(rec%lon, 1, my_mpi_real, &
            mr, MPI_COMM_WORLD, mpi_err)
        CALL MPI_BCAST(rec%depth, 1, my_mpi_real, &
            mr, MPI_COMM_WORLD, mpi_err)

        CALL MPI_BCAST(rec%network, LEN(rec%network), MPI_CHARACTER, &
            mr, MPI_COMM_WORLD, mpi_err)
        CALL MPI_BCAST(rec%station, LEN(rec%station), MPI_CHARACTER, &
            mr, MPI_COMM_WORLD, mpi_err)
        CALL MPI_BCAST(rec%location, LEN(rec%location), MPI_CHARACTER, &
            mr, MPI_COMM_WORLD, mpi_err)

        DO k=1, SIZE(rec%attributes(:))
            attribute = rec%attributes(k)
            CALL MPI_BCAST(attribute, 8, MPI_CHARACTER, &
                mr, MPI_COMM_WORLD, mpi_err)
            rec%attributes(k) = attribute
        END DO

        CALL MPI_BCAST(rec%override, 1, MPI_LOGICAL, &
            mr, MPI_COMM_WORLD, mpi_err)
        CALL MPI_BCAST(rec%prefix, LEN(rec%prefix), MPI_CHARACTER, &
            mr, MPI_COMM_WORLD, mpi_err)

    END SUBROUTINE bcast_receiver_typ

    SUBROUTINE bcast_integer(i, mr)
        INTEGER, INTENT(INOUT) :: i
        INTEGER, INTENT(IN) :: mr
        INTEGER :: mpi_err

        CALL MPI_Bcast(i, 1, my_mpi_integer, &
            mr, MPI_COMM_WORLD, mpi_err)

    END SUBROUTINE bcast_integer

    INTEGER FUNCTION reduce_sum(i, mr) RESULT(s)
        INTEGER, INTENT(IN) :: i
        INTEGER, INTENT(IN) :: mr
        INTEGER :: mpi_err

        CALL MPI_REDUCE( i, s, 1, my_mpi_integer, MPI_SUM, &
            mr, MPI_COMM_WORLD, mpi_err )

    END FUNCTION reduce_sum

END MODULE bcast_mod
