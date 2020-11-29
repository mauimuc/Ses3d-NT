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
!! $Date: 2013-12-02 18:27:28 +0100 (Mon, 02 Dec 2013) $
!! $Author: mauerberger $
!! $Revision: 833 $
!! @copyright GNU General Public License version 3 or later



!> Grid module
MODULE grid_mod
    USE error_mod
    USE parameters_mod

    IMPLICIT NONE

    PRIVATE

    ! Class, handling parallel and grid parameters
    TYPE :: grid_typ
        INTEGER :: pp_(3) ! Number of cpus px, py and pz
        INTEGER :: mpi_rank_
        INTEGER :: n_cpus_
        INTEGER :: mpi_comm_cart_
        INTEGER :: mpi_cart_coords_(3)
        INTEGER :: mpi_north_
        INTEGER :: mpi_south_
        INTEGER :: mpi_west_
        INTEGER :: mpi_east_
        INTEGER :: mpi_bottom_
        INTEGER :: mpi_top_
        INTEGER :: nx_glb_ !< Number of elements in x-direction
        INTEGER :: ny_glb_ !< Number of elements in x-direction
        INTEGER :: nz_glb_ !< Number of elements in x-direction
        INTEGER :: lpd_ !< Lagrange-Polynomial-Degree for x, y and z-directions
        INTEGER :: taper_width
        REAL(real_kind):: taper_slope
        CHARACTER(fnl) :: log_file_name_
        CHARACTER(fnl) :: log_file_dir_
        CHARACTER(16) :: event_name_
        CHARACTER(16) :: workflow_
        INTEGER :: log_file_unit
        REAL(real_kind) :: lat_min_, lat_max_
        REAL(real_kind) :: lon_min_, lon_max_
        REAL(real_kind) :: z_glb_min_, z_glb_max_
        INTEGER :: model_type_
        CHARACTER(fnl) :: rhoinv_
        CHARACTER(fnl) :: lambda_
        CHARACTER(fnl) :: mu_
        CHARACTER(fnl) :: a_
        CHARACTER(fnl) :: b_
        CHARACTER(fnl) :: c_
        CHARACTER(fnl) :: q_
        LOGICAL :: override_
        INTEGER :: nt_ !< number of time-steps
        REAL(real_kind) :: dt_ !< time-increment
        REAL(real_kind) :: fc_loc_max_ = 0.0
        REAL(real_kind) :: fc_glb_max_ = 0.0
        REAL(real_kind) :: cfl_loc_ = 0.0
        REAL(real_kind) :: cfl_glb_ = 0.0
        INTEGER :: date_time_(8)
    CONTAINS
        PROCEDURE, NON_OVERRIDABLE :: init_grid
        PROCEDURE :: init => init_grid
        PROCEDURE :: print
        PROCEDURE, NON_OVERRIDABLE, PRIVATE :: init_mpi_cart_topo
        PROCEDURE, NON_OVERRIDABLE :: mpi_rank
        PROCEDURE, NON_OVERRIDABLE :: mpi_size
        PROCEDURE, NON_OVERRIDABLE :: mpi_comm_cart
        PROCEDURE, NON_OVERRIDABLE :: px
        PROCEDURE, NON_OVERRIDABLE :: py
        PROCEDURE, NON_OVERRIDABLE :: pz
        PROCEDURE, NON_OVERRIDABLE :: mpi_x_coord
        PROCEDURE, NON_OVERRIDABLE :: mpi_y_coord
        PROCEDURE, NON_OVERRIDABLE :: mpi_z_coord
        PROCEDURE, NON_OVERRIDABLE :: mpi_east
        PROCEDURE, NON_OVERRIDABLE :: mpi_west
        PROCEDURE, NON_OVERRIDABLE :: mpi_north
        PROCEDURE, NON_OVERRIDABLE :: mpi_south
        PROCEDURE, NON_OVERRIDABLE :: mpi_top
        PROCEDURE, NON_OVERRIDABLE :: mpi_bottom
        PROCEDURE, NON_OVERRIDABLE :: n_glb
        PROCEDURE, NON_OVERRIDABLE :: n_loc
        PROCEDURE, NON_OVERRIDABLE :: lpd
        PROCEDURE, NON_OVERRIDABLE :: nx_loc
        PROCEDURE, NON_OVERRIDABLE :: ny_loc
        PROCEDURE, NON_OVERRIDABLE :: nz_loc
        PROCEDURE, NON_OVERRIDABLE :: nx_glb
        PROCEDURE, NON_OVERRIDABLE :: ny_glb
        PROCEDURE, NON_OVERRIDABLE :: nz_glb
        PROCEDURE, NON_OVERRIDABLE :: shape_glb
        PROCEDURE, NON_OVERRIDABLE :: shape_loc
    END TYPE grid_typ

    PUBLIC :: grid_typ, nx_loc, ny_loc, nz_loc, lpd

CONTAINS


    !> This function calculates the global number of elements
    ELEMENTAL INTEGER FUNCTION n_glb( self )
        CLASS(grid_typ), INTENT(IN) :: self
        n_glb = self%nx_glb_ * self%ny_glb_ * self%nz_glb_
    END FUNCTION n_glb



    !> This function calculates the local number of elements
    ELEMENTAL INTEGER FUNCTION n_loc( self )
        CLASS(grid_typ), INTENT(IN) :: self
        n_loc = (self%nx_loc()+1)*(self%ny_loc()+1)*(self%nz_loc()+1)
    END FUNCTION n_loc



    !> This function returns mpi_rank
    ELEMENTAL INTEGER FUNCTION mpi_rank( self )
        CLASS(grid_typ), INTENT(IN) :: self
        mpi_rank = self%mpi_rank_
    END FUNCTION mpi_rank



    !> This function returns size of MPI_COMM_WORLD
    ELEMENTAL INTEGER FUNCTION mpi_size( self )
        CLASS(grid_typ), INTENT(IN) :: self
        mpi_size = self%n_cpus_
    END FUNCTION mpi_size



    !> This function returns MPI Cartesian communicator
    ELEMENTAL INTEGER FUNCTION mpi_comm_cart( self )
        CLASS(grid_typ), INTENT(IN) :: self
        mpi_comm_cart = self%mpi_comm_cart_
    END FUNCTION mpi_comm_cart



    !> This function calculates the Cartesian MPI-Topology coordinates
    !! and returns mpi_x-coord
    ELEMENTAL INTEGER FUNCTION mpi_x_coord( self )
        CLASS( grid_typ ), INTENT( IN ) :: self
        mpi_x_coord = self%mpi_cart_coords_(1)
    END FUNCTION mpi_x_coord


    !> This function calculates the Cartesian MPI-Topology coordinates
    !! and returns mpi_y-coord
    ELEMENTAL INTEGER FUNCTION mpi_y_coord( self )
        CLASS( grid_typ ), INTENT( IN ) :: self
        mpi_y_coord = self%mpi_cart_coords_(2)
    END FUNCTION mpi_y_coord


    !> This function calculates the Cartesian MPI-Topology coordinates
    !! and returns mpi_z-coord
    ELEMENTAL INTEGER FUNCTION mpi_z_coord( self )
        CLASS( grid_typ ), INTENT( IN ) :: self
        mpi_z_coord = self%mpi_cart_coords_(3)
    END FUNCTION mpi_z_coord

    !> This function returns Lagrange polynomial degree (Number of collocation points: lpd+1)
    ELEMENTAL INTEGER FUNCTION lpd( self )
        CLASS(grid_typ), INTENT(IN) :: self
        lpd = self%lpd_
    END FUNCTION lpd


    !> This function returns number of processes in x-direction
    ELEMENTAL INTEGER FUNCTION px( self )
        CLASS(grid_typ), INTENT(IN) :: self
        px = self%pp_(1)
    END FUNCTION px


    !> This function returns number of processes in y-direction
    ELEMENTAL INTEGER FUNCTION py( self )
        CLASS(grid_typ), INTENT(IN) :: self
        py = self%pp_(2)
    END FUNCTION py


    !> This function returns number of processes in z-direction
    ELEMENTAL INTEGER FUNCTION pz( self )
        CLASS(grid_typ), INTENT(IN) :: self
        pz = self%pp_(3)
    END FUNCTION pz


    !> This function returns global number of elements in x-direction
    ELEMENTAL INTEGER FUNCTION nx_glb( self )
        CLASS(grid_typ), INTENT(IN) :: self
        nx_glb = self%nx_glb_
    END FUNCTION nx_glb


    !> This function returns global number of elements in y-direction
    ELEMENTAL INTEGER FUNCTION ny_glb( self )
        CLASS(grid_typ), INTENT(IN) :: self
        ny_glb = self%ny_glb_
    END FUNCTION ny_glb


    !> This function returns global number of elements in z-direction

    ELEMENTAL INTEGER FUNCTION nz_glb( self )
        CLASS(grid_typ), INTENT(IN) :: self
        nz_glb = self%nz_glb_
    END FUNCTION nz_glb



    !> This function returns local number of elements in x-direction
    ELEMENTAL INTEGER FUNCTION nx_loc( self )
        CLASS(grid_typ), INTENT(IN) :: self
        nx_loc = self%nx_glb_ / self%pp_(1) - 1
    END FUNCTION nx_loc


    !> This function returns local number of elements in y-direction
    ELEMENTAL INTEGER FUNCTION ny_loc( self )
        CLASS(grid_typ), INTENT(IN) :: self
        ny_loc = self%ny_glb_ / self%pp_(2) - 1
    END FUNCTION ny_loc


    !> This function returns local number of elements in z-direction
    ELEMENTAL INTEGER FUNCTION nz_loc( self )
        CLASS(grid_typ), INTENT(IN) :: self
        nz_loc = self%nz_glb_ / self%pp_(3) - 1
    END FUNCTION nz_loc


    PURE FUNCTION shape_glb( self )
        CLASS(grid_typ), INTENT(IN) :: self
        INTEGER :: shape_glb(6)
        shape_glb(:) = [ self%nx_glb_, self%ny_glb_, self%nz_glb_, &
                         self%lpd_+1,  self%lpd_+1,  self%lpd_+1  ]
    END FUNCTION shape_glb

    PURE FUNCTION shape_loc( self )
        CLASS(grid_typ), INTENT(IN) :: self
        INTEGER :: shape_loc(6)
        shape_loc(:) = [ self%nx_loc(), self%ny_loc(), self%nz_loc(), &
                         self%lpd_+1,   self%lpd_+1,   self%lpd_+1  ]
    END FUNCTION shape_loc


    ELEMENTAL INTEGER FUNCTION mpi_east( self )
        CLASS(grid_typ), INTENT(IN) :: self
        mpi_east = self%mpi_east_
    END FUNCTION mpi_east
    ELEMENTAL INTEGER FUNCTION mpi_west( self )
        CLASS(grid_typ), INTENT(IN) :: self
        mpi_west = self%mpi_west_
    END FUNCTION mpi_west
    ELEMENTAL INTEGER FUNCTION mpi_north( self )
        CLASS(grid_typ), INTENT(IN) :: self
        mpi_north = self%mpi_north_
    END FUNCTION mpi_north
    ELEMENTAL INTEGER FUNCTION mpi_south( self )
        CLASS(grid_typ), INTENT(IN) :: self
        mpi_south = self%mpi_south_
    END FUNCTION mpi_south
    ELEMENTAL INTEGER FUNCTION mpi_top( self )
        CLASS(grid_typ), INTENT(IN) :: self
        mpi_top = self%mpi_top_
    END FUNCTION mpi_top
    ELEMENTAL INTEGER FUNCTION mpi_bottom( self )
        CLASS(grid_typ), INTENT(IN) :: self
        mpi_bottom = self%mpi_bottom_
    END FUNCTION mpi_bottom



    !> This subroutine setups a MPI-Cartesian-communicator together with its resulting mpi-
    !! coordinates and adjacent ranks
    SUBROUTINE init_mpi_cart_topo( self )
#ifndef MPI_INCLUDE
        !! RECOMMENDED !!
        ! in case MPI_INCLUDE is NOT defined the mpi module will be used
        USE mpi
#else
        !! DEPRECATED !!
        ! when MPI_INCLUDE is defined, mpi header-files will be included
        INCLUDE 'mpif.h'
#endif
        CLASS(grid_typ), INTENT(INOUT) :: self
        INTEGER :: mpi_err

        ! Setup a Cartesian communicator with information on its Cartesian
        ! topology attached
        ! Partitioning has already been don in init-subroutine
        CALL MPI_CART_CREATE( MPI_COMM_WORLD, &
            3, &                           ! 3-dimensions
            self%pp_, &                    ! Extents px, py, pz
            (/.FALSE.,.FALSE.,.FALSE./), & ! Non periodic grid
            .FALSE., &                     ! No reordering
            self%mpi_comm_cart_, &         ! New communicator handle
            mpi_err )

        ! Determines MPI-Topology coordinates
        ! array of ( my_mpi_x-coord, my_mpi_y-coord, my_mpi_z-coord )
        CALL MPI_CART_COORDS( self%mpi_comm_cart_, &
            self%mpi_rank_, &
            3, &
            self%mpi_cart_coords_(:), &
            mpi_err )

        ! Why the hell are z-coordinates upside down?
        self%mpi_cart_coords_(3) = self%pp_(3) - 1 - self%mpi_cart_coords_(3)

        ! Find neighbouring ranks in x-direction
        CALL MPI_CART_SHIFT( self%mpi_comm_cart_, &
            0,  1, self%mpi_north_, self%mpi_south_, mpi_err )
        ! Find neighbouring ranks in y-direction
        CALL MPI_CART_SHIFT( self%mpi_comm_cart_, &
            1,  1, self%mpi_west_, self%mpi_east_, mpi_err )
        ! Find neighbouring ranks in z-direction
        CALL MPI_CART_SHIFT( self%mpi_comm_cart_, &
            2,  1, self%mpi_bottom_, self%mpi_top_, mpi_err )

    END SUBROUTINE init_mpi_cart_topo



    !> This wrapper subroutine printing all parameter-grout into a unit-number
    SUBROUTINE print( self, unit )
        CLASS(grid_typ), INTENT(IN) :: self
        INTEGER, INTENT(IN) :: unit

        ! Writes parallel-parameters
        CALL print_parallel(self, unit)
        ! Writes grid-parameters
        CALL print_grid(self, unit)

    END SUBROUTINE print


    !> This subroutine performs a formatted write of most parallel-parameters into a unit-number
    SUBROUTINE print_parallel( self, unit )
        CLASS(grid_typ), INTENT(IN) :: self
        INTEGER, INTENT(IN) :: unit
        CHARACTER(LEN=fsl) :: fmt

        fmt = '(A)'
        WRITE( UNIT=unit, FMT=fmt ) 'Parallel:'
        fmt = '( 4x, A, I0)'
        WRITE( UNIT=unit, FMT=fmt ) 'MPI rank = ', self%mpi_rank_
        fmt = '( 4x, "px = " I0, ", py = ", I0, ", pz = ", I0 )'
        WRITE( UNIT=unit, FMT=fmt ) self%pp_(:)
        fmt = '( 4x, "mpi x-coord = " I0, ", mpi y-coord = ", I0, ", mpi z-coord = ", I0 )'
        WRITE( UNIT=unit, FMT=fmt ) self%mpi_cart_coords_
        fmt = '( 4x, "northern rank = ", I0, ", southern rank = ", I0 )'
        WRITE( UNIT=unit, FMT=fmt ) self%mpi_north_, self%mpi_south_
        fmt = '( 4x, "western rank = ", I0, ", easter rank = ", I0 )'
        WRITE( UNIT=unit, FMT=fmt ) self%mpi_west_, self%mpi_east_
        fmt = '( 4x, "top rank = ", I0, ", bottom rank = ", I0 )'
        WRITE( UNIT=unit, FMT=fmt ) self%mpi_top_, self%mpi_bottom_

    END SUBROUTINE print_parallel


    !> This subroutine performs a formatted write of most parallel-parameters into a unit-number
    SUBROUTINE print_grid( self, unit )
        CLASS(grid_typ), INTENT(IN) :: self
        INTEGER, INTENT(IN) :: unit
        CHARACTER(LEN=fsl) :: fmt

        fmt = '(A)'
        WRITE( UNIT=unit, FMT=fmt ) 'Grid:'
        fmt = '( 4x, 3((A, " = ", I0),:,",  ") )'
        WRITE( UNIT=unit, FMT=fmt ) 'nx_glb_', self%nx_glb_, &
                                    'ny_glb_', self%ny_glb_, &
                                    'nz_glb_', self%nz_glb_
        WRITE( UNIT=unit, FMT=fmt ) 'nx_loc', self%nx_loc(), &
                                    'ny_loc', self%ny_loc(), &
                                    'nz_loc', self%nz_loc()
        fmt = '( 4x, A, " = ", I0)'
        WRITE( UNIT=unit, FMT=fmt ) 'total number of elements', self%n_glb()
        WRITE( UNIT=unit, FMT=fmt ) 'local number of elements', self%n_loc()
        WRITE( UNIT=unit, FMT=fmt ) 'lagrange polynomial degree', self%lpd_
        WRITE( UNIT=unit, FMT=fmt ) 'Width of tapering boundaries', self%taper_width
        fmt = '( 4x, A, " = ", G0)'
        WRITE( UNIT=unit, FMT=fmt ) 'Slope of tapering boundaries', self%taper_slope

    END SUBROUTINE print_grid



    !> This subroutine initializes: partition, parallelization and double-check values
    SUBROUTINE init_grid( self )
        USE string_utilities_mod, ONLY : i2c
        USE auto_partitionize_mod, ONLY : auto_partitionize
#ifndef MPI_INCLUDE
        !! RECOMMENDED !!
        ! in case MPI_INCLUDE is NOT defined the mpi module will be used
        USE mpi
#else
        !! DEPRECATED !!
        ! when MPI_INCLUDE is defined, mpi header-files will be included
        INCLUDE 'mpif.h'
#endif
        CLASS(grid_typ), INTENT(INOUT) :: self
        INTEGER :: mpi_err

        ! Returns the size of the group
        CALL MPI_COMM_SIZE( MPI_COMM_WORLD, self%n_cpus_, mpi_err )

        ! Determines the rank of the calling process
        CALL MPI_COMM_RANK( MPI_COMM_WORLD, self%mpi_rank_, mpi_err )

        ! Try to auto-partitionize grid
        self%pp_(:) = auto_partitionize( self%nx_glb_, &
                                         self%ny_glb_, &
                                         self%nz_glb_, &
                                         self%n_cpus_ )

        ! Check if properly partitioned
        IF ( MOD( self%nx_glb_, self%pp_(1) ) /= 0 ) &
            CALL abort( 'The ratio nx/px must be a whole number!' )
        IF ( MOD( self%ny_glb_, self%pp_(2) ) /= 0 ) &
            CALL abort( 'The ratio ny/py must be a whole number!' )
        IF ( MOD( self%nz_glb_, self%pp_(3) ) /= 0 ) &
            CALL abort( 'The ratio nz/pz must be a whole number!' )
        IF ( PRODUCT( self%pp_(:) ) /= self%n_cpus_ ) &
            CALL abort( 'Not properly partitionized! ' // i2c(self%n_cpus_ ) &
                        // '!='  // i2c( PRODUCT( self%pp_(:) ) ) )

        ! Check if there are at least two elements assigned into each direction
        IF ( self%nx_loc() < 2 ) &
            CALL abort( 'Local number of elements in x-direction must be &
                         & greater than 2!' )
        IF ( self%ny_loc() < 2 ) &
            CALL abort( 'Local number of elements in x-direction must be &
                         & greater than 2!' )
        IF ( self%nz_loc() < 2 ) &
            CALL abort( 'Local number of elements in x-direction must be &
                         & greater than 2!' )

        ! Initialize mpi Cartesian topology
        CALL self%init_mpi_cart_topo()

    END SUBROUTINE init_grid




END MODULE grid_mod
