!> @file
!! Ses3d-NT - simulation of elastic wave propagation in spherical sections
!!
!! (c) by Stefan Mauerberger
!!    and Maksym Melnyk
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
!! $Date: 2013-09-16 15:37:13 +0200 (Mon, 16 Sep 2013) $
!! $Author: mauerberger $
!! $Revision: 686 $
!! @copyright GNU General Public License version 3 or later



!> Auto partitionize module
MODULE auto_partitionize_mod

    IMPLICIT NONE

    PRIVATE

    PUBLIC :: auto_partitionize

CONTAINS

    !> This function partitionizes a grid (nx,ny,nz) into n_cpus pieces
    !> @param nx Number of elements in theta-direction
    !> @param ny Number of elements in phi-direction
    !> @param nz Number of elements in r-direction
    !> @param n_cpus Number of CPUs
    !!
    !> @return an array [px,py,pz].
    FUNCTION auto_partitionize( nx, ny, nz, n_cpus ) RESULT( pp )
        USE error_mod, ONLY : abort
        INTEGER, INTENT(IN) :: n_cpus, nx, ny, nz
        INTEGER :: pp(3)
        INTEGER :: px, py, pz, p_yz, cs

        ! Initialize cross-section to HUGE
        cs = HUGE(1)
        !
        pp(:) = -1

        ! Obtain number of partitions in px-direction
        ! Using brute-force factorization by division
        ! Iterate through all possible px values
        DO px = 1, n_cpus
            ! If px does not divide n_cpus as a whole number => cycle
            IF ( MOD( n_cpus, px ) /= 0 ) &
                CYCLE
            ! Remaining number of cpus; will partitionize (ny,nz)
            p_yz = n_cpus / px
            ! Iterate through all remaining py values
            DO py = 1, p_yz
                ! If py does not divide p_yz as a whole number => cycle
                IF ( MOD( p_yz, py ) /= 0 ) &
                    CYCLE
                ! Determine pz
                pz = p_yz / py
                ! There must be more than one element per process
                IF ( (px==nx) .OR. (py==ny) .OR. (pz==nz) ) &
                    CYCLE
                ! Check if factorization splits grid into whole-numbers
                IF ( ( MOD( nx, px ) == 0 ) .AND. &
                     ( MOD( ny, py ) == 0 ) .AND. &
                     ( MOD( nz, pz ) == 0 ) .AND. &
                     ! Check if cross-section is smaller than previous
                     ( px*ny*nz + py*nx*nz + pz*nx*ny < cs ) ) THEN
                    ! Save new, smaller cross-section value
                    cs = px*ny*nz + py*nx*nz + pz*nx*ny
                    pp(:) = [ px, py, pz ]
                END IF
            END DO
        END DO

        IF ( ANY( pp(:) < 1 ) ) &
            CALL abort( 'Auto-partition failed! Choos differnt number of CPUs &
                         & or alter nx, ny or nz values.' )

    END FUNCTION auto_partitionize

END MODULE auto_partitionize_mod
