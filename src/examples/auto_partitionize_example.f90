!> @example auto_partitionize_example.f90
!!
!! Example program demonstrating how Ses3d-NT's auto parallelization
!! scheme works. The program expects the number of elements assigned to
!! each direction theta `nx`, phi `nz` and r `nz` and the total number of
!! processes `cpus`. It that grid can be split into `cpus` partitions the
!! parallelization parameters are written to stdout otherwise an error is
!! reported.
!!
!! Compile by simply using the Makefile `make auto_partitionize_example`.
!!
!! Usage: `./auto_partitionize [nx] [ny] [nz] [n_cpus]`. If you want to
!! iterate over an interval of `cpus` you can do the following:
!!
!!     for i in {[start]..[stop]}; do ./auto_partitionize [nx] [ny] [nz] $i 2> /dev/null; done
!!
!! This will only list possible configurations. Error messages are redirected
!! `/dev/null`.
PROGRAM auto_partitionize_example
    USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY : OUTPUT_UNIT, ERROR_UNIT
    USE auto_partitionize_mod, ONLY : auto_partitionize

    IMPLICIT NONE

    INTEGER :: nx, ny, nz, n_cpus, pp(3)
    CHARACTER(LEN=8) :: arg(4)
    CHARACTER(LEN=*), PARAMETER :: fmt = '("px=", I0, " py=", I0, " pz=", I0, &
        &" - nx_loc=", I0, " ny_loc=", I0, " nz_loc=", I0, " - cpus=", I0)'

    IF ( COMMAND_ARGUMENT_COUNT() < 4 ) THEN
        WRITE(UNIT=ERROR_UNIT, FMT='(A)') 'Missing parameters!'
        WRITE(UNIT=ERROR_UNIT, FMT='(A)') 'Usage: auto_partitionize &
            &[nx] [ny] [nz] [n_cpus]'
        STOP 1
    END IF

    CALL GET_COMMAND_ARGUMENT(1, arg(1))
    CALL GET_COMMAND_ARGUMENT(2, arg(2))
    CALL GET_COMMAND_ARGUMENT(3, arg(3))
    CALL GET_COMMAND_ARGUMENT(4, arg(4))

    READ(UNIT=arg(:), FMT=*) nx, ny, nz, n_cpus

    pp = auto_partitionize(nx, ny, nz, n_cpus)

    WRITE(UNIT=OUTPUT_UNIT, FMT=fmt) pp(:), nx/pp(1), ny/pp(2), nz/pp(3), PRODUCT(pp)

END PROGRAM auto_partitionize_example

!> @file
!! $Date: 2013-11-27 13:56:48 +0100 (Wed, 27 Nov 2013) $
!! $Author: mauerberger $
!! $Revision: 819 $
!! @copyright GNU General Public License version 3 or later
