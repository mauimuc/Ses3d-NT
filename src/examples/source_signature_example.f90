!> @example source_signature_example.f90
!!
!! Example program which demonstrates how to use the source_signature module
!!
!! Compile by simply using the Makefile:
!!
!!     make source_signature_example
!!
!! You can use Python to visualize the generated file
!! @code{.py} 
!!import numpy as np
!!import matplotlib.pyplot as plt
!!
!!brune = np.loadtxt('brune.dat', delimiter=',')
!!gauss = np.loadtxt('gauss.dat', delimiter=',')
!!delta = np.loadtxt('delta.dat', delimiter=',')
!!sin_3 = np.loadtxt('sin_3.dat', delimiter=',')
!!ricker = np.loadtxt('ricker.dat', delimiter=',')
!!
!!plt.plot(brune[:,0],brune[:,1], label='brune')
!!plt.plot(gauss[:,0],gauss[:,1], label='gauss')
!!plt.plot(delta[:,0],delta[:,1], label='delta')
!!plt.plot(sin_3[:,0],sin_3[:,1], label='sin_3')
!!plt.plot(ricker[:,0],ricker[:,1], label='ricker')
!!plt.legend()
!!plt.show()
!! @endcode
!! 
!! @todo Annotation
PROGRAM source_signature_example
    USE parameters_mod, ONLY : real_kind
    USE source_signature_mod, ONLY : ricker, brune, gauss, delta, sin_3

    IMPLICIT NONE

    INTEGER, PARAMETER :: nt = 1500
    REAL(real_kind), PARAMETER :: dt = 0.02
    INTEGER :: i
    REAL(real_kind) :: time(nt), stf(nt)

    ! Make time vector
    time(:) = REAL( [(i, i=0, nt-1)], real_kind ) * dt

    ! Brune pulse
    stf(:) = brune(time=time(:), onset=dt, width=5.0_real_kind)
    CALL write_stf('brune.dat', dt, stf )

    ! Gaussian wavelet 
    stf(:) = gauss(time=time(:), onset=10.0_real_kind, width=4.0_real_kind)
    CALL write_stf('gauss.dat', dt, stf )

    ! 'Delta' pulse
    stf(:) = delta(time=time(:), onset=10*dt)
    CALL write_stf('delta.dat', dt, stf )

    ! Sin**3 wavelet 
    stf(:) = sin_3(time=time(:), onset=1.0_real_kind, width=10.0_real_kind)
    CALL write_stf('sin_3.dat', dt, stf )

    ! Mexican hat wavelet 
    stf(:) = ricker(time=time(:), onset=10.0_real_kind, width=5.0_real_kind)
    CALL write_stf('ricker.dat', dt, stf )

CONTAINS

    SUBROUTINE write_stf(file, dt, stf)
        CHARACTER(LEN=*) :: file
        REAL(real_kind) :: stf(:), dt
        INTEGER :: lun

        OPEN(NEWUNIT=lun, FILE=file)
        WRITE(UNIT=lun, FMT='("# ", A, I0)') 'nt = ', SIZE(stf)
        WRITE(UNIT=lun, FMT='("# ", A, G0)') 'dt = ', dt
        WRITE(UNIT=lun, FMT='("# time, stf")')
        DO i=1, SIZE(time)
            WRITE(UNIT=lun, FMT='(G0,", ",G0)') (i-1)*dt, stf(i)
        END DO
        CLOSE(UNIT=lun)

    END SUBROUTINE 

END PROGRAM source_signature_example


!> @file
!! $Date: 2013-10-30 16:25:25 +0100 (Wed, 30 Oct 2013) $
!! $Author: mauerberger $
!! $Revision: 743 $
!! @copyright GNU General Public License version 3 or later 
