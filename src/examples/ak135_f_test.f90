!> @example ak135_f_test.f90
!!
!! Example program which tests procedures of the module ak135_f_mod. 
!! The following checks are carried out:
!!  * Testing for the Earth's surface
!!  * Testing for Lehmann discontinuity
!!  * Testing for 410 km transition zone
!!  * Testing for 660 km transition zone
!!  * Testing for CMB
!!  * Testing for ICB
!!
!!  The expected values of the returned parameters according to the passed
!!  radius (depth) you can find on the <a href="http://rses.anu.edu.au/seismology/ak135/ak135f.html">ak135_f model homepage</a> 
!!  
!! Compile by simply using the Makefile:
!!
!!     make ak135_f_test
!!
!! Run the program by executing:
!!
!!     ./ak135_f_test
PROGRAM ak135_f_test
    USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY : OUTPUT_UNIT, ERROR_UNIT
    USE parameters_mod, ONLY : real_kind
    USE ak135_f_mod, ONLY : ak135_f
    
    IMPLICIT NONE

    REAL(real_kind) :: x, rho, vp, vs, qmu, qkappa, rho_expect, vp_expect, vs_expect, qmu_expect, qkappa_expect
    CHARACTER(LEN=128) :: msg
    LOGICAL :: l


    ! Testing for the Earth's surface x=6371.00
    msg = "Testing subroutine 'ak135_f()' for the Earth's surface (radius x=6371.00):"
    WRITE(UNIT=OUTPUT_UNIT, FMT='(/, A)') TRIM(msg)
    WRITE(UNIT=OUTPUT_UNIT, FMT='(A)') REPEAT('-',LEN_TRIM(msg))

    x = 6371.00_real_kind

    rho_expect = 1.02_real_kind
    vp_expect = 1.45_real_kind
    vs_expect = 0.00_real_kind
    qmu_expect = 0.00_real_kind
    qkappa_expect = 57822.00_real_kind
    CALL ak135_f( x, rho, vp, vs, qmu, qkappa )

    WRITE(UNIT=msg, FMT='(G0," = ak135_f( ",G0,", rho, vp, vs, qmu, qkappa ) = ",G0 )') rho_expect, x, rho
    WRITE(*,*) "checking for `rho`"
    l = ( rho == rho_expect )
    CALL check(l,msg)

    WRITE(UNIT=msg, FMT='(G0," = ak135_f( ",G0,", rho, vp, vs, qmu, qkappa ) = ",G0 )') vp_expect, x, vp
    WRITE(*,*) "checking for `vp`"
    l = ( vp == vp_expect )
    CALL check(l,msg)

    WRITE(UNIT=msg, FMT='(G0," = ak135_f( ",G0,", rho, vp, vs, qmu, qkappa ) = ",G0 )') vs_expect, x, vs
    WRITE(*,*) "checking for `vs`"
    l = ( vs == vs_expect )
    CALL check(l,msg)

    WRITE(UNIT=msg, FMT='(G0," = ak135_f( ",G0,", rho, vp, vs, qmu, qkappa ) = ",G0 )') qmu_expect, x, qmu
    WRITE(*,*) "checking for `qmu`"
    l = ( qmu == qmu_expect )
    CALL check(l,msg)

    WRITE(UNIT=msg, FMT='(G0," = ak135_f( ",G0,", rho, vp, vs, qmu, qkappa ) = ",G0 )') qkappa_expect, x, qkappa
    WRITE(*,*) "checking for `qkappa`"
    l = ( qkappa == qkappa_expect )
    CALL check(l,msg)



    ! Testing for Lehmann discontinuity depth=210 km => x=6161.00 km
    msg = "Testing subroutine 'ak135_f()' for Lehmann discontinuity (radius x=6161.00):"
    WRITE(UNIT=OUTPUT_UNIT, FMT='(/, A)') TRIM(msg)
    WRITE(UNIT=OUTPUT_UNIT, FMT='(A)') REPEAT('-',LEN_TRIM(msg))

    x =  6161.00_real_kind
   
    rho_expect = 3.3243_real_kind
    vp_expect = 8.3007_real_kind
    vs_expect = 4.5184_real_kind
    qmu_expect = 133.72_real_kind
    qkappa_expect = 338.47_real_kind
    CALL ak135_f( x, rho, vp, vs, qmu, qkappa )

    WRITE(UNIT=msg, FMT='(G0," = ak135_f( ",G0,", rho, vp, vs, qmu, qkappa ) = ",G0 )') rho_expect, x, rho
    WRITE(*,*) "checking for `rho`"
    l = ( rho == rho_expect )
    CALL check(l,msg)

    WRITE(UNIT=msg, FMT='(G0," = ak135_f( ",G0,", rho, vp, vs, qmu, qkappa ) = ",G0 )') vp_expect, x, vp
    WRITE(*,*) "checking for `vp`"
    l = ( vp == vp_expect )
    CALL check(l,msg)

    WRITE(UNIT=msg, FMT='(G0," = ak135_f( ",G0,", rho, vp, vs, qmu, qkappa ) = ",G0 )') vs_expect, x, vs
    WRITE(*,*) "checking for `vs`"
    l = ( vs == vs_expect )
    CALL check(l,msg)

    WRITE(UNIT=msg, FMT='(G0," = ak135_f( ",G0,", rho, vp, vs, qmu, qkappa ) = ",G0 )') qmu_expect, x, qmu
    WRITE(*,*) "checking for `qmu`"
    l = ( qmu == qmu_expect )
    CALL check(l,msg)

    WRITE(UNIT=msg, FMT='(G0," = ak135_f( ",G0,", rho, vp, vs, qmu, qkappa ) = ",G0 )') qkappa_expect, x, qkappa
    WRITE(*,*) "checking for `qkappa`"
    l = ( qkappa == qkappa_expect )
    CALL check(l,msg)



    ! Testing for 410 km transition zone x=5961.00 km
    msg = "Testing subroutine 'ak135_f()' for 410 km transition zone (radius x=5961.00):"
    WRITE(UNIT=OUTPUT_UNIT, FMT='(/, A)') TRIM(msg)
    WRITE(UNIT=OUTPUT_UNIT, FMT='(A)') REPEAT('-',LEN_TRIM(msg))

    x =  5961.00_real_kind
   
    rho_expect = 3.9317_real_kind
    vp_expect = 9.3601_real_kind
    vs_expect = 5.0806_real_kind
    qmu_expect = 162.50_real_kind
    qkappa_expect = 413.66_real_kind
    CALL ak135_f( x, rho, vp, vs, qmu, qkappa )

    WRITE(UNIT=msg, FMT='(G0," = ak135_f( ",G0,", rho, vp, vs, qmu, qkappa ) = ",G0 )') rho_expect, x, rho
    WRITE(*,*) "checking for `rho`"
    l = ( rho == rho_expect )
    CALL check(l,msg)

    WRITE(UNIT=msg, FMT='(G0," = ak135_f( ",G0,", rho, vp, vs, qmu, qkappa ) = ",G0 )') vp_expect, x, vp
    WRITE(*,*) "checking for `vp`"
    l = ( vp == vp_expect )
    CALL check(l,msg)

    WRITE(UNIT=msg, FMT='(G0," = ak135_f( ",G0,", rho, vp, vs, qmu, qkappa ) = ",G0 )') vs_expect, x, vs
    WRITE(*,*) "checking for `vs`"
    l = ( vs == vs_expect )
    CALL check(l,msg)

    WRITE(UNIT=msg, FMT='(G0," = ak135_f( ",G0,", rho, vp, vs, qmu, qkappa ) = ",G0 )') qmu_expect, x, qmu
    WRITE(*,*) "checking for `qmu`"
    l = ( qmu == qmu_expect )
    CALL check(l,msg)

    WRITE(UNIT=msg, FMT='(G0," = ak135_f( ",G0,", rho, vp, vs, qmu, qkappa ) = ",G0 )') qkappa_expect, x, qkappa
    WRITE(*,*) "checking for `qkappa`"
    l = ( qkappa == qkappa_expect )
    CALL check(l,msg)



    ! Testing for 660 km transition zone x=5711.00 km
    msg = "Testing subroutine 'ak135_f()' for 660 km transition zone (radius x=5711.00):"
    WRITE(UNIT=OUTPUT_UNIT, FMT='(/, A)') TRIM(msg)
    WRITE(UNIT=OUTPUT_UNIT, FMT='(A)') REPEAT('-',LEN_TRIM(msg))

    x =  5711.00_real_kind
   
    rho_expect = 4.2387_real_kind
    vp_expect = 10.7909_real_kind
    vs_expect = 5.9607_real_kind
    qmu_expect = 549.45_real_kind
    qkappa_expect = 1350.54_real_kind
    CALL ak135_f( x, rho, vp, vs, qmu, qkappa )

    WRITE(UNIT=msg, FMT='(G0," = ak135_f( ",G0,", rho, vp, vs, qmu, qkappa ) = ",G0 )') rho_expect, x, rho
    WRITE(*,*) "checking for `rho`"
    l = ( rho == rho_expect )
    CALL check(l,msg)

    WRITE(UNIT=msg, FMT='(G0," = ak135_f( ",G0,", rho, vp, vs, qmu, qkappa ) = ",G0 )') vp_expect, x, vp
    WRITE(*,*) "checking for `vp`"
    l = ( vp == vp_expect )
    CALL check(l,msg)

    WRITE(UNIT=msg, FMT='(G0," = ak135_f( ",G0,", rho, vp, vs, qmu, qkappa ) = ",G0 )') vs_expect, x, vs
    WRITE(*,*) "checking for `vs`"
    l = ( vs == vs_expect )
    CALL check(l,msg)

    WRITE(UNIT=msg, FMT='(G0," = ak135_f( ",G0,", rho, vp, vs, qmu, qkappa ) = ",G0 )') qmu_expect, x, qmu
    WRITE(*,*) "checking for `qmu`"
    l = ( qmu == qmu_expect )
    CALL check(l,msg)

    WRITE(UNIT=msg, FMT='(G0," = ak135_f( ",G0,", rho, vp, vs, qmu, qkappa ) = ",G0 )') qkappa_expect, x, qkappa
    WRITE(*,*) "checking for `qkappa`"
    l = ( qkappa == qkappa_expect )
    CALL check(l,msg)



    ! Testing for CMB depth=2891.50 km => x=3479.50 km
    msg = "Testing subroutine 'ak135_f()' for CMB (radius x=3479.50):"
    WRITE(UNIT=OUTPUT_UNIT, FMT='(/, A)') TRIM(msg)
    WRITE(UNIT=OUTPUT_UNIT, FMT='(A)') REPEAT('-',LEN_TRIM(msg))

    x =  3479.50_real_kind
   
    rho_expect = 9.9145_real_kind
    vp_expect = 8.00_real_kind
    vs_expect = 0.00_real_kind
    qmu_expect = 0.00_real_kind
    qkappa_expect = 57822.00_real_kind
    CALL ak135_f( x, rho, vp, vs, qmu, qkappa )

    WRITE(UNIT=msg, FMT='(G0," = ak135_f( ",G0,", rho, vp, vs, qmu, qkappa ) = ",G0 )') rho_expect, x, rho
    WRITE(*,*) "checking for `rho`"
    l = ( rho == rho_expect )
    CALL check(l,msg)

    WRITE(UNIT=msg, FMT='(G0," = ak135_f( ",G0,", rho, vp, vs, qmu, qkappa ) = ",G0 )') vp_expect, x, vp
    WRITE(*,*) "checking for `vp`"
    l = ( vp == vp_expect )
    CALL check(l,msg)

    WRITE(UNIT=msg, FMT='(G0," = ak135_f( ",G0,", rho, vp, vs, qmu, qkappa ) = ",G0 )') vs_expect, x, vs
    WRITE(*,*) "checking for `vs`"
    l = ( vs == vs_expect )
    CALL check(l,msg)

    WRITE(UNIT=msg, FMT='(G0," = ak135_f( ",G0,", rho, vp, vs, qmu, qkappa ) = ",G0 )') qmu_expect, x, qmu
    WRITE(*,*) "checking for `qmu`"
    l = ( qmu == qmu_expect )
    CALL check(l,msg)

    WRITE(UNIT=msg, FMT='(G0," = ak135_f( ",G0,", rho, vp, vs, qmu, qkappa ) = ",G0 )') qkappa_expect, x, qkappa
    WRITE(*,*) "checking for `qkappa`"
    l = ( qkappa == qkappa_expect )
    CALL check(l,msg)



    ! Testing for Lehmann–Bullen discontinuity (inner core boundary) depth=5153.50 km => x=1217.50 km
    msg = "Testing subroutine 'ak135_f()' for ICB (radius x=1217.50):"
    WRITE(UNIT=OUTPUT_UNIT, FMT='(/, A)') TRIM(msg)
    WRITE(UNIT=OUTPUT_UNIT, FMT='(A)') REPEAT('-',LEN_TRIM(msg))

    x =  1217.50_real_kind
   
    rho_expect = 12.7037_real_kind
    vp_expect = 11.0427_real_kind
    vs_expect = 3.5043_real_kind
    qmu_expect = 85.03_real_kind
    qkappa_expect = 633.26_real_kind
    CALL ak135_f( x, rho, vp, vs, qmu, qkappa )

    WRITE(UNIT=msg, FMT='(G0," = ak135_f( ",G0,", rho, vp, vs, qmu, qkappa ) = ",G0 )') rho_expect, x, rho
    WRITE(*,*) "checking for `rho`"
    l = ( rho == rho_expect )
    CALL check(l,msg)

    WRITE(UNIT=msg, FMT='(G0," = ak135_f( ",G0,", rho, vp, vs, qmu, qkappa ) = ",G0 )') vp_expect, x, vp
    WRITE(*,*) "checking for `vp`"
    l = ( vp == vp_expect )
    CALL check(l,msg)

    WRITE(UNIT=msg, FMT='(G0," = ak135_f( ",G0,", rho, vp, vs, qmu, qkappa ) = ",G0 )') vs_expect, x, vs
    WRITE(*,*) "checking for `vs`"
    l = ( vs == vs_expect )
    CALL check(l,msg)

    WRITE(UNIT=msg, FMT='(G0," = ak135_f( ",G0,", rho, vp, vs, qmu, qkappa ) = ",G0 )') qmu_expect, x, qmu
    WRITE(*,*) "checking for `qmu`"
    l = ( qmu == qmu_expect )
    CALL check(l,msg)

    WRITE(UNIT=msg, FMT='(G0," = ak135_f( ",G0,", rho, vp, vs, qmu, qkappa ) = ",G0 )') qkappa_expect, x, qkappa
    WRITE(*,*) "checking for `qkappa`"
    l = ( qkappa == qkappa_expect )
    CALL check(l,msg)




CONTAINS
    
    SUBROUTINE check(l,s)
        LOGICAL :: l
        CHARACTER(LEN=*) :: s

        IF ( l ) THEN
            WRITE(UNIT=OUTPUT_UNIT, FMT='(A,A)') '[PASS] ', TRIM(s)
        ELSE
            WRITE(UNIT=ERROR_UNIT, FMT='(A,A)') '[FAIL] ', TRIM(s)
  !          STOP 1
        END IF

    END SUBROUTINE


END PROGRAM ak135_f_test
