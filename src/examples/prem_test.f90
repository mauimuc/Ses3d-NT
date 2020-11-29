!> @example prem_test.f90
!!
!! Example program which tests procedures of the module prem_mod. 
!! The following checks are carried out:
!!  * Testing for Inner Core
!!  * Testing for ICB 
!!  * Testing for CMB
!!  * Testing for D"
!!  * Testing for 670 km transition zone
!!  * Testing for MOHO
!!  * Testing for the Earth's surface
!!
!!  The expected values of the returned parameters according to the passed
!!  radius (depth) you can find on the <a href="http://www.iris.edu/dms/products/emc-prem/">PREM</a> 
!!  
!! Compile by simply using the Makefile:
!!
!!     prem_test
!!
!! Run the program by executing:
!!
!!     ./prem_test
PROGRAM prem_test
    USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY : OUTPUT_UNIT, ERROR_UNIT
    USE parameters_mod, ONLY : real_kind
    USE prem_mod, ONLY : prem
    
    IMPLICIT NONE

REAL(real_kind) :: r, drhodr, rho, vpv, vph, vsv, vsh,  eta, qmu, qkappa
REAL(real_kind) :: rho_, vpv_, vph_, vsv_, vsh_, eta_, qmu_, qkappa_
REAL(real_kind) :: rho_1, rho_2, vpv_1, vpv_2, vph_1, vph_2, vsv_1, vsv_2, vsh_1, vsh_2, eta_1, eta_2

    CHARACTER(LEN=128) :: msg
    LOGICAL :: l

!    LOGICAL :: noocean!, nocrust

    
    ! Testing for the Earth's center (r=0.00 m)
    msg = "Testing subroutine `prem` for the Earth's center (r=0.00):"
    WRITE(UNIT=OUTPUT_UNIT, FMT='(/, A)') TRIM(msg)
    WRITE(UNIT=OUTPUT_UNIT, FMT='(A)') REPEAT('-',LEN_TRIM(msg))

    r = 0.00_real_kind

    rho_ = 13.0885_real_kind
    vpv_ = 11.2622_real_kind
    vph_ = 11.2622_real_kind
    vsv_ = 3.6678_real_kind
    vsh_ = 3.6678_real_kind
    eta_ = 1.00_real_kind
    qkappa_ = 1327.7_real_kind
    qmu_ = 84.6_real_kind

    CALL prem( r, drhodr, rho, vpv, vph, vsv, vsh, eta, qkappa, qmu )
 
    WRITE(UNIT=msg, FMT='(G0," = prem( ",G0,", drhodr, rho, vpv, vph, vsv,&
                         &vsh, eta, qkappa, qmu ) = ",G0)') rho_, r, rho
    WRITE(*,*) "checking for `rho`"
    l = ( rho == rho_ )
    CALL check(l,msg)

    WRITE(UNIT=msg, FMT='(G0," = prem( ",G0,", drhodr, rho, vpv, vph, vsv,&
                         &vsh, eta, qkappa, qmu ) = ",G0)') vpv_, r, vpv
    WRITE(*,*) "checking for `vpv`"
    l = ( vpv == vpv_ )
    CALL check(l,msg)

    WRITE(UNIT=msg, FMT='(G0," = prem( ",G0,", drhodr, rho, vpv, vph, vsv,&
                         &vsh, eta, qkappa, qmu ) = ",G0)') vph_, r, vph
    WRITE(*,*) "checking for `vph`"
    l = ( vph == vph_ )
    CALL check(l,msg)

    WRITE(UNIT=msg, FMT='(G0," = prem( ",G0,", drhodr, rho, vpv, vph, vsv,&
                         &vsh, eta, qkappa, qmu ) = ",G0)') vsv_, r, vsv
    WRITE(*,*) "checking for `vsv`"
    l = ( vsv == vsv_ )
    CALL check(l,msg)

    WRITE(UNIT=msg, FMT='(G0," = prem( ",G0,", drhodr, rho, vpv, vph, vsv,&
                         &vsh, eta, qkappa, qmu ) = ",G0)') vsh_, r, vsh
    WRITE(*,*) "checking for `vsh`"
    l = ( vsh == vsh_ )
    CALL check(l,msg)

    WRITE(UNIT=msg, FMT='(G0," = prem( ",G0,", drhodr, rho, vpv, vph, vsv,&
                         &vsh, eta, qkappa, qmu ) = ",G0)') eta_, r, eta
    WRITE(*,*) "checking for `eta`"
    l = ( eta == eta_ )
    CALL check(l,msg)

    WRITE(UNIT=msg, FMT='(G0," = prem( ",G0,", drhodr, rho, vpv, vph, vsv,&
                         &vsh, eta, qkappa, qmu ) = ",G0)') qkappa_, r, qkappa
    WRITE(*,*) "checking for `qkappa`"
    l = ( qkappa == qkappa_ )
    CALL check(l,msg)

    WRITE(UNIT=msg, FMT='(G0," = prem( ",G0,", drhodr, rho, vpv, vph, vsv,&
                         &vsh, eta, qkappa, qmu ) = ",G0)') qmu_, r, qmu
    WRITE(*,*) "checking for `qmu`"
    l = ( qmu == qmu_ )
    CALL check(l,msg)


    ! Testing 10 m  below CMB (r=3479990.00 m)
    msg = "Testing subroutine `prem` for 10 m above CMB (r=3479990.00 m):"
    WRITE(UNIT=OUTPUT_UNIT, FMT='(/, A)') TRIM(msg)
    WRITE(UNIT=OUTPUT_UNIT, FMT='(A)') REPEAT('-',LEN_TRIM(msg))

    r = 3479990.00_real_kind
   
    rho_1 = 9.90348_real_kind
    rho_2 = 9.90344_real_kind
    vpv_1 = 8.06482_real_kind
    vpv_2 = 8.06479_real_kind
    vph_1 = 8.06482_real_kind
    vph_2 = 8.06479_real_kind
    vsv_1 = 0.00_real_kind
    vsv_2 = 0.00_real_kind
    vsh_1 = 0.00_real_kind
    vsh_2 = 0.00_real_kind
    eta_ = 1.00
    qkappa_ = 57823.00_real_kind
    qmu_ = 0.00_real_kind

    CALL prem( r, drhodr, rho, vpv, vph, vsv, vsh, eta, qkappa, qmu )

    WRITE(UNIT=msg, FMT='("checking if `rho` in the range of [",G0,":",G0,"] &
                          &rho=",G0)') rho_2, rho_1, rho
    l = ( rho_2 <= rho .AND. rho <= rho_1 )
    CALL check(l,msg)

    WRITE(UNIT=msg, FMT='("checking if `vpv` in the range of [",G0,":",G0,"] &
                          &vpv=",G0)') vpv_2, vpv_1, vpv
    l = ( vpv_2 <= vpv .AND. vpv <= vpv_1 )
    CALL check(l,msg)

    WRITE(UNIT=msg, FMT='("checking if `vph` in the range of [",G0,":",G0,"] &
                          &vph=",G0)') vph_2, vph_1, vph
    l = ( vph_2 <= vph .AND. vph <= vph_1 )
    CALL check(l,msg)

    WRITE(UNIT=msg, FMT='("checking for `vsv`: ",G0," = vsv = ",G0)') vsv_2, vsv
    l = ( vsv_2 <= vsv .AND. vsv <= vsv_1 )
    CALL check(l,msg)

    WRITE(UNIT=msg, FMT='("checking for `vsh`: ",G0," = vsh = ",G0)') vsh_2, vsh
    l = ( vsh_2 <= vsh .AND. vsh <= vsh_1 )
    CALL check(l,msg)

    WRITE(UNIT=msg, FMT='("checking for `eta`: ",G0," = eta = ",G0)') eta_, eta
    l = ( eta == eta_ )
    CALL check(l,msg)

    WRITE(UNIT=msg, FMT='("checking for `qkappa`: ",G0," = qkappa = ",G0)') &
                          &qkappa_, qkappa
    l = ( qkappa == qkappa_ )
    CALL check(l,msg)

    WRITE(UNIT=msg, FMT='("checking for `qmu`: ",G0," = qmu  = ",G0)') qmu_, qmu
    l = ( qmu == qmu_ )
    CALL check(l,msg)


    ! Testing 10 m  above CMB (r=3480000.00 m)
    msg = "Testing subroutine `prem` for 10 m above CMB (r=3480010.00 m):"
    WRITE(UNIT=OUTPUT_UNIT, FMT='(/, A)') TRIM(msg)
    WRITE(UNIT=OUTPUT_UNIT, FMT='(A)') REPEAT('-',LEN_TRIM(msg))

    r = 3480010.00_real_kind
   
    rho_1 = 5.56646_real_kind
    rho_2 = 5.56608_real_kind
    vpv_1 = 13.71662_real_kind
    vpv_2 = 13.71570_real_kind
    vph_1 = 13.71662_real_kind
    vph_2 = 13.71570_real_kind
    vsv_1 = 7.26465_real_kind
    vsv_2 = 7.26513_real_kind
    vsh_1 = 7.26465_real_kind
    vsh_2 = 7.26513_real_kind
    eta_ = 1.00
    qkappa_ = 57823.00_real_kind
    qmu_ = 312.00_real_kind

    CALL prem( r, drhodr, rho, vpv, vph, vsv, vsh, eta, qkappa, qmu )

    WRITE(UNIT=msg, FMT='("checking if `rho` in the range of [",G0,":",G0,"] &
                          &rho=",G0)') rho_2, rho_1, rho
    l = ( rho_2 <= rho .AND. rho <= rho_1 )
    CALL check(l,msg)

    WRITE(UNIT=msg, FMT='("checking if `vpv` in the range of [",G0,":",G0,"] &
                          &vpv=",G0)') vpv_2, vpv_1, vpv
    l = ( vpv_2 <= vpv .AND. vpv <= vpv_1 )
    CALL check(l,msg)

    WRITE(UNIT=msg, FMT='("checking if `vph` in the range of [",G0,":",G0,"] &
                          &vph=",G0)') vph_2, vph_1, vph
    l = ( vph_2 <= vph .AND. vph <= vph_1 )
    CALL check(l,msg)

    WRITE(UNIT=msg, FMT='("checking if `vsv` in the range of [",G0,":",G0,"] &
                          &vsv=",G0)') vsv_1, vsv_2, vsv
    l = ( vsv_1 <= vsv .AND. vsv <= vsv_2 )
    CALL check(l,msg)

    WRITE(UNIT=msg, FMT='("checking if `vsh` in the range of [",G0,":",G0,"] &
                          &vsh=",G0)') vsh_1, vsh_2, vsh
    l = ( vsh_1 <= vsh .AND. vsh <= vsh_2 )
    CALL check(l,msg)

    WRITE(UNIT=msg, FMT='("checking for `eta`: ",G0," = eta = ",G0)') eta_, eta
    l = ( eta == eta_ )
    CALL check(l,msg)

    WRITE(UNIT=msg, FMT='("checking for `qkappa`: ",G0," = qkappa = ",G0)') &
                          &qkappa_, qkappa
    l = ( qkappa == qkappa_ )
    CALL check(l,msg)

    WRITE(UNIT=msg, FMT='("checking for `qmu`: ",G0," = qmu  = ",G0)') qmu_, qmu
    l = ( qmu == qmu_ )
    CALL check(l,msg)


    ! Testing 10 m below 670 km transition zone (r=5700990.00 m)
    msg = 'Testing subroutine `prem` for 10 m below 670 km transition zone (r=5700990.00 m):'
    WRITE(UNIT=OUTPUT_UNIT, FMT='(/, A)') TRIM(msg)
    WRITE(UNIT=OUTPUT_UNIT, FMT='(A)') REPEAT('-',LEN_TRIM(msg))

    r = 5700990.00_real_kind

    rho_1 = 4.38117_real_kind
    rho_2 = 4.38074_real_kind
    vpv_1 = 10.75328_real_kind
    vpv_2 = 10.75132_real_kind
    vph_1 = 10.75328_real_kind
    vph_2 = 10.75132_real_kind
    vsv_1 = 5.94571_real_kind
    vsv_2 = 5.94513_real_kind
    vsh_1 = 5.94571_real_kind
    vsh_2 = 5.94513_real_kind
    eta_ = 1.00
    qkappa_ = 57823.00_real_kind
    qmu_ = 312.00_real_kind


    CALL prem( r, drhodr, rho, vpv, vph, vsv, vsh, eta, qkappa, qmu )

    WRITE(UNIT=msg, FMT='("checking if `rho` in the range of [",G0,":",G0,"] &
                          &rho=",G0)') rho_2, rho_1, rho
    l = ( rho_2 <= rho .AND. rho <= rho_1 )
    CALL check(l,msg)

    WRITE(UNIT=msg, FMT='("checking if `vpv` in the range of [",G0,":",G0,"] &
                          &vpv=",G0)') vpv_2, vpv_1, vpv
    l = ( vpv_2 <= vpv .AND. vpv <= vpv_1 )
    CALL check(l,msg)

    WRITE(UNIT=msg, FMT='("checking if `vph` in the range of [",G0,":",G0,"] &
                          &vph=",G0)') vph_2, vph_1, vph
    l = ( vph_2 <= vph .AND. vph <= vph_1 )
    CALL check(l,msg)

    WRITE(UNIT=msg, FMT='("checking if `vsv` in the range of [",G0,":",G0,"] &
                          &vsv=",G0)') vsv_2, vsv_1, vsv
    l = ( vsv_2 <= vsv .AND. vsv <= vsv_1 )
    CALL check(l,msg)

    WRITE(UNIT=msg, FMT='("checking if `vsh` in the range of [",G0,":",G0,"] &
                          &vsh=",G0)') vsh_2, vsh_1, vsh
    l = ( vsh_2 <= vsh .AND. vsh <= vsh_1 )
    CALL check(l,msg)

    WRITE(UNIT=msg, FMT='("checking for `eta`: ",G0," = eta = ",G0)') eta_, eta
    l = ( eta == eta_ )
    CALL check(l,msg)

    WRITE(UNIT=msg, FMT='("checking for `qkappa`: ",G0," = qkappa = ",G0)') &
                          &qkappa_, qkappa
    l = ( qkappa == qkappa_ )
    CALL check(l,msg)

    WRITE(UNIT=msg, FMT='("checking for `qmu`: ",G0," = qmu  = ",G0)') qmu_, qmu
    l = ( qmu == qmu_ )
    CALL check(l,msg)



    ! Testing 10 m above 670 km transition zone (r=5701010.00 m)
    msg = 'Testing subroutine `prem` for 10 m above 670 km transition zone (r=5701010.00 m):'
    WRITE(UNIT=OUTPUT_UNIT, FMT='(/, A)') TRIM(msg)
    WRITE(UNIT=OUTPUT_UNIT, FMT='(A)') REPEAT('-',LEN_TRIM(msg))

    r = 5701010.00_real_kind

    rho_1 = 3.99212_real_kind
    rho_2 = 3.99172_real_kind
    vpv_1 = 10.26617_real_kind
    vpv_2 = 10.26515_real_kind
    vph_1 = 10.26617_real_kind
    vph_2 = 10.26515_real_kind
    vsv_1 = 5.57021_real_kind
    vsv_2 = 5.56965_real_kind
    vsh_1 = 5.57021_real_kind
    vsh_2 = 5.56965_real_kind
    eta_ = 1.00
    qkappa_ = 57823.00_real_kind
    qmu_ = 143.00_real_kind


    CALL prem( r, drhodr, rho, vpv, vph, vsv, vsh, eta, qkappa, qmu )

    WRITE(UNIT=msg, FMT='("checking if `rho` in the range of [",G0,":",G0,"] &
                          &rho=",G0)') rho_2, rho_1, rho
    l = ( rho_2 <= rho .AND. rho <= rho_1 )
    CALL check(l,msg)

    WRITE(UNIT=msg, FMT='("checking if `vpv` in the range of [",G0,":",G0,"] &
                          &vpv=",G0)') vpv_2, vpv_1, vpv
    l = ( vpv_2 <= vpv .AND. vpv <= vpv_1 )
    CALL check(l,msg)

    WRITE(UNIT=msg, FMT='("checking if `vph` in the range of [",G0,":",G0,"] &
                          &vph=",G0)') vph_2, vph_1, vph
    l = ( vph_2 <= vph .AND. vph <= vph_1 )
    CALL check(l,msg)

    WRITE(UNIT=msg, FMT='("checking if `vsv` in the range of [",G0,":",G0,"] &
                          &vsv=",G0)') vsv_2, vsv_1, vsv
    l = ( vsv_2 <= vsv .AND. vsv <= vsv_1 )
    CALL check(l,msg)

    WRITE(UNIT=msg, FMT='("checking if `vsh` in the range of [",G0,":",G0,"] &
                          &vsh=",G0)') vsh_2, vsh_1, vsh
    l = ( vsh_2 <= vsh .AND. vsh <= vsh_1 )
    CALL check(l,msg)

    WRITE(UNIT=msg, FMT='("checking for `eta`: ",G0," = eta = ",G0)') eta_, eta
    l = ( eta == eta_ )
    CALL check(l,msg)

    WRITE(UNIT=msg, FMT='("checking for `qkappa`: ",G0," = qkappa = ",G0)') &
                          &qkappa_, qkappa
    l = ( qkappa == qkappa_ )
    CALL check(l,msg)

    WRITE(UNIT=msg, FMT='("checking for `qmu`: ",G0," = qmu  = ",G0)') qmu_, qmu
    l = ( qmu == qmu_ )
    CALL check(l,msg)


    ! Testing 10 m below MOHO (r=6346590.00 m)
    msg = 'Testing subroutine `prem` for 10 m below MOHO (r=6346590.00 m):'
    WRITE(UNIT=OUTPUT_UNIT, FMT='(/, A)') TRIM(msg)
    WRITE(UNIT=OUTPUT_UNIT, FMT='(A)') REPEAT('-',LEN_TRIM(msg))

    r = 6346590.00_real_kind

    rho_1 = 3.38014_real_kind
    rho_2 = 3.38075_real_kind
    vpv_1 = 8.02061_real_kind
    vpv_2 = 8.02206_real_kind
    vph_1 = 8.18884_real_kind
    vph_2 = 8.19032_real_kind
    vsv_1 = 4.39681_real_kind
    vsv_2 = 4.39602_real_kind
    vsh_1 = 4.61097_real_kind
    vsh_2 = 4.61180_real_kind
    eta_1 = 0.90055_real_kind
    eta_2 = 0.90039_real_kind 
    qkappa_ = 57823.00_real_kind
    qmu_ = 600.00_real_kind


    CALL prem( r, drhodr, rho, vpv, vph, vsv, vsh, eta, qkappa, qmu )

    WRITE(UNIT=msg, FMT='("checking if `rho` in the range of [",G0,":",G0,"] &
                          &rho=",G0)') rho_1, rho_2, rho
    l = ( rho_1 <= rho .AND. rho <= rho_2 )
    CALL check(l,msg)

    WRITE(UNIT=msg, FMT='("checking if `vpv` in the range of [",G0,":",G0,"] &
                          &vpv=",G0)') vpv_1, vpv_2, vpv
    l = ( vpv_1 <= vpv .AND. vpv <= vpv_2 )
    CALL check(l,msg)

    WRITE(UNIT=msg, FMT='("checking if `vph` in the range of [",G0,":",G0,"] &
                          &vph=",G0)') vph_1, vph_2, vph
    l = ( vph_1 <= vph .AND. vph <= vph_2 )
    CALL check(l,msg)

    WRITE(UNIT=msg, FMT='("checking if `vsv` in the range of [",G0,":",G0,"] &
                          &vsv=",G0)') vsv_2, vsv_1, vsv
    l = ( vsv_2 <= vsv .AND. vsv <= vsv_1 )
    CALL check(l,msg)

    WRITE(UNIT=msg, FMT='("checking if `vsh` in the range of [",G0,":",G0,"] &
                          &vsh=",G0)') vsh_1, vsh_2, vsh
    l = ( vsh_1 <= vsh .AND. vsh <= vsh_2 )
    CALL check(l,msg)

    WRITE(UNIT=msg, FMT='("checking if `eta` in the range of [",G0,":",G0,"] &
                          &eta=",G0)') eta_2, eta_1, eta
    l = ( eta_2 <= eta .AND. eta <= eta_1 )
    CALL check(l,msg)

    WRITE(UNIT=msg, FMT='("checking for `qkappa`: ",G0," = qkappa = ",G0)') &
                          &qkappa_, qkappa
    l = ( qkappa == qkappa_ )
    CALL check(l,msg)

    WRITE(UNIT=msg, FMT='("checking for `qmu`: ",G0," = qmu  = ",G0)') qmu_, qmu
    l = ( qmu == qmu_ )
    CALL check(l,msg)


    ! Testing 10 m above MOHO (r=6346610.00 m)
    msg = 'Testing subroutine `prem` for 10 m above MOHO (r=6346610.00 m):'
    WRITE(UNIT=OUTPUT_UNIT, FMT='(/, A)') TRIM(msg)
    WRITE(UNIT=OUTPUT_UNIT, FMT='(A)') REPEAT('-',LEN_TRIM(msg))

    r = 6346610.00_real_kind

    rho_1 = 2.9_real_kind
    rho_2 = 2.9_real_kind
    vpv_1 = 6.8_real_kind
    vpv_2 = 6.8_real_kind
    vph_1 = 6.8_real_kind
    vph_2 = 6.8_real_kind
    vsv_1 = 3.9_real_kind
    vsv_2 = 3.9_real_kind
    vsh_1 = 3.9_real_kind
    vsh_2 = 3.9_real_kind
    eta_ = 1.00
    qkappa_ = 57823.00_real_kind
    qmu_ = 600.00_real_kind


    CALL prem( r, drhodr, rho, vpv, vph, vsv, vsh, eta, qkappa, qmu )

    WRITE(UNIT=msg, FMT='("checking for `rho`: ",G0," = rho = ",G0)') rho_2, rho
    l = ( rho_2 <= rho .AND. rho <= rho_1 )
    CALL check(l,msg)

    WRITE(UNIT=msg, FMT='("checking for `vpv`: ",G0," = vpv = ",G0)') vpv_2, vpv
    l = ( vpv_2 <= vpv .AND. vpv <= vpv_1 )
    CALL check(l,msg)

    WRITE(UNIT=msg, FMT='("checking for `vph`: ",G0," = vph = ",G0)') vph_2, vph
    l = ( vph_2 <= vph .AND. vph <= vph_1 )
    CALL check(l,msg)

    WRITE(UNIT=msg, FMT='("checking for `vsv`: ",G0," = vsv = ",G0)') vsv_2, vsv
    l = ( vsv_2 <= vsv .AND. vsv <= vsv_1 )
    CALL check(l,msg)

    WRITE(UNIT=msg, FMT='("checking for `vsh`: ",G0," = vsh = ",G0)') vsh_2, vsh
    l = ( vsh_2 <= vsh .AND. vsh <= vsh_1 )
    CALL check(l,msg)

    WRITE(UNIT=msg, FMT='("checking for `eta`: ",G0," = eta = ",G0)') eta_, eta
    l = ( eta == eta_ )
    CALL check(l,msg)

    WRITE(UNIT=msg, FMT='("checking for `qkappa`: ",G0," = qkappa = ",G0)') qkappa_, qkappa
    l = ( qkappa == qkappa_ )
    CALL check(l,msg)

    WRITE(UNIT=msg, FMT='("checking for `qmu`: ",G0," = qmu  = ",G0)') qmu_, qmu
    l = ( qmu == qmu_ )
    CALL check(l,msg)

    
    ! Testing for the Earth's surface (noocean = .TRUE.)
    msg = "Testing subroutine `prem` for thr Earth's surface (r=6371000.00 m), nooceam = .TRUE.:"
    WRITE(UNIT=OUTPUT_UNIT, FMT='(/, A)') TRIM(msg)
    WRITE(UNIT=OUTPUT_UNIT, FMT='(A)') REPEAT('-',LEN_TRIM(msg))

    r = 6371000.00_real_kind

    rho_1 = 2.6_real_kind
    rho_2 = 2.6_real_kind
    vpv_1 = 5.8_real_kind
    vpv_2 = 5.8_real_kind
    vph_1 = 5.8_real_kind
    vph_2 = 5.8_real_kind
    vsv_1 = 3.2_real_kind
    vsv_2 = 3.2_real_kind
    vsh_1 = 3.2_real_kind
    vsh_2 = 3.2_real_kind
    eta_ = 1.00
    qkappa_ = 57823.00_real_kind
    qmu_ = 600.00_real_kind

    

    CALL prem( r, drhodr, rho, vpv, vph, vsv, vsh, eta, qkappa, qmu )

    WRITE(UNIT=msg, FMT='("checking for `rho`: ",G0," = rho = ",G0)') rho_2, rho
    l = ( rho_2 <= rho .AND. rho <= rho_1 )
    CALL check(l,msg)

    WRITE(UNIT=msg, FMT='("checking for `vpv`: ",G0," = vpv = ",G0)') vpv_2, vpv
    l = ( vpv_2 <= vpv .AND. vpv <= vpv_1 )
    CALL check(l,msg)

    WRITE(UNIT=msg, FMT='("checking for `vph`: ",G0," = vph = ",G0)') vph_2, vph
    l = ( vph_2 <= vph .AND. vph <= vph_1 )
    CALL check(l,msg)

    WRITE(UNIT=msg, FMT='("checking for `vsv`: ",G0," = vsv = ",G0)') vsv_2, vsv
    l = ( vsv_2 <= vsv .AND. vsv <= vsv_1 )
    CALL check(l,msg)

    WRITE(UNIT=msg, FMT='("checking for `vsh`: ",G0," = vsh = ",G0)') vsh_2, vsh
    l = ( vsh_2 <= vsh .AND. vsh <= vsh_1 )
    CALL check(l,msg)

    WRITE(UNIT=msg, FMT='("checking for `eta`: ",G0," = eta = ",G0)') eta_, eta
    l = ( eta == eta_ )
    CALL check(l,msg)

    WRITE(UNIT=msg, FMT='("checking for `qkappa`: ",G0," = qkappa = ",G0)') qkappa_, qkappa
    l = ( qkappa == qkappa_ )
    CALL check(l,msg)

    WRITE(UNIT=msg, FMT='("checking for `qmu`: ",G0," = qmu  = ",G0)') qmu_, qmu
    l = ( qmu == qmu_ )
    CALL check(l,msg)

 ! Testing for the Earth's surface (noocean = .FALSE.)
!........



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


END PROGRAM prem_test
