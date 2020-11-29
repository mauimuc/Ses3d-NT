!> @example coordinate_utilities_test.f90
!!
!! Example program which tests procedures of the module coordinate_utilities_mod. 
!! The following checks are carried out:
!!  * Testing if depth2radius() and radius2depth() are reverse
!!  * Testing if lat2colat() and colat2lat() are reverse
!!  * Testing if rad2deg() and deg2rad() are reverse
!!  * Testing if xyz2tpr() and tpr2xyz() are reverse
!!  * Testing different cases conversion Cartesian to spherical coordinates
!!  * Testing different cases conversion spherical to Cartesian coordinates 
!!
!! Compile by simply using the Makefile:
!!
!!     make coordinate_utilities_test
!!
!! Run the program by executing:
!!
!!     ./coordinate_utilities_test
PROGRAM coordinate_utilities_test
    USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY : OUTPUT_UNIT, ERROR_UNIT
    USE parameters_mod, ONLY : real_kind, pi
    USE coordinate_utilities_mod, ONLY : deg2rad, rad2deg, &
                                         colat2lat, lat2colat, &
                                         radius2depth, depth2radius, &
                                         xyz2tpr, tpr2xyz

    IMPLICIT NONE 

    REAL(real_kind) :: t, r, p, x, y, z, x_expect, y_expect, z_expect, t_expect, p_expect, r_expect
    CHARACTER(LEN=128) :: msg
    LOGICAL :: l

  

    ! Conversion Cartesian (x,y,z) to spherical coordinates (theta,phi,r)
    WRITE(*,*) 'Different cases convertion Cartesian to spherical coordinates:'


    ! 1. Point in the Earth's interior (x=40.0, y=160.0, z=70.0)
    msg = "Testing if xyz2tpr() and tpr2xyz() are reverse"
    WRITE(UNIT=OUTPUT_UNIT, FMT='(/, A)') TRIM(msg)
    WRITE(UNIT=OUTPUT_UNIT, FMT='(A)') REPEAT('-',LEN_TRIM(msg))

    x_expect = 40.0_real_kind
    y_expect = 160.0_real_kind
    z_expect = 70.0_real_kind

    CALL xyz2tpr(x_expect, y_expect, z_expect, t, p, r)
    CALL tpr2xyz(t, p, r, x, y, z)

    WRITE(UNIT=msg, FMT='(G0," = tpr2xyz(xyz2tpr(",G0,")")') x_expect, x
    WRITE(*,*) "checking for `x`"
    l = ( x == x_expect )
    CALL check(l,msg)

    WRITE(UNIT=msg, FMT='(G0," = tpr2xyz(xyz2tpr(",G0,")")') y_expect, y
    WRITE(*,*) "checking for `y`"
    l = ( y == y_expect )
    CALL check(l,msg)

    WRITE(UNIT=msg, FMT='(G0," = tpr2xyz(xyz2tpr(",G0,")")') z_expect, z
    WRITE(*,*) "checking for `z`"
    l = ( z == z_expect )
    CALL check(l,msg)


    ! 2. Point on the z axis Cartesian coordinates (x=0.0, y=0.0, z=3453.0)
    msg = "Testing conversion Cartesian (x=0.0, y=0.0 ,z=3453.0) to spherical coordinates (theta,phi,r)"
    WRITE(UNIT=OUTPUT_UNIT, FMT='(/, A)') TRIM(msg)
    WRITE(UNIT=OUTPUT_UNIT, FMT='(A)') REPEAT('-',LEN_TRIM(msg))

    CALL xyz2tpr(x=0.0_real_kind, &
                 y=0.0_real_kind, &
                 z=3453.0_real_kind, &
                 theta=t, phi=p, r=r)
    t_expect = 0.0_real_kind
    p_expect = 0.0_real_kind
    r_expect = 3453.0_real_kind

    WRITE(UNIT=msg, FMT='(G0," = xyz2tpr() = ",G0)') t, t_expect
    WRITE(*,*) "checking for `t`"
    l = ( t == t_expect )
    CALL check(l,msg)

    WRITE(UNIT=msg, FMT='(G0," = xyz2tpr() = ",G0)') p, p_expect
    WRITE(*,*) "checking for `p`"
    l = ( p == p_expect )
    CALL check(l,msg)

    WRITE(UNIT=msg, FMT='(G0," = xyz2tpr() = ",G0)') r, r_expect
    WRITE(*,*) "checking for `r`"
    l = ( r == r_expect )
    CALL check(l,msg)


    ! 3. Point in the center of the Earth (x=0.0, y=0.0, z=0.0)
    msg = "Testing conversion Cartesian (x=0.0, y=0.0 ,z=0.0) to spherical coordinates (theta,phi,r)"
    WRITE(UNIT=OUTPUT_UNIT, FMT='(/, A)') TRIM(msg)
    WRITE(UNIT=OUTPUT_UNIT, FMT='(A)') REPEAT('-',LEN_TRIM(msg))

    CALL xyz2tpr(x=0.0_real_kind, &
                 y=0.0_real_kind, &
                 z=0.0_real_kind, &
                 theta=t, phi=p, r=r)
    t_expect = 90.0_real_kind
    p_expect = 0.0_real_kind
    r_expect = 0.0_real_kind


    WRITE(UNIT=msg, FMT='(G0," = xyz2tpr() = ",G0)') t, t_expect
    WRITE(*,*) "checking for `t`"
    l = ( t == t_expect )
    CALL check(l,msg)

    WRITE(UNIT=msg, FMT='(G0," = xyz2tpr() = ",G0)') p, p_expect
    WRITE(*,*) "checking for `p`"
    l = ( p == p_expect )
    CALL check(l,msg)

    WRITE(UNIT=msg, FMT='(G0," = xyz2tpr() = ",G0)') r, r_expect
    WRITE(*,*) "checking for `r`"
    l = ( r == r_expect )
    CALL check(l,msg)
   

    ! Conversion spherical (theta,phi,r) to Cartesian coordinates (x,y,z)
    WRITE(*,*) 'Different cases convertion spherical to Cartesian coordinates'
    WRITE(*,*) 'using functions tpr2xyz(),lat2colat() deg2rad(), depth2rad():'


    ! 1.  Point in the center of the Earth (theta=62.0, phi=35.0, r=0.0)
    msg = "Testing conversion spherical (theta=62.0, phi=35.0, r=0.0) to Cartesian coordinates (x, y, z)"
    WRITE(UNIT=OUTPUT_UNIT, FMT='(/, A)') TRIM(msg)
    WRITE(UNIT=OUTPUT_UNIT, FMT='(A)') REPEAT('-',LEN_TRIM(msg))

    ! 1.1 Conversion theta and phi angles from degree to radians with function deg2rad
    t = deg2rad(deg=62.0_real_kind)
    p = deg2rad(deg=35.0_real_kind)
    ! 1.2 Conversion depth=6371000.0 to radius (Earth's center)
    r = depth2radius(depth=6371000.0_real_kind)
    ! 1.3 Conversion spherical (theta=62.0, phi=35.0, r=0.0) to Cartesian coordinates (x, y, z)
    CALL tpr2xyz(t, p, r, x, y, z)

    x_expect = 0.0_real_kind
    y_expect = 0.0_real_kind
    z_expect = 0.0_real_kind

    WRITE(UNIT=msg, FMT='(G0," = tpr2xyz() = ",G0)') x, x_expect
    WRITE(*,*) "checking for `x`"
    l = ( x == x_expect )
    CALL check(l,msg)

    WRITE(UNIT=msg, FMT='(G0," = tpr2xyz() = ",G0)') y, y_expect
    WRITE(*,*) "checking for `y`"
    l = ( y == y_expect )
    CALL check(l,msg)

    WRITE(UNIT=msg, FMT='(G0," = tpr2xyz() = ",G0)') z, z_expect
    WRITE(*,*) "checking for `z`"
    l = ( z == z_expect )
    CALL check(l,msg)


    ! 2. Point on the equator (theta=90, phi=90, r=6371000.0)
    msg = "Testing conversion spherical (theta=90.0, phi=90.0, r=6371000.0) to Cartesian coordinates (x, y, z)"
    WRITE(UNIT=OUTPUT_UNIT, FMT='(/, A)') TRIM(msg)
    WRITE(UNIT=OUTPUT_UNIT, FMT='(A)') REPEAT('-',LEN_TRIM(msg))
    ! 2.1 Conversion theta and phi angles from degree to radians with function deg2rad
    t = deg2rad(deg=90.0_real_kind)
    p = deg2rad(deg=90.0_real_kind)
    ! 2.2 Conversion depth=0.0 to radius (Earth's surface)
    r = depth2radius(depth=0.0_real_kind)
    ! 2.3 Conversion spherical (theta=90.0, phi=90.0, r=6371000.0) to Cartesian coordinates (x, y, z)
    CALL tpr2xyz(t, p, r, x, y, z)

    x_expect = 0.0_real_kind
    y_expect = 6371000.0_real_kind
    z_expect = 0.0_real_kind

    WRITE(UNIT=msg, FMT='(G0," = tpr2xyz() = ",G0)') x, x_expect
    WRITE(*,*) "checking for `x`"
    l = ( x == x_expect )
    CALL check(l,msg)

    WRITE(UNIT=msg, FMT='(G0," = tpr2xyz() = ",G0)') y, y_expect
    WRITE(*,*) "checking for `y`"
    l = ( y == y_expect )
    CALL check(l,msg)

    WRITE(UNIT=msg, FMT='(G0," = tpr2xyz() = ",G0)') z, z_expect
    WRITE(*,*) "checking for `z`"
    l = ( z == z_expect )
    CALL check(l,msg)

   




    msg = "Testing if depth2radius() and radius2depth() are reverse"
    WRITE(UNIT=OUTPUT_UNIT, FMT='(/, A)') TRIM(msg)
    WRITE(UNIT=OUTPUT_UNIT, FMT='(A)') REPEAT('-',LEN_TRIM(msg))

    r = -1.0
    y = depth2radius( radius2depth(r) )
    WRITE(UNIT=msg, FMT='(G0," = depth2radius(radius2depth(",G0,"))")') r, y
    l = ( r == y )
    CALL check(l, msg)

    r = -0.0
    y = depth2radius( radius2depth(r) )
    WRITE(UNIT=msg, FMT='(G0," = depth2radius(radius2depth(",G0,"))")') r, y
    l = ( r == y )
    CALL check(l, msg)

    r = 1.0
    y = depth2radius( radius2depth(r) )
    WRITE(UNIT=msg, FMT='(G0," = depth2radius(radius2depth(",G0,"))")') r, y
    l = ( r == y )
    CALL check(l, msg)

   
    msg = "Testing if lat2colat() and colat2lat() are reverse"
    WRITE(UNIT=OUTPUT_UNIT, FMT='(/, A)') TRIM(msg)
    WRITE(UNIT=OUTPUT_UNIT, FMT='(A)') REPEAT('-',LEN_TRIM(msg))

    r = -1.0
    y = lat2colat( colat2lat(r) )
    WRITE(UNIT=msg, FMT='(G0," = lat2colat(colat2lat(",G0,"))")') r, y
    l = ( r == y )
    CALL check(l, msg)

    r = -0.0
    y = lat2colat( colat2lat(r) )
    WRITE(UNIT=msg, FMT='(G0," = lat2colat(colat2lat(",G0,"))")') r, y
    l = ( r == y )
    CALL check(l, msg)

    r = 1.0
    y = lat2colat( colat2lat(r) )
    WRITE(UNIT=msg, FMT='(G0," = lat2colat(colat2lat(",G0,"))")') r, y
    l = ( r == y )
    CALL check(l, msg)

   
    msg = "Testing if rad2deg() and deg2rad() are reverse"
    WRITE(UNIT=OUTPUT_UNIT, FMT='(/, A)') TRIM(msg)
    WRITE(UNIT=OUTPUT_UNIT, FMT='(A)') REPEAT('-',LEN_TRIM(msg))

    r = -1.0
    y = rad2deg( deg2rad(r) )
    WRITE(UNIT=msg, FMT='(G0," = rad2deg(deg2rad(",G0,"))")') r, y
    l = ( r == y )
    CALL check(l, msg)

    r = -0.0
    y = rad2deg( deg2rad(r) )
    WRITE(UNIT=msg, FMT='(G0," = rad2deg(deg2rad(",G0,"))")') r, y
    l = ( r == y )
    CALL check(l, msg)

    r = 1.0
    y = rad2deg( deg2rad(r) )
    WRITE(UNIT=msg, FMT='(G0," = rad2deg(deg2rad(",G0,"))")') r, y
    l = ( r == y )
    CALL check(l, msg)

    WRITE(UNIT=OUTPUT_UNIT, FMT='()') 


    !d = 6371000.1 
    ! Well, one difficulty of depth2radius() is:
    !  * Depth > 6371000 
    ! What can we do about it?
    !  * A solution might be to take the absolute value 


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

END PROGRAM coordinate_utilities_test


!> @file
!! $Date: 2013-10-30 16:25:25 +0100 (Wed, 30 Oct 2013) $
!! $Author: mauerberger $
!! $Revision: 743 $
!! @copyright GNU General Public License version 3 or later
