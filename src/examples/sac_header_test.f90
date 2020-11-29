!> @example sac_header_test.f90
!!
!! Example program which tests procedures of sac_header_mod
!!
!! @todo Describe what it does ...
!!
!! Compile simply useing the Makefile:
!!
!!     make string_utilities_test
!!
!! ...
PROGRAM sac_header_test
    USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY : OUTPUT_UNIT, ERROR_UNIT
    USE sac_header_mod, ONLY : enum, enumerated_header_variables, &
                               header, header_variables, &
                               rk, header_cls, iday

    IMPLICIT NONE
    
    INTEGER :: n
    CHARACTER(LEN=128) :: msg

    INTEGER ::   i0(4), i(4)
    REAL(rk)     r0(4), r(4)
    LOGICAL ::   l0(4), l(4)
    CHARACTER(LEN=8)  :: c0(4), c(4)

    CHARACTER(LEN=13), PARAMETER :: OoB = 'OUT OF BOUNDS'
    CHARACTER(LEN=LEN(OoB)) :: hdr_var
    CHARACTER(LEN=LEN(OoB)), ALLOCATABLE :: vals(:)
    INTEGER, ALLOCATABLE :: words(:)
   

    ! Storage size
    WRITE(UNIT=OUTPUT_UNIT, FMT='(/,A)') 'Testing storage size'

    WRITE(msg, "('Storage size of logical(', I0, ') is ', I0, '-bit!')") &
        KIND(.true.), STORAGE_SIZE(.true.)
    IF ( STORAGE_SIZE(.true.) /= 32 ) THEN
        WRITE(UNIT=ERROR_UNIT, FMT=*) '[FAIL] ', TRIM(msg), &
            ' SAC expects 32-bit.'
        STOP 
    ELSE
        WRITE(UNIT=OUTPUT_UNIT, FMT=*) '[PASS] ', TRIM(msg)
    END IF




    ! Test enum()
    WRITE(UNIT=OUTPUT_UNIT, FMT='(/,A)') 'Testing function enum()'


    ! 1st Test
    n=1
    DO 
        hdr_var = enum(n)
        IF ( hdr_var==OoB ) &
            EXIT
        n = n + 1
    END DO
    n = n - 1 ! We started counting from 1 so we have to subtract 1

    IF ( n /= SIZE(enumerated_header_variables) ) THEN
        STOP "Counting SAC's enumerated header variables!"
    ELSE
        WRITE(UNIT=OUTPUT_UNIT, FMT=*) "[PASS] Counting SAC's enumerated &
            &header variables!"
    END IF


    ! Check bounds 
    IF ( enum(0) /= OoB .OR. enum(n+1) /= OoB ) THEN
        STOP "Upper or lower bound of enum() invalid!"
    ELSE 
        WRITE(UNIT=OUTPUT_UNIT, FMT=*) "[PASS] Checking enum()'s upper and &
            &lower bounds!"
    END IF


    ! Check for unknown variables
    IF ( enum('') /= -1 .OR. enum(' ') /= -1 .OR. &
         enum('123456789') /= -1 ) THEN
        STOP 'enum() could not detect unknown variables!'
    ELSE
        WRITE(UNIT=OUTPUT_UNIT, FMT=*) '[PASS] enum() properly detected &
            &unknown variables!'
    END IF


    ! 2ed Test
    ALLOCATE( vals(n) )

    DO n=1, SIZE(vals)
        vals(n) = enum(n)
    END DO

    IF ( .NOT. ALL(vals==enumerated_header_variables) ) THEN
        STOP 'enum() returned unknown variable!'
    ELSE
        WRITE(UNIT=OUTPUT_UNIT, FMT=*) '[PASS] All varibales returned by &
            &enum() are valid!'
    END IF
        

    ! 3rd Test 
    IF ( .NOT. ALL(enum(vals) == &
                   [(n, n=1, SIZE(enumerated_header_variables))] ) ) THEN
        STOP 'enum() returned wrong position! '
    ELSE
        WRITE(UNIT=OUTPUT_UNIT, FMT=*) '[PASS] All positions returned by &
            &enum() are matching!'
    END IF


    

    ! Test header()
    WRITE(UNIT=OUTPUT_UNIT, FMT='(/,A)') 'Testing function header()'
        
    ! Test keywords UNUSED and INTERNAL
    IF ( header('UNUSED') /= -2 .OR. header('INTERNAL') /= -3 ) THEN
        STOP 'header() could not detect UNUSED or INTERNAL!'
    ELSE
        WRITE(UNIT=OUTPUT_UNIT, FMT=*) '[PASS] Header keywords UNUSED and &
            &INTERNAL properly detected!'
    END IF

    ! Test for unknown keywords 
    IF ( header('') /= -1 .OR. header(' ') /= -1 .OR. &
         header('123456789') /= -1 ) THEN
        STOP 'header() could not detect unknown variables!'
    ELSE
        WRITE(UNIT=OUTPUT_UNIT, FMT=*) '[PASS] header() properly detected &
            &unknown variables!'
    END IF

    ! Check bounds 
    IF ( header(0) /= OoB .OR. header(159) /= OoB ) THEN
        STOP "Upper or lower bound of header() invalid!"
    ELSE 
        WRITE(UNIT=OUTPUT_UNIT, FMT=*) "[PASS] Checking header()'s upper &
            &and lower bounds!"
    END IF

    ! Check returned header keywords
    DEALLOCATE(vals) 
    ALLOCATE(vals(158))

    DO n=1, 158
        vals(n) = header(n)
    END DO

    IF ( .NOT. ALL(vals==header_variables) ) THEN
        STOP 'header() returned unknown variable!'
    ELSE
        WRITE(UNIT=OUTPUT_UNIT, FMT=*) '[PASS] All varibales returned by &
            &header() are valid!'
    END IF
        
    ! Check if returned position are matching
    words = [(n,n=1,158)]
    WHERE ( header_variables == 'INTERNAL' )
        words = -3
    ELSE WHERE ( header_variables == 'UNUSED' )
        words = -2
    END WHERE

    IF ( .NOT. ALL(header(vals) == words)) THEN
        STOP 'header() returned wrong position! '
    ELSE
        WRITE(UNIT=OUTPUT_UNIT, FMT=*) '[PASS] All positions returned by &
            &header() are matching!'
    END IF


    ! Tesing class header_cls()
    WRITE(UNIT=OUTPUT_UNIT, FMT='(/,A)') 'Tesing header_cls'

    ! Check piping reals through 
    r0(:) = [TINY(1.0),0.0,1.0,HUGE(1.0)]

    r(:) = header_cls(r0)
    r(:) = header_cls(r)
    WRITE(msg, "('Piping real(', I0, ') throug real(', I0, ')!')") &
        KIND(r0), KIND(r)
    IF ( .NOT. ALL( r == r0) ) THEN
        WRITE(UNIT=ERROR_UNIT, FMT=*) '[FAIL] ', TRIM(msg)
        STOP 
    ELSE
        WRITE(UNIT=OUTPUT_UNIT, FMT=*) '[PASS] ', TRIM(msg)
    END IF

    i(:) = header_cls(r0)
    r(:) = header_cls(i)
    WRITE(msg, "('Piping real(', I0, ') throug integer(', I0, ')!')") &
        KIND(r0), KIND(i)
    IF ( STORAGE_SIZE(i)<STORAGE_SIZE(r0) ) THEN
        WRITE(UNIT=OUTPUT_UNIT, FMT=*) '[SKIP] ', TRIM(msg)
    ELSE IF ( .NOT. ALL( r == r0) ) THEN
        WRITE(UNIT=OUTPUT_UNIT, FMT=*) '[FAIL] ', TRIM(msg)
        STOP 
    ELSE
        WRITE(UNIT=OUTPUT_UNIT, FMT=*) '[PASS] ', TRIM(msg)
    END IF

    l(:) = header_cls(r0)
    r(:) = header_cls(l)
    WRITE(msg, "('Piping real(', I0, ') throug logical(', I0, ')!')") &
        KIND(r0), KIND(l)
    IF ( STORAGE_SIZE(l)<STORAGE_SIZE(r0) ) THEN
        WRITE(UNIT=OUTPUT_UNIT, FMT=*) '[SKIP] ', TRIM(msg)
    ELSE IF ( .NOT. ALL( r == r0) ) THEN
        WRITE(UNIT=OUTPUT_UNIT, FMT=*) '[FAIL] ', TRIM(msg)
        STOP 
    ELSE
        WRITE(UNIT=OUTPUT_UNIT, FMT=*) '[PASS] ', TRIM(msg)
    END IF

    c(:) = header_cls(r0)
    r(:) = header_cls(c)
    WRITE(msg, "('Piping real(', I0, ') throug character(', I0, ')!')") &
        KIND(r0), LEN(c)
    IF ( STORAGE_SIZE(c)<STORAGE_SIZE(r0) ) THEN
        WRITE(UNIT=OUTPUT_UNIT, FMT=*) '[SKIP] ', TRIM(msg)
    ELSE IF ( .NOT. ALL( r == r0) ) THEN
        WRITE(UNIT=OUTPUT_UNIT, FMT=*) '[FAIL] ', TRIM(msg)
        STOP 
    ELSE
        WRITE(UNIT=OUTPUT_UNIT, FMT=*) '[PASS] ', TRIM(msg)
    END IF


    ! Check piping logicals through 
    l0(:) = [.true.,.false.,.false.,.true.]

    r(:) = header_cls(l0)
    l(:) = header_cls(r)
    WRITE(msg, "('Piping logical(', I0, ') throug real(', I0, ')!')") &
        KIND(l0), KIND(r)
    IF ( STORAGE_SIZE(r)<STORAGE_SIZE(l0) ) THEN
        WRITE(UNIT=OUTPUT_UNIT, FMT=*) '[SKIP] ', TRIM(msg)
    ELSE IF ( .NOT. ALL( l .EQV. l0) ) THEN
        WRITE(UNIT=OUTPUT_UNIT, FMT=*) '[FAIL] ', TRIM(msg)
        STOP 
    ELSE
        WRITE(UNIT=OUTPUT_UNIT, FMT=*) '[PASS] ', TRIM(msg)
    END IF

    i(:) = header_cls(l0)
    l(:) = header_cls(i)
    WRITE(msg, "('Piping logical(', I0, ') throug integer(', I0, ')!')") &
        KIND(l0), KIND(i)
    IF ( STORAGE_SIZE(i)<STORAGE_SIZE(l0) ) THEN
        WRITE(UNIT=OUTPUT_UNIT, FMT=*) '[SKIP] ', TRIM(msg)
    ELSE IF ( .NOT. ALL( l .EQV. l0) ) THEN
        WRITE(UNIT=OUTPUT_UNIT, FMT=*) '[FAIL] ', TRIM(msg)
        STOP 
    ELSE
        WRITE(UNIT=OUTPUT_UNIT, FMT=*) '[PASS] ', TRIM(msg)
    END IF

    l = header_cls(l0)
    l(:) = header_cls(l)
    WRITE(msg, "('Piping logical(', I0, ') throug logical(', I0, ')!')") &
        KIND(l0), KIND(l)
    IF ( .NOT. ALL( l .EQV. l0) ) THEN
        WRITE(UNIT=OUTPUT_UNIT, FMT=*) '[FAIL] ', TRIM(msg)
        STOP 
    ELSE
        WRITE(UNIT=OUTPUT_UNIT, FMT=*) '[PASS] ', TRIM(msg)
    END IF

    c(:) = header_cls(l0)
    l(:) = header_cls(c)
    WRITE(msg, "('Piping logical(', I0, ') throug character(', I0, ')!')") &
        KIND(l0), LEN(c)
    IF ( STORAGE_SIZE(c)<STORAGE_SIZE(l0) ) THEN
        WRITE(UNIT=OUTPUT_UNIT, FMT=*) '[SKIP] ', TRIM(msg)
    ELSE IF ( .NOT. ALL( l .EQV. l0) ) THEN
        WRITE(UNIT=OUTPUT_UNIT, FMT=*) '[FAIL] ', TRIM(msg)
        STOP 
    ELSE
        WRITE(UNIT=OUTPUT_UNIT, FMT=*) '[PASS] ', TRIM(msg)
    END IF


    ! Check piping characters through
    c0(:) = ['12345678','90abcdef','ghijklmn','opqrstuv' ]

    r(:) = header_cls(c0)
    c(:) = header_cls(r)
    WRITE(msg, "('Piping character(', I0, ') throug real(', I0, ')!')") &
        LEN(c0), KIND(r)
    IF ( STORAGE_SIZE(r)<STORAGE_SIZE(c0) ) THEN
        WRITE(UNIT=OUTPUT_UNIT, FMT=*) '[SKIP] ', TRIM(msg)
    ELSE IF ( .NOT. ALL( c == c0) ) THEN
        WRITE(UNIT=OUTPUT_UNIT, FMT=*) '[FAIL] ', TRIM(msg)
        STOP 
    ELSE
        WRITE(UNIT=OUTPUT_UNIT, FMT=*) '[PASS] ', TRIM(msg)
    END IF

    i(:) = header_cls(c0)
    c(:) = header_cls(i)
    WRITE(msg, "('Piping character(', I0, ') throug integer(', I0, ')!')") &
        LEN(c0), KIND(i)
    IF ( STORAGE_SIZE(i)<STORAGE_SIZE(c0) ) THEN
        WRITE(UNIT=OUTPUT_UNIT, FMT=*) '[SKIP] ', TRIM(msg)
    ELSE IF ( .NOT. ALL( c == c0) ) THEN
        WRITE(UNIT=OUTPUT_UNIT, FMT=*) '[FAIL] ', TRIM(msg)
        STOP 
    ELSE
        WRITE(UNIT=OUTPUT_UNIT, FMT=*) '[PASS] ', TRIM(msg)
    END IF

    l(:) = header_cls(c0)
    c(:) = header_cls(l)
    WRITE(msg, "('Piping character(', I0, ') throug logical(', I0, ')!')") &
        LEN(c0), KIND(l)
    IF ( STORAGE_SIZE(l)<STORAGE_SIZE(c0) ) THEN
        WRITE(UNIT=OUTPUT_UNIT, FMT=*) '[SKIP] ', TRIM(msg)
    ELSE IF ( .NOT. ALL( c == c0) ) THEN
        WRITE(UNIT=OUTPUT_UNIT, FMT=*) '[FAIL] ', TRIM(msg)
        STOP 
    ELSE
        WRITE(UNIT=OUTPUT_UNIT, FMT=*) '[PASS] ', TRIM(msg)
    END IF

    c(:) = header_cls(c0)
    c(:) = header_cls(c)
    WRITE(msg, "('Piping character(', I0, ') throug character(', I0, ')!')") &
        LEN(c0), LEN(c)
    IF ( .NOT. ALL( c == c0) ) THEN
        WRITE(UNIT=OUTPUT_UNIT, FMT=*) '[FAIL] ', TRIM(msg)
        STOP 
    ELSE
        WRITE(UNIT=OUTPUT_UNIT, FMT=*) '[PASS] ', TRIM(msg)
    END IF


    ! Check piping a integer through 
    i0(:) = [-HUGE(1), 0, 1, HUGE(1)]

    r(:) = header_cls(i0)
    i(:) = header_cls(r)
    WRITE(msg, "('Piping integer(', I0, ') throug real(', I0, ')!')") &
        KIND(i0), KIND(r)
    IF ( STORAGE_SIZE(r)<STORAGE_SIZE(i0) ) THEN
        WRITE(UNIT=OUTPUT_UNIT, FMT=*) '[SKIP] ', TRIM(msg)
    ELSE IF ( .NOT. ALL( i == i0) ) THEN
        WRITE(UNIT=ERROR_UNIT, FMT=*) '[FAIL] ', TRIM(msg)
        STOP 
    ELSE
        WRITE(UNIT=OUTPUT_UNIT, FMT=*) '[PASS] ', TRIM(msg)
    END IF

    i(:) = header_cls(i0)
    i(:) = header_cls(i)
    WRITE(msg, "('Piping integer(', I0, ') throug integer(', I0, ')!')") &
        KIND(i0), KIND(i)
    IF ( .NOT. ALL( i == i0) ) THEN
        WRITE(UNIT=OUTPUT_UNIT, FMT=*) '[FAIL] ', TRIM(msg)
        STOP 
    ELSE
        WRITE(UNIT=OUTPUT_UNIT, FMT=*) '[PASS] ', TRIM(msg)
    END IF

    l(:) = header_cls(i0)
    i(:) = header_cls(l)
    WRITE(msg, "('Piping integer(', I0, ') throug logical(', I0, ')!')") &
        KIND(i0), KIND(l)
    IF ( STORAGE_SIZE(l)<STORAGE_SIZE(i0) ) THEN
        WRITE(UNIT=OUTPUT_UNIT, FMT=*) '[SKIP] ', TRIM(msg)
    ELSE IF ( .NOT. ALL( i == i0) ) THEN
        WRITE(UNIT=OUTPUT_UNIT, FMT=*) '[FAIL] ', TRIM(msg)
        STOP 
    ELSE
        WRITE(UNIT=OUTPUT_UNIT, FMT=*) '[PASS] ', TRIM(msg)
    END IF

    c(:) = header_cls(i0)
    i(:) = header_cls(c)
    WRITE(msg, "('Piping integer(', I0, ') throug character(', I0, ')!')") &
        KIND(i0), LEN(c)
    IF ( STORAGE_SIZE(c)<STORAGE_SIZE(i0) ) THEN
        WRITE(UNIT=OUTPUT_UNIT, FMT=*) '[SKIP] ', TRIM(msg)
    ELSE IF ( .NOT. ALL( i == i0) ) THEN
        WRITE(UNIT=OUTPUT_UNIT, FMT=*) '[FAIL] ', TRIM(msg)
        STOP 
    ELSE
        WRITE(UNIT=OUTPUT_UNIT, FMT=*) '[PASS] ', TRIM(msg)
    END IF



    
    ! Testing function iday()
    WRITE(UNIT=OUTPUT_UNIT, FMT='(/,A)') 'Tesing function iday()'

    ! 1st of January
    call check_iday(y=2000, m=01, d=01, j=001)
    ! 31st of December 2000
    call check_iday(y=2000, m=12, d=31, j=366) 
    ! 31st of December 2001
    call check_iday(y=2001, m=12, d=31, j=365) 
    ! 31st of December 2012
    call check_iday(y=2012, m=12, d=31, j=366) 
    ! 11th of October 2012
    call check_iday(y=2012, m=10, d=11, j=285) 
    

    WRITE(UNIT=OUTPUT_UNIT, FMT='()') 

CONTAINS
    
    ! Tests function iday()
    SUBROUTINE check_iday(y, m, d, j)
        INTEGER :: y, m, d, j
        LOGICAL :: l
        CHARACTER(LEN=128) :: msg

        l = j /= iday(year=y, month=m, day=d)

        WRITE(msg, "('iday(', I4, '/', I2, '/',  I2,') == ', I0 )") y, m, d, j

        IF ( l ) THEN
            WRITE(UNIT=ERROR_UNIT, FMT=*) '[FAIL] ', TRIM(msg)
            STOP 
        ELSE 
            WRITE(UNIT=OUTPUT_UNIT, FMT=*) '[PASS] ', TRIM(msg)
        END IF

    END SUBROUTINE

END PROGRAM sac_header_test


!> @file
!! $Date: 2013-10-30 16:25:25 +0100 (Wed, 30 Oct 2013) $
!! $Author: mauerberger $
!! $Revision: 743 $
!! @copyright GNU General Public License version 3 or later 
