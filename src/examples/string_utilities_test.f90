!> @example string_utilities_test.f90
!!
!! Program which gives examples how to use procedures of the module coordinate_utilities_mod. 
!! The following examples are carried out:
!!  * Mapping strings to upper case with `raise_case()`
!!  * Mapping strings to lower case with `lower_case()`
!!  * Printing the integer N as a string with `i2c()`
!!  * Replacing single chars in the strings with `replace_char()`
!!  * Sweeping single chars in the strings with `sweep_char()`
!!  * Empty entries dropped and whitespaces are removed with `strip()`
!!  * Concatenates an array strings into a single string with `array_to_string()`
!!
!! Compile with simply using the Makefile:
!!
!!     make string_utilities_test
!!
!! Run the program by executing:
!!
!!     ./string_utilities_test
PROGRAM string_utilities_test
    USE string_utilities_mod, ONLY : i2c, raise_case, lower_case, &
                                     sweep_char, replace_char, &
                                     array_to_string, strip

    IMPLICIT NONE

    CHARACTER (LEN=*), PARAMETER :: string = 'Ses3d-NT'

    CHARACTER, PARAMETER :: old = 'd', new = 'D'

    INTEGER, PARAMETER :: n = 12345


    ! example  mapping strings to upper case
    IF (raise_case(string) == 'SES3D-NT') THEN
        WRITE(*,*) "PASS: raise_case('Ses3d-NT') returns string 'SES3D-NT'"
    ELSE
        STOP "FAIL: raise_case('Ses3d-NT') /= 'SES3D-NT'"
    END IF


    ! example  mapping strings to lower case
    IF (lower_case(string) == 'ses3d-nt') THEN
        WRITE(*,*) "PASS: lower_case('Ses3d-NT') returns string 'ses3d-nt'"
    ELSE
        STOP "FAIL: raise_case('Ses3d-NT') /= 'ses3d-nt'"
    END IF


    ! example printed form of the integer N as a string
    IF (i2c(n) == '12345') THEN
        WRITE(*,*) "PASS: i2c('12345') returns integer n=12345 as a string '12345'"
    ELSE
        STOP "FAIL: i2c('12345') /= '12345'"
    END IF


    ! example replacing single chars in the strings
    IF (replace_char(string,old,new) == 'Ses3D-NT') THEN
        WRITE(*,*) "PASS: replace_char('Ses3d-NT','d','D') returns string 'Ses3D-NT'"
    ELSE
        STOP "replace_char('Ses3d-nt','d','D') /= 'Ses3D-NT'"
    END IF


    ! example sweeping single chars in the strings
    IF (sweep_char(string,old) == 'Ses3-NT') THEN
        WRITE(*,*) "PASS: sweep_char('Ses3d-NT','d') returns string 'Ses3-NT'"
    ELSE
        STOP "sweep_char('Ses3d-nt','d') /= 'Ses3-NT'"
    END IF


    ! example how to use strip
    IF ( SIZE(strip([" "," "])) /= 0 ) THEN
        STOP 'FAIL: SIZE(strip([" "," "])) /= 0'
    ELSE
        WRITE(*,*) 'PASS: strip() of [" "," "] returns an empty array. '
    END IF

    IF ( ALL(strip( ["a ", "  ", " b", "cd"] ) /= ["a ", "b ", "cd"]) ) THEN
       STOP 'FAIL: strip( ["a ", "  ", " b", "cd"] ) /= ["a ", "b ", "cd"] '
    ELSE
        WRITE(*,*) 'PASS: strip() of ["a ", "  ", " b", "cd"] returns a &
            &size 3 array ["a ", "b ", "cd"].' 
    END IF


    ! exampel how to use array_to_string()
    IF ( array_to_string(['a ',' b','cd']) .NE. "'a', ' b', 'cd'" ) THEN
        STOP "FAIL: array_to_string(['a ',' b','cd'] .NE. ""'a', ' b', 'cd'"" "
    ELSE
        WRITE(*,*) "PASS: array_to_string() of ['a ',' b','cd'] returns ""'a', ' b', 'cd'"". "
    END IF

END PROGRAM string_utilities_test


!> @file
!! $Date: 2013-11-19 14:49:58 +0100 (Tue, 19 Nov 2013) $
!! $Author: mmelnyk $
!! $Revision: 798 $
!! @copyright GNU General Public License version 3 or later 
