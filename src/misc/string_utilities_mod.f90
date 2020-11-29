!> @file
!! Ses3d-NT - simulation of elastic wave propagation in spherical sections
!! 
!! This material was produced under U.S. Government contract DE-AC52-06NA25396
!! for Los Alamos National Laboratory (LANL), which is operated by Los Alamos
!! National Security, LLC for the U.S. Department of Energy. The U.S. Government
!! has rights to use, reproduce, and distribute this software.  NEITHER THE
!! GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY, LLC MAKES ANY WARRANTY, EXPRESS
!! OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.  If 
!! software is modified to produce derivative works, such modified software 
!! should be clearly marked, so as not to confuse it with the version available 
!! from LANL.
!! 
!! Additionally, this program is free software; you can redistribute it and/or
!! modify it under the terms of the GNU General Public License as published by
!! the Free Software Foundation; either version 2 of the License, or (at your
!! option) any later version. Accordingly, this program is distributed in the
!! hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
!! implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!! See the GNU General Public License for more details.
!!
!! $Date: 2013-11-01 18:22:04 +0100 (Fri, 01 Nov 2013) $
!! $Author: mmelnyk $
!! $Revision: 755 $
!> @copyright Copyright 2012.  Los Alamos National Security, LLC.



!! This module provides a collection of useful character string-related
!! utilities
module string_utilities_mod

  implicit none

  private
  
  public :: raise_case, lower_case, i2c, &
            replace_char, sweep_char, array_to_string, strip
  
contains

   
!> Returns a string of the same length and value as CS, but with each
!! lowercase letter replaced by its uppercase counterpart
!> @param cs string
  elemental function raise_case (cs)
    character(len=*), intent(in) :: cs
    character(len=len(cs)) :: raise_case
    integer :: i
    do i = 1, len(cs)
      if (iachar(cs(i:i)) >= iachar('a') .and. iachar(cs(i:i)) <= iachar('z')) then
        raise_case(i:i) = achar(iachar(cs(i:i)) - iachar('a') + iachar('A'))
      else
        raise_case(i:i) = cs(i:i)
      end if
    end do
  end function raise_case

!> Returns a string of the same length and value as CS, but with each
!! uppercase letter replaced by its lowercase counterpart  
!> @param cs string
  elemental function lower_case (cs)
    character(len=*), intent(in) :: cs
    character(len=len(cs)) :: lower_case
    integer :: i
    do i = 1, len(cs)
      if (iachar(cs(i:i)) >= iachar('A') .and. iachar(cs(i:i)) <= iachar('Z')) then
        lower_case(i:i) = achar(iachar(cs(i:i)) - iachar('A') + iachar('a'))
      else
        lower_case(i:i) = cs(i:i)
      end if
    end do
  end function lower_case

   
!> @brief Returns the printed form of the integer N as a string.
!!
!> @details The length of the string is precisely the number of characters
!! required to express the integer without leading or trailing blanks
  pure function i2c (n) result (s)
  
    integer, intent(in) :: n
    character(len=output_length(n)) :: s
    
    integer :: absn
    integer :: pos, digit
    
    absn = abs(n)
    if (absn < 0) then ! Can't treat -huge(0)-1 properly
      s = repeat('*', len(s))
      return
    end if
    
    if (n < 0) then
      s(1:1) = '-'
    else if (n == 0) then
      s = '0'
      return
    end if
    
    pos = len(s)
    do while (absn > 0)
      digit = mod(absn, 10) + 1
      s(pos:pos) = '0123456789'(digit:digit)
      pos = pos - 1
      absn = absn / 10
    end do
    
  end function i2c


!> Returns the required length  
  pure function output_length (n) result (l)
  
    integer, intent(in) :: n
    integer :: l
    
    select case (n)
    case (:-1)  ! negative
      l = int(log10(-dble(n)+0.1d0) + 2)
    case (0)
      l = 1
    case (1:)   ! positive
      l = int(log10(n+0.1d0) + 1)
    end select
    
  end function output_length



!> Returns a string of the same length as STR, but with each
!! character OLD replaced by NEW.
!> @param str string
!> @param old single character
!> @param new single character
  elemental function replace_char( str, old, new ) result( res )
    character(*), intent(in) :: str
    character, intent(in) :: old, new
    character :: tmp(len(str))
    character(len=len(str)) :: res
    tmp = transfer( str, tmp )
    tmp = merge( tmp, new, tmp /= old )
    res = transfer(tmp, str) 
  end function replace_char


!> Returns a string STR with all characters CHR removed
!> @param str string
!> @param chr single character
  pure function sweep_char( str, chr ) result( res )
    character(len=*), intent(in) :: str
    character, intent(in) :: chr
    character :: tmp(len(str))
    character(len=str_len(str,chr)) :: res
    ! Transfer str into an array of characters 
    tmp = transfer( str, tmp )
    ! Sweep all characters "chr" and transfer back into a string
    res = transfer( pack( tmp, tmp /= chr ), repeat('a',str_len(str,chr) ) ) 
  end function sweep_char


  !> Returns length of the string
  pure function str_len( str, chr ) 
    character(len=*), intent(in) :: str
    character, intent(in) :: chr
    character :: tmp(len(str))
    integer :: str_len
    ! Transfer str into an array of characters 
    tmp = transfer( str, tmp )
    ! Determine length of resulting string
    str_len = size( pack( tmp, tmp /= chr ) )
  end function str_len


  !> Returns a copy of an array of strings 'str' with empty entries 
  !! dropped. In addition all leading whitespaces are removed. 
  PURE FUNCTION strip(str)
    CHARACTER(LEN=*), INTENT(IN) :: str(:)
    CHARACTER(LEN=LEN(str(1))), ALLOCATABLE :: strip(:)
    CHARACTER(LEN=LEN(str(1))) :: str_(SIZE(str))
    
    ! Get a copy of str with leading whitespace removed 
    str_(:) = ADJUSTL(str)

    ! Allocate, assign and return non-empty entities
    strip = PACK( str_, str_ /= '' ) 
    
  END FUNCTION strip

  !> Concatenates an array strings into a single string with 
  !! trailing whitespaces removed. 
  !! @todo make sep and qte optional
  !! @bug Crash if sep="'" 
  PURE FUNCTION array_to_string(array) RESULT(str)
    CHARACTER(LEN=*), INTENT(IN) :: array(:)
    CHARACTER(LEN=*), PARAMETER :: sep=", ", qte="'"
    CHARACTER(LEN=SUM(LEN_TRIM(array)) + 2*LEN(qte)*SIZE(array) + LEN(sep)*(SIZE(array) - 1)) :: str
    INTEGER :: i
    CHARACTER(LEN=30) :: fmt_

    IF ( SIZE(array) == 0 ) &
        RETURN

    WRITE(fmt_, FMT="(A,I0,A,A,A)") "(", SIZE(array), "(A,:,'", sep, "'))"
    WRITE(str, FMT=fmt_) ( qte//TRIM(array(i))//qte, i=1, SIZE(array) )

  END FUNCTION array_to_string

end module string_utilities_mod
