! NLSL Version 1.9.0 beta 2/11/15
!*********************************************************************
!
!    STRUTL2       Fortran 90/95 does not have intrinsic functions for
!                  converting a string to uppercase, nor for finding
!                  the index of the first blank character in a string.
!                  Therefore, subroutine touppr and function itrim are
!                  being supplied here. They were previously included
!                  in the strutl.f collection.
!
!*********************************************************************
!
!----------------------------------------------------------------------
!                    =========================
!                        subroutine TOUPPR
!                    =========================
!
!   Converts string argument to uppercase
!----------------------------------------------------------------------
!
      subroutine touppr(string,lth)
      implicit none
      character(len=*) :: string
      character(len=1) :: chr
      integer :: i, ich, lth
!
      do i=1,lth
         chr=string(i:i)
         ich=iachar(chr)
         if (ich.ge.97 .and. ich.le.122) string(i:i)=achar(ich-32)
      end do
      return
      end subroutine touppr
!
!----------------------------------------------------------------------
!                    =========================
!                         function ITRIM
!                    =========================
!
!   Returns position of first blank in a string
!----------------------------------------------------------------------
!
      function itrim( string )
      implicit none
      integer :: itrim, j, lth
      character(len=*) :: string
      lth=len(string)
      do j=1,lth
         if (string(j:j).eq.' ') then
            itrim=j-1
            return
         end if
      end do
!     If here: no blanks in string, return full length
!     Differs from F77 index() intrinsic, which would return 0 
      itrim=lth
      return
      end function itrim
