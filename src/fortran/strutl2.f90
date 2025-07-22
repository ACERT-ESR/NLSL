c NLSL Version 1.9.0 beta 2/11/15
c*********************************************************************
c
c    STRUTL2       Fortran 90/95 does not have intrinsic functions for
c                  converting a string to uppercase, nor for finding
c                  the index of the first blank character in a string.
c                  Therefore, subroutine touppr and function itrim are
c                  being supplied here. They were previously included
c                  in the strutl.f collection.
c
c*********************************************************************
c
c----------------------------------------------------------------------
c                    =========================
c                        subroutine TOUPPR
c                    =========================
c
c   Converts string argument to uppercase
c----------------------------------------------------------------------
c
      subroutine touppr(string,lth)
      implicit none
      character(len=*) :: string
      character(len=1) :: chr
      integer :: i, ich, lth
c
      do i=1,lth
         chr=string(i:i)
         ich=iachar(chr)
         if (ich.ge.97 .and. ich.le.122) string(i:i)=achar(ich-32)
      end do
      return
      end subroutine touppr
c
c----------------------------------------------------------------------
c                    =========================
c                         function ITRIM
c                    =========================
c
c   Returns position of first blank in a string
c----------------------------------------------------------------------
c
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
c     If here: no blanks in string, return full length
c     Differs from F77 index() intrinsic, which would return 0 
      itrim=lth
      return
      end function itrim
