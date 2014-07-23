c       VERSION 1.0     2/5/99
c*********************************************************************
c
c                         SUBROUTINE IPAR
c                         ===============
c
c        This integer function subroutine returns the value -1 
c        if its argument is odd, and +1 if it is even.
c
c*********************************************************************
c
      integer function ipar(num)
c
      integer num
c
c#####################################################################
c
      if (mod(num,2).eq.0) then
        ipar=1
      else
        ipar=-1
      end if
      return
      end
