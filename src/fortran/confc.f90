! Version 1.4 1/1/95
!----------------------------------------------------------------------
!                    =========================
!                      subroutine CONFC
!                    =========================
!----------------------------------------------------------------------
      subroutine confc( line )
!
      use iterat
      use stdio
!
      implicit none
      character*80 line
!
      character*30 token
      double precision fval
      integer lth
!
      logical ftoken
      external ftoken
!
      call gettkn(line,token,lth)
      if (lth.eq.0) then
         write (luttyo,1000)
      end if
      if ( ftoken(token,lth,fval) ) then
         if (fval.lt.0.0d0 .or. fval.gt.1.0D0) then
            write(luttyo,1002)
         else
            confid=fval
         end if
      else
         write (luttyo,1003) token(:lth)
      end if
      return
!
! ##### Formats ##################################################
!
 1000 format(' *** Confidence level missing ***')
 1002 format(' *** Confidence level must be between 0 and 1 ***')
 1003 format(' *** Illegal floating point value: ''',a,''' ***')
      end
