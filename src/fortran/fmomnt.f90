! NLSL Version 1.4 10/10/94
!----------------------------------------------------------------------
!                    =========================
!                       function FMOMNT  
!                    =========================
!
!  Returns the first moment of a tabulated function in index units  
!  (first x-value is 1, second is 2, etc.). If f is tabulated as zeroth
!  derivative data, the first moment is defined as
!
!    intgrl[ x*f(x) dx ] / intgrl[ f(x) dx ]
!
!  If f is first derivative, the first moment returned is
!
!    intgrl[ x* intgrl[ f(x) dx ] dx ] / intgrl[ intgrl[ f(x) dx ] dx ]
!
!  Inputs:
!    ay   Tabulated function values (equal x-spacing assumed)
!    n       Number of points in ay
!    ideriv  Derivative flag (0=zeroth, 1=first)
!----------------------------------------------------------------------
!
      function fmomnt( ay, n, ideriv )
      implicit none
      integer ideriv,n
      double precision ay(n),fmomnt 
!
      integer i
      double precision ai,ain,am1,sglsum
!
      double precision ZERO
      parameter(ZERO=0.0d0)
!
! ----
!
      am1=ZERO
      ain=ZERO
!                                 ***     0th derivative data
      if (ideriv.eq.0) then
         do 2 i=1,n
            ai=dfloat(i)
            am1=am1+ai*ay(i)
            ain=ain+ay(i)
 2       continue
      else
!                                 ***     1st derivative data
         sglsum=ZERO
         do 4 i=1,n
            ai=dfloat(i)
            sglsum=sglsum+ay(i)
            am1=am1+ai*sglsum
            ain=ain+sglsum
 4       continue
      end if
!
      if (ain.ne.ZERO) then
         fmomnt = am1/ain
      else 
         fmomnt=ZERO
      end if
      return
      end
