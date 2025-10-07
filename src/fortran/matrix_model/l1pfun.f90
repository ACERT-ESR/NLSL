!  NLSL Version 1.5b 11/5/95
!----------------------------------------------------------------------
!                    =========================
!                        function L1PFUN
!                    =========================
!
!     Single-parameter function for minimization of experimental-calc'd
!     residuals using Brent's method minimization.
!
!----------------------------------------------------------------------
      function l1pfun( parm, bflag )
!
      use nlsdim
      use parcom
      use expdat
      use iterat
      use lmcom
!
      implicit none
      double precision l1pfun,parm(1)
      integer bflag,isp,tmpflg
!
      integer IPRINT
      double precision ZERO
      parameter (IPRINT=0,ZERO=0.0d0)
!
      double precision enorm
      external enorm
!
      iter=iter+1
      call lfun(ndatot,1,parm,fvec,fjac,MXPT,bflag)
      fnorm=enorm(ndatot,fvec)

!     -------------------------------------------
!      Update shifts when a new minimum is found
!     -------------------------------------------
      if (fnorm.lt.fnmin) then
         fnmin=fnorm
         do isp=1,nspc
            shft(isp)=shft(isp)+tmpshft(isp)
            tmpshft(isp)=ZERO
         end do
      end if

      tmpflg=IPRINT
      call lfun(ndatot,1,parm,fvec,fjac,MXPT,tmpflg)
!
      l1pfun=fnorm
!         
      return
      end
