c  NLSL Version 1.5b 11/5/95
c----------------------------------------------------------------------
c                    =========================
c                        function L1PFUN
c                    =========================
c
c     Single-parameter function for minimization of experimental-calc'd
c     residuals using Brent's method minimization.
c
c----------------------------------------------------------------------
      function l1pfun( parm, bflag )
c
      use nlsdim
      use parcom
      use expdat
      use iterat
      use lmcom
c
      implicit none
      double precision l1pfun,parm(1)
      integer bflag,isp,tmpflg
c
      integer IPRINT
      double precision ZERO
      parameter (IPRINT=0,ZERO=0.0d0)
c
      double precision enorm
      external enorm
c
      iter=iter+1
      call lfun(ndatot,1,parm,fvec,fjac,MXPT,bflag)
      fnorm=enorm(ndatot,fvec)

c     -------------------------------------------
c      Update shifts when a new minimum is found
c     -------------------------------------------
      if (fnorm.lt.fnmin) then
         fnmin=fnorm
         do isp=1,nspc
            shft(isp)=shft(isp)+tmpshft(isp)
            tmpshft(isp)=ZERO
         end do
      end if

      tmpflg=IPRINT
      call lfun(ndatot,1,parm,fvec,fjac,MXPT,tmpflg)
c
      l1pfun=fnorm
c         
      return
      end
