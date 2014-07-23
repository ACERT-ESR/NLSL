c  VERSION 1.0  (NLSPMC version)   2/5/99
c----------------------------------------------------------------------
c                    =========================
c                        function P1PFUN
c                    =========================
c----------------------------------------------------------------------
      function p1pfun( parm,iflag )
c
      implicit none
c
      include 'limits.inc'
      include 'parms.inc'
      include 'datas.inc'
      include 'lmcomm.inc'
c
      integer iflag
      double precision p1pfun,parm,enorm
      external enorm
c
      iflag=1
      call pfun(nptot,1,parm,fvec,fjac,mxpt,iflag)
      if (ihltcmd.ne.0) return
      if (iflag.ge.0) p1pfun = enorm(nptot,fvec)/dsqrt(dfloat(nptot))
      return
      end
