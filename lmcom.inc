c Version 1.4  12/29/31
c----------------------------------------------------------------------
c                       =========================
c                          common block LMCOM
c                       =========================
c     Work arrays for Levenberg-Marquardt least squares fitting
c     procedure and parameters controlling termination criteria.
c
c     NOTE: File 'nlsdim.inc' should be included before this file
c           in the program source code.
c
c----------------------------------------------------------------------
      double precision corr,diag,fjac,fvec,gnvec,gradf,qtf,tcov,
     #                 work1,work2,work3,work4,x,ftol,xtol,gtol,factor
      integer ipvt,maxev,maxitr,nprint,mode,info,lmflag
c
      common /lmcom/ fjac(MXPT,MXJCOL),fvec(MXPT),x(MXJCOL),
     #               diag(MXJCOL),qtf(MXJCOL),corr(MXJCOL,MXJCOL),
     #               work1(MXJCOL),work2(MXJCOL),work3(MXJCOL),
     #               work4(MXPT),gnvec(MXJCOL),gradf(MXJCOL),
     #               tcov(MXTV,MXTV),ftol,gtol,xtol,factor,
     #               ipvt(MXJCOL),nprint,maxev,maxitr,mode,info,lmflag
c
      integer NFLMPR,NILMPR
      parameter (NFLMPR=4,NILMPR=4)
c
      double precision flmprm(NFLMPR)
      integer ilmprm(NILMPR)
      equivalence (flmprm(1),ftol),(ilmprm(1),maxev)
c
