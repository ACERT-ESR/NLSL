c     Version 1.5  5/2/94
c----------------------------------------------------------------------
c                    =========================
c                       common  BASWT
c                    =========================
c
c   Storage for tridiagonal matrix and vectors used in Lanczos and
c   conjugate gradients routines
c
c   NOTE: File 'stddim.inc' must be included before this file
c----------------------------------------------------------------------
      integer ixfld,keep
      double precision bsswt,totwt
      common /baswt/ bsswt(MXDIM),totwt(MXDIM),ixfld(MXDIM),keep(MXDIM)
