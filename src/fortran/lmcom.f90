c NLSL Version 1.9.0 beta 2/13/15
c----------------------------------------------------------------------
c                       =========================
c                             module LMCOM
c                       =========================
c     Work arrays for Levenberg-Marquardt least squares fitting
c     procedure and parameters controlling termination criteria.
c
c     NOTE: This module uses the nlsdim.inc module.
c
c----------------------------------------------------------------------
c
      module lmcom
      use nlsdim
      implicit none
c
      double precision, save ::
     #                 fjac(MXPT,MXJCOL), fvec(MXPT), x(MXJCOL),
     #                 diag(MXJCOL), qtf(MXJCOL), corr(MXJCOL,MXJCOL),
     #                 work1(MXJCOL), work2(MXJCOL), work3(MXJCOL),
     #                 work4(MXPT), gnvec(MXJCOL), gradf(MXJCOL),
     #                 tcov(MXTV,MXTV)
c
      integer, save :: ipvt(MXJCOL), nprint, lmflag
c
      integer, parameter :: NFLMPR=4, NILMPR=4
c
      double precision, target, save  :: flmprm(NFLMPR)
      double precision, pointer, save :: ftol, gtol, xtol, factor
c
      integer, target, save  :: ilmprm(NILMPR)
      integer, pointer, save :: maxev, maxitr, mode, info
      character*10, dimension(NFLMPR), save :: flmprm_name
c
      contains
c
      subroutine lmcom_init
      implicit none
c
      ftol   => flmprm(1)
      flmprm_name(1) = "ftol"
      gtol   => flmprm(2)
      flmprm_name(2) = "gtol"
      xtol   => flmprm(3)
      flmprm_name(3) = "xtol"
      factor => flmprm(4)
      flmprm_name(4) = "factor"
c
      maxev  => ilmprm(1)
      maxitr => ilmprm(2)
      mode   => ilmprm(3)
      info   => ilmprm(4)
c
      end subroutine lmcom_init
c
      end module lmcom
