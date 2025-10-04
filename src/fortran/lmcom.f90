! NLSL Version 1.9.0 beta 2/13/15
!----------------------------------------------------------------------
!                       =========================
!                             module LMCOM
!                       =========================
!     Work arrays for Levenberg-Marquardt least squares fitting
!     procedure and parameters controlling termination criteria.
!
!     NOTE: This module uses the nlsdim.inc module.
!
!----------------------------------------------------------------------
!
      module lmcom
      use nlsdim
      implicit none
!
      double precision, save :: fjac(MXPT,MXJCOL)
      double precision, save :: fvec(MXPT)
      double precision, save :: x(MXJCOL)
      double precision, save :: diag(MXJCOL)
      double precision, save :: qtf(MXJCOL)
      double precision, save :: corr(MXJCOL,MXJCOL)
      double precision, save :: work1(MXJCOL)
      double precision, save :: work2(MXJCOL)
      double precision, save :: work3(MXJCOL)
      double precision, save :: work4(MXPT)
      double precision, save :: gnvec(MXJCOL)
      double precision, save :: gradf(MXJCOL)
      double precision, save :: tcov(MXTV,MXTV)
!
      integer, save :: ipvt(MXJCOL)
      integer, save :: nprint
      integer, save :: lmflag
!
      double precision, target, save  :: flmprm(NFLMPR)
      double precision, pointer, save :: ftol, gtol, xtol, factor
!
      integer, target, save  :: ilmprm(NILMPR)
      integer, pointer, save :: maxev, maxitr, mode, info
      character*10, dimension(NFLMPR), save :: flmprm_name
      character*10, dimension(NILMPR), save :: ilmprm_name
!
      contains
!
      subroutine lmcom_init
      implicit none
!
      ftol   => flmprm(1)
      flmprm_name(1) = "ftol"
      gtol   => flmprm(2)
      flmprm_name(2) = "gtol"
      xtol   => flmprm(3)
      flmprm_name(3) = "xtol"
      factor => flmprm(4)
      flmprm_name(4) = "bound"
!
      maxev  => ilmprm(1)
      ilmprm_name(1) = "maxfun"
      maxitr => ilmprm(2)
      ilmprm_name(2) = "maxitr"
      mode   => ilmprm(3)
      ilmprm_name(3) = "mode"
      info   => ilmprm(4)
      ilmprm_name(4) = "info"
!
      end subroutine lmcom_init
!
      end module lmcom
