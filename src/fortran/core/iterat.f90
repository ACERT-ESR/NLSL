! NLSL Version 1.9.0 beta 2/12/15
!----------------------------------------------------------------------
!                    =========================
!                          module ITERAT
!                    =========================
!
!   Iteration counter, residual norm, and related statistical quantities
!   for nonlinear least-squares procedure
!----------------------------------------------------------------------
!
      module iterat
      implicit none
!
      integer, save :: iter, iwflag
      double precision , save :: fnorm, chisqr, rdchsq, 
     #                ch2bnd, chnbnd, confid, delchi1,
     #                qfit, f2bnd, fnbnd, tbound, fnmin
      logical, save :: newitr, covarOK, xreset
!
      end module iterat






