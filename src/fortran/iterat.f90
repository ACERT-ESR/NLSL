c NLSL Version 1.9.0 beta 2/12/15
c----------------------------------------------------------------------
c                    =========================
c                          module ITERAT
c                    =========================
c
c   Iteration counter, residual norm, and related statistical quantities
c   for nonlinear least-squares procedure
c----------------------------------------------------------------------
c
      module iterat
      implicit none
c
      integer, save :: iter, weighted_flag
      double precision , save :: fnorm, chisqr, rdchsq, 
     #                ch2bnd, chnbnd, confid, delchi1,
     #                qfit, f2bnd, fnbnd, tbound, fnmin
      logical, save :: newitr, covarOK, xreset
c
      end module iterat






