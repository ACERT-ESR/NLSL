! NLSL Version 1.9.0 beta 2/4/15
!----------------------------------------------------------------------
!                    =========================
!                          module PARCOM
!                    =========================
!
!>@brief a list of quantities (1) associated with each of the fitting
!> parameters being varied (i.e. each parameter in the x array for
!> LMDER-family programs) as well as (2) a set of parameters that don't
!> fall into this category
!
!     NOTE: this module uses nlsdim, so compile nlsdim.f90 first.
!
! prmax  : Vector containing upper limits for each parameter in x array
! prmin  : Vector containing lower limits for each parameter in x array
! prscl  : Vector containing desired absolute accuracy for each parameter
!          in x array
! xfdstp : Vector containing size of forward-difference step for each
!          parameter in x array
! xerr   : Vector containing uncertainty estimates for each parameter in x
! serval : List of values for the parameter being varied in a series
! wlb    : List of line-broadening widths for each spectrum in a series
! ibnd   : Flag for boundaries imposed on each parameter:
!          0=none, 1=minimum, 2=maximum, 3=both
! ixpr   : Index of each variable parameter appearing in the x array 
!          into the fepr array in module eprprm
! ixst   : Secondary index of each parameter appearing in the x array 
!          identifying which site or spectrum in a series the parameter
!          is associated with
! ixx    : Index of each parameter into the variable parameter array, x
!          (0 if parameter is not being varied)
! iser   : Index of the parameter being varied in a series of spectra
! nser   : Number of values given for parameter <iser> in the series
!          (should equal number of spectra in the series)
! nsite  : Number of sites defined for a given spectrum
! nprm   : Number of parameters being varied
! ptol   : Parameter convergence tolerance for 1-parameter searches
! pftol  : Function convergence tolerance for 1-parameter searches
! pbound : Search bound for 1-parameter searches
! srange : Allowed shifting range
!----------------------------------------------------------------------
!
       module parcom
       use nlsdim
       implicit none
!
!      double precision prmax,prmin,prscl,serval,xfdstp,xerr,fparm,
!     #                 ctol,ptol,pftol,pstep,pbound,srange
!      integer iparm,ibnd,ixpr,ixst,ixx,iser,nser,nsite,nprm,njcol,
!     #        nshift,noneg,ixp1p,ixs1p,mxpitr,itridg,iitrfl,jacobi,
!     #        output
!      logical mtxclc
!      character*9 tag
!
      double precision, target, save :: fparm(NFPRM,MXSITE)
!
      double precision, save :: prmax(MXVAR)
      double precision, save :: prmin(MXVAR)
      double precision, save :: prscl(MXVAR)
      double precision, save :: xfdstp(MXVAR)
      double precision, save :: xerr(MXJCOL),serval(MXSPC),ctol,ptol,
     #                pftol,pstep,pbound,srange
      integer, target, save :: iparm(NIPRM,MXSITE)
!
      integer, save :: ixx(NFPRM,MXSITE)
      integer, save :: ibnd(MXVAR)
      integer, save :: ixpr(MXVAR)
      integer, save :: ixst(MXVAR)
      integer, save :: iser,nser,njcol,nshift,noneg,itridg,
     #                iitrfl,jacobi,ixp1p,ixs1p,mxpitr,output
      integer, save :: nprm
      integer, target, save :: nsite
      integer, pointer, save :: pnsite => nsite
      logical, save :: mtxclc
      character*9, save :: tag(MXJCOL)
!
      end module parcom
