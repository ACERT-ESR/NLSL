c Version 1.5.1b 11/6/95
c----------------------------------------------------------------------
c                    =========================
c                         common PARCOM
c                    =========================
c
c     Common block containing a list of quantities associated with
c     each of the fitting parameters being varied (i.e. each parameter
c     in the x array for LMDER-family programs).
c
c     NOTE: the file 'nlsdim.inc' must be included before this file.
c
c prmax  : Vector containing upper limits for each parameter in x array
c prmin  : Vector containing lower limits for each parameter in x array
c prscl  : Vector containing desired absolute accuracy for each parameter
c          in x array
c xfdstp : Vector containing size of forward-difference step for each
c          parameter in x array
c xerr   : Vector containing uncertainty estimates for each parameter in x
c serval : List of values for the parameter being varied in a series
c wlb    : List of line-broadening widths for each spectrum in a series
c ibnd   : Flag for boundaries imposed on each parameter:
c          0=none, 1=minimum, 2=maximum, 3=both
c ixpr   : Index of each variable parameter appearing in the x array 
c          into the fepr array in /eprprm/
c ixst   : Secondary index of each parameter appearing in the x array 
c          identifying which site or spectrum in a series the parameter
c          is associated with
c ixx    : Index of each parameter into the variable parameter array, x
c          (0 if parameter is not being varied)
c iser   : Index of the parameter being varied in a series of spectra
c nser   : Number of values given for parameter <iser> in the series
c          (should equal number of spectra in the series)
c nsite  : Number of sites defined for a given spectrum
c nprm   : Number of parameters being varied
c ptol   : Parameter convergence tolerance for 1-parameter searches
c pftol  : Function convergence tolerance for 1-parameter searches
c pbound : Search bound for 1-parameter searches
c srange : Allowed shifting range
c----------------------------------------------------------------------
c

       module parcom
       use nlsdim

c      double precision prmax,prmin,prscl,serval,xfdstp,xerr,fparm,
c     #                 ctol,ptol,pftol,pstep,pbound,srange
c      integer iparm,ibnd,ixpr,ixst,ixx,iser,nser,nsite,nprm,njcol,
c     #        nshift,noneg,ixp1p,ixs1p,mxpitr,itridg,iitrfl,jacobi,
c     #        output
c      logical mtxclc
c      character*9 tag
c
      double precision, target, save :: fparm(NFPRM,MXSITE),
     #                prmax(MXVAR),prmin(MXVAR),
     #                prscl(MXVAR),xfdstp(MXVAR),xerr(MXJCOL),
     #                serval(MXSPC),ctol,ptol,pftol,pstep,pbound,
     #                srange
      integer, target, save :: iparm(NIPRM,MXSITE),ixx(NFPRM,MXSITE),
     #                ibnd(MXVAR),ixpr(MXVAR),ixst(MXVAR),
     #                iser,nser,nsite,nprm,njcol,nshift,noneg,itridg,
     #                iitrfl,jacobi,ixp1p,ixs1p,mxpitr,output,mtxclc,
     #                tag(MXJCOL)
c

      end module parcom
