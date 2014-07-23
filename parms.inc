c----------------------------------------------------------------------
c                    =======================
c                         common PARMS
c                    =======================
c
c     Common block containing a list of quantities associated with
c     each of the fitting parameters being varied (i.e. each parameter
c     in the x array for LMDER-family programs).
c
c     NOTE: the file 'limits.inc' must be included before this file.
c
c prmax  : Vector containing upper limits for each parameter in x array
c prmin  : Vector containing lower limits for each parameter in x array
c prscl  : Vector containing desired absolute accuracy for each parameter
c          in x array
c xfdstp : Vector containing size of forward-difference step for each
c          parameter in x array
c xerr   : Vector containing uncertainty estimates for each parameter in x
c ibnd   : Flag for boundaries imposed on each parameter:
c          0=none, 1=minimum, 2=maximum, 3=both
c ixpr   : Index of each variable parameter into the fepr array in /eprprm/
c nser   : Number of spectra in the series
c nsparm : Number of series parameters set
c          (should be 8 or 9 before data command)
c nprm   : Number of parameters being varied
c matrxok: Define what is needed to give spectra.  0=need matrix, 
c		1=have eval,evec, 2=have xspec, 3=have spec.
c----------------------------------------------------------------------
c
      double precision fparm,prmax,prmin,prscl,xfdstp,xerr,
     #                 stol,sstep,smn,smx,stval,fnorm
      integer iparm,idebug,ievec,iter,ihltcmd,ihltfit,
     #	ixx,ibnd,ixpr,nser,nsparm,nprm,negegv

      logical gflag,aflag,rotflg,potflg,uniflg,dataOK
      integer matrxok
      character*6 tag
      common /parms/ fparm(NFPRM,MXSPEC,MXSITE),prmax(MXVAR),
     #	prmin(MXVAR),
     #                prscl(MXVAR),xfdstp(MXVAR),xerr(MXVAR),
     #                stol,sstep,smn,smx,stval,fnorm,
     # 		      iparm(NIPRM,MXSPEC,MXSITE),idebug,ievec,
     #                gflag,aflag,rotflg,potflg,uniflg,
     #		      iter,ihltcmd,ihltfit,
     #                ixx(NFPRM,MXSPEC,MXSITE),ibnd(MXVAR),
     #                ixpr(MXVAR),nser,nsparm,nprm,negegv,
     #                matrxok(MXSITE),
     #                dataOK,tag(MXVAR)

