c NLSL Version 1.5.1 beta 11/25/95
c----------------------------------------------------------------------
c                         ====================
c                          common block EXPDAT
c                         ====================
c     Experimental data used for least-squares fitting and associated
c     parameters. 
c
c     NOTE: File 'nlsdim.inc' should be included before this file
c           in the program source code.
c
c     data  : Experimental data points
c     rmsn  : Noise level for each spectrum
c     sbi   : Starting fields for each spectrum
c     sdb   : Field step size for each spectrum
c     srng  : Field range for each spectrum
c     sb0   : Resonance field (for average g-value) for each spectrum
c     sphs  : Microwave phase for each spectrum
c     spsi  : Director tilt for each spectrum
c     scale : Scale factor for each spectrum (calculated during fit)
c     iform : Input format of each datafile
c     ibase : Baseline correction method for each datafile
c     nrmlz : "Normalized" flag for each spectrum (0=no)
c     npts  : Number of data points in each spectrum
c     nft   : Smallest power of 2 larger than npts
c     idrv  : Flag=1 if data are in 1st derivative form
c     ixsp  : Starting index of each spectrum in data array
c     nspc  : Number of spectra
c     nwin  : Number of open windows
c     ndatot: Total number of points (all spectra)
c     ishglb: Flag=1 if shifting enabled for all spectra
c     dataid: Identification string for each datafile
c     wndoid: Identification string for each window
c     written: Flags whether current calculation has been written to disk
c
c----------------------------------------------------------------------
      double precision data,rmsn,sbi,sdb,sb0,sphs,spsi,shft,slb,spltmp,
     #                 srng,tmpshft
      integer iform,ibase,npts,nft,nrmlz,ishft,idrv,ixsp,nspc,nwin,
     #        ndatot,ishglb,inform,bcmode,drmode,nspline,shftflg,
     #        normflg,written
c
      character*30 dataid
      character*20 wndoid
c
      common /expdat/ data(MXPT),spltmp(MXINP,3),rmsn(MXSPC),
     #                sbi(MXSPC),sdb(MXSPC),srng(MXSPC),shft(MXSPC),
     #                sb0(MXSPC),sphs(MXSPC),spsi(MXSPC),slb(MXSPC),
     #                tmpshft(MXSPC),
     #                iform(MXSPC),ibase(MXSPC),nft(MXSPC),
     #                npts(MXSPC),ishft(MXSPC),idrv(MXSPC),ixsp(MXSPC),
     #                nrmlz(MXSPC),nspc,nwin,ndatot,ishglb,
     #                dataid(MXSPC),wndoid(MXSPC),inform,bcmode,drmode,
     #                nspline,shftflg,normflg,written
c
