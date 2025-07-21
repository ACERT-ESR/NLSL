c NLSL Version 1.9.0 beta 2/9/15
c----------------------------------------------------------------------
c                         ====================
c                            module EXPDAT
c                         ====================
c     Experimental data used for least-squares fitting and associated
c     parameters. 
c
c     NOTE: module nlsdim is used by this module.
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
      module expdat
      use nlsdim
      implicit none
c
      double precision, save :: data(MXPT),spltmp(MXINP,3),rmsn(MXSPC),
     #                sbi(MXSPC),sdb(MXSPC),srng(MXSPC),shft(MXSPC),
     #                sb0(MXSPC),sphs(MXSPC),spsi(MXSPC),slb(MXSPC),
     #                tmpshft(MXSPC)
c
      integer, save :: iform(MXSPC),ibase(MXSPC),nft(MXSPC),
     #                npts(MXSPC),ishft(MXSPC),idrv(MXSPC),ixsp(MXSPC),
     #                nrmlz(MXSPC),nspc,nwin,ndatot,ishglb,
     #                inform,bcmode,drmode,
     #                nspline,shftflg,normflg,written
c
      character*30, save :: dataid(MXSPC)
      character*20, save :: wndoid(MXSPC)
c
      end module expdat
