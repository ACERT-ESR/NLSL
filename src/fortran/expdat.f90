! NLSL Version 1.9.0 beta 2/9/15
!----------------------------------------------------------------------
!                         ====================
!                            module EXPDAT
!                         ====================
!     Experimental data used for least-squares fitting and associated
!     parameters. 
!
!     NOTE: module nlsdim is used by this module.
!
!     data  : Experimental data points
!     rmsn  : Noise level for each spectrum
!     sbi   : Starting fields for each spectrum
!     sdb   : Field step size for each spectrum
!     srng  : Field range for each spectrum
!     sb0   : Resonance field (for average g-value) for each spectrum
!     sphs  : Microwave phase for each spectrum
!     spsi  : Director tilt for each spectrum
!     scale : Scale factor for each spectrum (calculated during fit)
!     iform : Input format of each datafile
!     ibase : Baseline correction method for each datafile
!     nrmlz : "Normalized" flag for each spectrum (0=no)
!     npts  : Number of data points in each spectrum
!     nft   : Smallest power of 2 larger than npts
!     idrv  : Flag=1 if data are in 1st derivative form
!     ixsp  : Starting index of each spectrum in data array
!     nspc  : Number of spectra
!     nwin  : Number of open windows
!     ndatot: Total number of points (all spectra)
!     ishglb: Flag=1 if shifting enabled for all spectra
!     dataid: Identification string for each datafile
!     wndoid: Identification string for each window
!     written: Flags whether current calculation has been written to disk
!
!----------------------------------------------------------------------
      module expdat
      use nlsdim
      implicit none
!
      double precision, save :: data(MXPT)
      double precision, save :: spltmp(MXINP,3)
      double precision, save :: rmsn(MXSPC)
      double precision, save :: sbi(MXSPC)
      double precision, save :: sdb(MXSPC)
      double precision, save :: srng(MXSPC)
      double precision, save :: shft(MXSPC)
      double precision, save :: sb0(MXSPC)
      double precision, save :: sphs(MXSPC)
      double precision, save :: spsi(MXSPC)
      double precision, save :: slb(MXSPC)
      double precision, save :: tmpshft(MXSPC)
!
      integer, save :: iform(MXSPC)
      integer, save :: ibase(MXSPC)
      integer, save :: nft(MXSPC)
      integer, save :: npts(MXSPC)
      integer, save :: ishft(MXSPC)
      integer, save :: idrv(MXSPC)
      integer, save :: ixsp(MXSPC)
      integer, save :: nrmlz(MXSPC)
      integer, save :: nspc
      integer, save :: nwin
      integer, save :: ndatot
      integer, save :: ishglb
      integer, save :: inform
      integer, save :: bcmode
      integer, save :: drmode
      integer, save :: nspline
      integer, save :: shftflg
      integer, save :: normflg
      integer, save :: written
!
      character*30, save :: dataid(MXSPC)
      character*20, save :: wndoid(MXSPC)
!
      end module expdat
