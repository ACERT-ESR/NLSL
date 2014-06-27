c NLSL Version 1.5.1 beta 1/20/96
c----------------------------------------------------------------------
c                    =========================
c                      common block MSPCTR
c                    =========================
c Defines common block for working arrays used for multiple-site 
c spectral fitting
c
c Note: file "nlsdim.inc" must be included in the source code before
c this file
c
c    spectr    Storage for calculated spectra for each site, spectrum
c
c    wspec     Work array for temporary storage of calculated spectra
c              Used by sshift to store correlation function for each
c              spectrum  
c              Used before call to sscale (array is replaced by QR
c              decomposition)
c
c    sfac      Scale factors for each site and spectrum
c
c    iscal     Flag indicating whether each site is to be automatically
c              scaled in the function calculation 1=auto scale
c
c    iscglb    Flag indicating whether all sites are to be automatically
c              scaled (1=auto scale)
c
c----------------------------------------------------------------------
c
      integer iscal,iscglb
      double precision spectr,wspec,sfac
      common /mspctr/ spectr(MXPT,MXSITE),wspec(MXPT,MXSITE),
     #                sfac(MXSITE,MXSPC),iscal(MXSITE),iscglb

