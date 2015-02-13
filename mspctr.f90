c NLSL Version 1.9.0 beta 2/12/15
c----------------------------------------------------------------------
c                    =========================
c                          module MSPCTR
c                    =========================
c Defines working arrays used for multiple-site spectral fitting
c
c Note: module nlsdim is used by this module
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
c              scaled in the function calculation (1=auto scale)
c
c    iscglb    Flag indicating whether all sites are to be automatically
c              scaled (1=auto scale)
c
c----------------------------------------------------------------------
c
      module mspctr
      use nlsdim
      implicit none
c
      integer, save :: iscal(MXSITE), iscglb
      double precision, save ::
     #   spectr(MXPT,MXSITE), wspec(MXPT,MXSITE), sfac(MXSITE,MXSPC)
c
      end module mspctr

