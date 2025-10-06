! NLSL Version 1.9.0 beta 2/12/15
!----------------------------------------------------------------------
!                    =========================
!                          module MSPCTR
!                    =========================
! Defines working arrays used for multiple-site spectral fitting
!
! Note: module nlsdim is used by this module
!
!    spectr    Storage for calculated spectra for each site, spectrum
!
!    wspec     Work array for temporary storage of calculated spectra
!              Used by sshift to store correlation function for each
!              spectrum  
!              Used before call to sscale (array is replaced by QR
!              decomposition)
!
!    sfac      Scale factors for each site and spectrum
!
!    iscal     Flag indicating whether each site is to be automatically
!              scaled in the function calculation (1=auto scale)
!
!    iscglb    Flag indicating whether all sites are to be automatically
!              scaled (1=auto scale)
!
!----------------------------------------------------------------------
!
      module mspctr
      use nlsdim
      implicit none
!
      integer, save :: iscal(MXSITE), iscglb
      double precision, target, save :: spectr(MXPT,MXSITE)
      double precision, target, save :: wspec(MXPT,MXSITE)
      double precision, target, save :: sfac(MXSITE,MXSPC)
!
      end module mspctr

