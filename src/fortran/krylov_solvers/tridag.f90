! NLSL Version 1.9.0 beta 2/13/15
!----------------------------------------------------------------------
!                    =========================
!                          module TRIDAG
!                    =========================
!
!     Defines arrays for storage of Lanczos tridiagonal matrix and
!     calculated spectra. Also includes storage for individual spectra
!     calculated during the process of nonlinear least-squares fitting.
!
!     alpha : Array containing diagonal of tridiagonal matrix for all
!             spectral calculations
!     beta  : Array containing off-diagonal of tridiag. matrix for all
!             spectral calculations
!     ixtd  : Starting index in alpha and beta for each spectral
!             calculation
!     ltd   : Dimension of tridiagonal matrix for each spectral calculation
!
!     NOTE: This module uses the nlsdim module.
!
!     David Budil 13 May 1992 Cornell University
!----------------------------------------------------------------------
!
      module tridag
      use nlsdim
      implicit none
!
      integer, save :: ixtd(MXSITE,MXSPC), ltd(MXSITE,MXSPC),
     #                 modtd(MXSITE,MXSPC), tdspec(MXTDM),
     #                 tdsite(MXTDM), nexttd, ntd
!
      double complex, save :: alpha(MXTDG), beta(MXTDG),
     #                        stv(MXDIM),y(MXDIM)
!
      end module tridag

