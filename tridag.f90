c NLSL Version 1.9.0 beta 2/13/15
c----------------------------------------------------------------------
c                    =========================
c                          module TRIDAG
c                    =========================
c
c     Defines arrays for storage of Lanczos tridiagonal matrix and
c     calculated spectra. Also includes storage for individual spectra
c     calculated during the process of nonlinear least-squares fitting.
c
c     alpha : Array containing diagonal of tridiagonal matrix for all
c             spectral calculations
c     beta  : Array containing off-diagonal of tridiag. matrix for all
c             spectral calculations
c     ixtd  : Starting index in alpha and beta for each spectral
c             calculation
c     ltd   : Dimension of tridiagonal matrix for each spectral calculation
c
c     NOTE: This module uses the nlsdim module.
c
c     David Budil 13 May 1992 Cornell University
c----------------------------------------------------------------------
c
      module tridag
      use nlsdim
      implicit none
c
      integer, save :: ixtd(MXSITE,MXSPC), ltd(MXSITE,MXSPC),
     #                 modtd(MXSITE,MXSPC), tdspec(MXTDM),
     #                 tdsite(MXTDM), nexttd, ntd
c
      double complex, save :: alpha(MXTDG), beta(MXTDG),
     #                        stv(MXDIM),y(MXDIM)
c
      end module tridag

