c NLSL Version 1.3
c----------------------------------------------------------------------
c                    =========================
c                     include file TRIDAG.INC
c                    =========================
c     Common block for storage of Lanczos tridiagonal matrix and
c     calculated spectra. Also includes storage for individual spectra
c     calculated during the process of nonlinear least-squares fitting.
c
c
c     alpha : Array containing diagonal of tridiagonal matrix for all
c             spectral calculations
c     beta  : Array containing off-diagonal of tridiag. matrix for all
c             spectral calculations
c     ixtd  : Starting index in alpha and beta for each spectral
c             calculation
c     ltd   : Dimension of tridiagonal matrix for each spectral calculation
c
c     NOTE: file nlsdim.inc should be included before including this file.
c     David Budil 13 May 1992 Cornell University
c----------------------------------------------------------------------
c
      integer ixtd,ltd,modtd,tdsite,tdspec,nexttd,ntd
      double complex alpha,beta,stv,y
c
      common /tridag/ alpha(MXTDG),beta(MXTDG),stv(MXDIM),y(MXDIM),
     #                ixtd(MXSITE,MXSPC),ltd(MXSITE,MXSPC),
     #                modtd(MXSITE,MXSPC),tdspec(MXTDM),tdsite(MXTDM),
     #                nexttd,ntd


