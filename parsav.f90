c NLSL Version 1.9.0 beta 2/13/15
c----------------------------------------------------------------------
c                    =========================
c                         module PARSAV
c                    =========================
c
c  Defines variables to be used for saving a copy of the working
c  parameters and NLS search vector.
c
c  Note: this module uses the nlsdim module.
c----------------------------------------------------------------------
c
      module parsav
      use nlsdim
      implicit none
c
      double precision, save :: fxsave(6,MXVAR), fprsav(NFPRM,MXSITE),
     #                          spsave(4,MXSPC)
c
      integer, save :: ixsave(3,MXVAR), iprsav(NIPRM,MXSITE),
     #                 ixxsav(NFPRM,MXSITE), nstsav, nspsav, nprsav
c
      logical, save :: xsaved, prsaved
c
      character*9, save :: tagsav(MXVAR)
c
      end module parsav
