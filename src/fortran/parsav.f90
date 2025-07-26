! NLSL Version 1.9.0 beta 2/13/15
!----------------------------------------------------------------------
!                    =========================
!                         module PARSAV
!                    =========================
!
!  Defines variables to be used for saving a copy of the working
!  parameters and NLS search vector.
!
!  Note: this module uses the nlsdim module.
!----------------------------------------------------------------------
!
      module parsav
      use nlsdim
      implicit none
!
      double precision, save :: fxsave(6,MXVAR), fprsav(NFPRM,MXSITE),
     #                          spsave(4,MXSPC)
!
      integer, save :: ixsave(3,MXVAR), iprsav(NIPRM,MXSITE),
     #                 ixxsav(NFPRM,MXSITE), nstsav, nspsav, nprsav
!
      logical, save :: xsaved, prsaved
!
      character*9, save :: tagsav(MXVAR)
!
      end module parsav
