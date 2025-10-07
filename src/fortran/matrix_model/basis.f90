! NLSL Version 1.9.0 beta 2/9/15
!----------------------------------------------------------------------
!                    =========================
!                          module BASIS
!                    =========================
!
!  Contains common block for storage of basis set index lists.
!  Module nlsdim is used by this module.
!
!  ibasis:  block containing basis indices
!  basno:   number of block with basis set for a given site,spectrum
!  ixbas:   Index of basis set block used by each set
!  ltbas:   Length of basis set block used by each set
!  bsused:  Flag indicating whether set is in use
!  nbas:    Number of currently defined basis sets
!  nextbs:  Next available element of storage block
!  basisID: Filenames from which pruned basis sets were read
!----------------------------------------------------------------------
      module basis
      use nlsdim
      implicit none
!
      integer, save :: ibasis(5,MXDIM),basno(MXSITE,MXSPC),
     #                 mts(9,MXMTS),ixbas(MXTDM),ltbas(MXTDM),
     #                 bsused(MXTDM),nbas,nextbs
!
      character*30, save :: basisID(MXTDM)
!
      end module basis
