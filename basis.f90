c NLSL Version 1.9.0 beta 2/9/15
c----------------------------------------------------------------------
c                    =========================
c                          module BASIS
c                    =========================
c
c  Contains common block for storage of basis set index lists.
c  Module nlsdim is used by this module.
c
c  ibasis:  block containing basis indices
c  basno:   number of block with basis set for a given site,spectrum
c  ixbas:   Index of basis set block used by each set
c  ltbas:   Length of basis set block used by each set
c  bsused:  Flag indicating whether set is in use
c  nbas:    Number of currently defined basis sets
c  nextbs:  Next available element of storage block
c  basisID: Filenames from which pruned basis sets were read
c----------------------------------------------------------------------
      module basis
      use nlsdim
      implicit none
c
      integer, save :: ibasis(5,MXDIM),basno(MXSITE,MXSPC),
     #                 mts(9,MXMTS),ixbas(MXTDM),ltbas(MXTDM),
     #                 bsused(MXTDM),nbas,nextbs
c
      character*30, save :: basisID(MXTDM)
c
      end module basis
