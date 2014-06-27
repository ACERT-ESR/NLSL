c NLSL Version 1.5 beta 11/23/95
c----------------------------------------------------------------------
c                    =========================
c                         common BASIS
c                    =========================
c
c  Contains common block for storage of basis set index lists
c  File "nlsdim.inc" should be included before this file 
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
      integer basno,ibasis,ixbas,ltbas,mts,nbas,nextbs,pruned,bsused
      character*30 basisID
c
      common /basis/ ibasis(5,MXDIM),basno(MXSITE,MXSPC),mts(9,MXMTS),
     #               ixbas(MXTDM),ltbas(MXTDM),bsused(MXTDM),
     #               nbas,nextbs,basisID(MXTDM)
