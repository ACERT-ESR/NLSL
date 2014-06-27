c Version 1.4 10/10/94
c----------------------------------------------------------------------
c                    =========================
c                         common PARSAV
c                    =========================
c
c  Defines common block to be used for saving a copy of the working
c  parameters and NLS search vector.
c----------------------------------------------------------------------
c
      integer ixsave,iprsav,ixxsav,nstsav,nsesav,nprsav
      double precision fxsave,fprsav,spsave
      character*9 tagsave
      logical xsaved,prsaved
c
      common /parsav/ fxsave(6,MXVAR),fprsav(NFPRM,MXSITE),
     #                spsave(4,MXSPC),
     #                ixsave(3,MXVAR),iprsav(NIPRM,MXSITE),
     #                ixxsav(NFPRM,MXSITE),nstsav,nspsav,nprsav,
     #                xsaved,prsaved,tagsav(MXVAR)
c
