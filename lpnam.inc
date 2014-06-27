c NLSL Version 1.3
c----------------------------------------------------------------------
c                    =========================
c                       common block LPNAM
c                    =========================
c
c   Definition of common block containing names of parameters,
c   aliases, and various symbols
c
c NOTE: the file nlsdim.inc must be included before this file
c----------------------------------------------------------------------
      character*6 parnam,alias1,alias2,iprnam
      character*10 symstr,symbol
      integer symval
      common /lpnam/ symval(NSYMBL),parnam(NFPRM),alias1(NALIAS),
     #  alias2(NALIAS),iprnam(NIPRM),symstr(NSYMTR),symbol(NSYMBL)
