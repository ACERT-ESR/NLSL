c     Version 1.5  5/2/94
c----------------------------------------------------------------------
c                    =========================
c                       common block INDEXL
c                    =========================
c
c  Common block containining a list of the basis set indices
c  for each element of the starting vector for the EPRLL family
c  of programs. The elements of this block should be initialized
c  by a call to subroutine LBASIX.
c
c  NOTE: File 'stddim.inc' must be included before this file
c
c----------------------------------------------------------------------
c
      integer l1,k1,m1,pi1,qi1
      common /indexl/ l1(MXDIM),k1(MXDIM),m1(MXDIM),pi1(MXDIM),
     #                qi1(MXDIM)









