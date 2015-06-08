c     Version 1.5  5/2/94
c----------------------------------------------------------------------
c                    =========================
c                      common block CGDATA
c                    =========================
c
c   Storage for temporary vectors used in CG algorithm.
c
c   NOTE: file 'stddim.inc' must be included before this file
c----------------------------------------------------------------------
      complex*16 p,r,w,z
      common /cgdata/ p(MXDIM),r(MXDIM),w(MXDIM),z(MXDIM)

