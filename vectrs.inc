c    Version 1.5   5/2/94
c----------------------------------------------------------------------
c                    =========================
c                         common VECTRS
c                    =========================
c
c   Storage for vectors used in Lanczos and conjugate gradients routines
c
c   Note: File stddim.inc must be included before this file.
c----------------------------------------------------------------------
      double precision x,y
      common /vectrs/ x(2,MXDIM),y(2,MXDIM)
