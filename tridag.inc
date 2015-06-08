c   Version 1.5  5/2/94
c----------------------------------------------------------------------
c                    =========================
c                       common TRIDAG
c                    =========================
c
c  Storage for tridiagonal matrix used in Lanczos and conjugate gradients
c  routines
c
c  Note: stddim.inc should be included before this file
c----------------------------------------------------------------------
      double precision alpha,beta
      common /tridag/  alpha(2,MXSTEP),beta(2,MXSTEP)
