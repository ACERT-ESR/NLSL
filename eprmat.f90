c NLSL Version 1.9.0 beta 2/12/15
c*********************************************************************
c                    =========================
c                          module EPRMAT
c                    =========================
c
c       Declarations of Z matrix for EPRLL
c
c       These declarations include the additional jzmat and kzmat
c       arrays used with the updated versions of MATRLL and SCMVM
c       subroutines.
c
c       Notes:
c               1) The dimensions of the arrays declared here are
c                  determined by parameters declared in the module 
c                  nlsdim, which is used below. 
c
c               2) The diagonal elements of Z are stored in the
c                  complex array zdiag. The upper diagonal of the 
c                  real and imaginary parts of Z are stored
c                  separately by row in the zmat array, skipping
c                  zero elements. Imaginary elements are packed
c                  starting at the beginning of zmat, and real
c                  elements starting at the upper limit of zmat. 
c                  Companion arrays jzmat and kzmat respectively
c                  give the the locations of the first imaginary
c                  and real elements of each row of Z within zmat,
c                  and izmat gives the column index corresponding
c                  to each element in zmat.
c
c               3) This was originally a common block named scmat,
c                  but all dependent subprograms simply included
c                  file eprmat.inc, so the module is named eprmat.
c
c         zdiag    : diagonal elements of Z
c         zmat     : real and imaginary elements of Z upper diagonal
c         jzmat(i) : index of 1st imaginary element of row i in zmat
c         kzmat(i) : index of 1st real element of row i in zmat
c         izmat(i) : Z column index of ith element of zmat
c
c       written by DJS 11-SEP-87
c
c*********************************************************************
c
      module eprmat
      use nlsdim
      implicit none
c
      double precision, save :: zmat(MXEL), zdiag(2,MXDIM)
      integer, save :: izmat(MXEL), jzmat(MXDIM+1), kzmat(MXDIM+1)
c
      end module eprmat
