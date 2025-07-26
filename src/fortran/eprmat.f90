! NLSL Version 1.9.0 beta 2/12/15
!*********************************************************************
!                    =========================
!                          module EPRMAT
!                    =========================
!
!       Declarations of Z matrix for EPRLL
!
!       These declarations include the additional jzmat and kzmat
!       arrays used with the updated versions of MATRLL and SCMVM
!       subroutines.
!
!       Notes:
!               1) The dimensions of the arrays declared here are
!                  determined by parameters declared in the module 
!                  nlsdim, which is used below. 
!
!               2) The diagonal elements of Z are stored in the
!                  complex array zdiag. The upper diagonal of the 
!                  real and imaginary parts of Z are stored
!                  separately by row in the zmat array, skipping
!                  zero elements. Imaginary elements are packed
!                  starting at the beginning of zmat, and real
!                  elements starting at the upper limit of zmat. 
!                  Companion arrays jzmat and kzmat respectively
!                  give the the locations of the first imaginary
!                  and real elements of each row of Z within zmat,
!                  and izmat gives the column index corresponding
!                  to each element in zmat.
!
!               3) This was originally a common block named scmat,
!                  but all dependent subprograms simply included
!                  file eprmat.inc, so the module is named eprmat.
!
!         zdiag    : diagonal elements of Z
!         zmat     : real and imaginary elements of Z upper diagonal
!         jzmat(i) : index of 1st imaginary element of row i in zmat
!         kzmat(i) : index of 1st real element of row i in zmat
!         izmat(i) : Z column index of ith element of zmat
!
!       written by DJS 11-SEP-87
!
!*********************************************************************
!
      module eprmat
      use nlsdim
      implicit none
!
      double precision, save :: zmat(MXEL), zdiag(2,MXDIM)
      integer, save :: izmat(MXEL), jzmat(MXDIM+1), kzmat(MXDIM+1)
!
      end module eprmat
