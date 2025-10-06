! NLSL Version 1.9.0 beta 2/12/15
!***********************************************************************
!                    =========================
!                           module MAXL
!                    =========================
!
!       Define the largest value of the L quantum number allowed 
!       in the basis set.
!
!       Notes:
!               This module is used by subroutine lcheck for the
!               purpose of verifying that the parameters input by the 
!               user are consistent with the arrays dimensioned in
!               the matrix element calculation subroutines. Those
!               calculation routines, which include matrll and pmatrl
!               as well as function w3j, also use this module.
!
!       written by DJS 11-SEP-87
!
!***********************************************************************
!
      module maxl
      implicit none
!
      integer, parameter :: mxlval=190
!
      end module maxl
