c NLSL Version 1.9.0 beta 2/12/15
c***********************************************************************
c                    =========================
c                           module MAXL
c                    =========================
c
c       Define the largest value of the L quantum number allowed 
c       in the basis set.
c
c       Notes:
c               This module is used by subroutine lcheck for the
c               purpose of verifying that the parameters input by the 
c               user are consistent with the arrays dimensioned in
c               the matrix element calculation subroutines. Those
c               calculation routines, which include matrll and pmatrl
c               as well as function w3j, also use this module.
c
c       written by DJS 11-SEP-87
c
c***********************************************************************
c
      module maxl
      implicit none
c
      integer, parameter :: mxlval=190
c
      end module maxl
