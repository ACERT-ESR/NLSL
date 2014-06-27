c       VERSION 1.0     12/17/88
c***********************************************************************
c
c       define the largest value of the L quantum number allowed 
c       in the basis set
c
c       Notes:
c               This file is included in the program lbll.f for the
c               purpose of verifying that the parameters input by the 
c               user are consistent with the arrays dimensioned in 
c               matrix element calculation subroutine, matrll.f, which 
c               also includes this file.
c
c       written by DJS 11-SEP-87
c
c***********************************************************************
c
        integer mxlval
        parameter (mxlval=190)
c
