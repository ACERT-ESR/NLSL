c  Version 1.5 5/2/94
c***********************************************************************
c
c	Define the largest value of the	L quantum number allowed 
c	in the basis set in the presence of a potential
c
c	Notes:
c		This file is included in the program lbll.f for	the
c		purpose	of verifying that the parameters input by the 
c		user are consistent with the arrays dimensioned	in 
c		matrix element calculation subroutine, matrll.f, which 
c		also includes this file.
c
c               In version 1.2, this file is also used by the bincof
c               routine to define the maximum order of the binomial
c               coefficients used to calculate Wigner 3-j symbols
c               in the asymptotic high-j limit (function wig3j in file
c               w3j.f)
c
c	written	by DJS 11-SEP-87
c       modified by DEB 22-May-92
c
c***********************************************************************
c
	integer	MXLVAL
	parameter (MXLVAL=80)
