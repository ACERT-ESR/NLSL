c      Version 1.5  5/4/94
c*********************************************************************
c
c	Declarations of common block to hold matrix for EPRLL
c
c       These declarations include the additional jzmat and kzmat
c       arrays used with the updated versions of MATRLL and SCMVM
c       subroutines.
c
c	Notes:
c		1) The dimensions of arrays used declared here are
c		   determined by parameters declared in	the include 
c		   file	stddim.inc, therefore the inclusion of the 
c		   file	must follow the	inclusion of stddim.inc.
c		2) For more information	on the usage of	these 
c		   arrays, see the header for the Lanczos algorithm 
c		   subroutine (cslnzs.f) and the sparse	matrix-vector
c		   multiplication routine (scmvm.f).
c
c	written	by DJS 11-SEP-87
c
c      NOTE: File 'stddim.inc' must be included before this file.
c*********************************************************************
c
	 integer izmat,jzmat,kzmat
	 double	precision zmat,zdiag
	 common	/scmat/	zmat(MXEL),zdiag(2,MXDIM),
     #		izmat(MXEL),jzmat(MXDIM1),kzmat(MXDIM1)
c
