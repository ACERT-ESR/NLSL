c  VERSION 1.0  (NLSPMC version)   2/5/99
c*********************************************************************
c
c     Declarations of common block to hold starting vector for off-
c     diagonal & diagonal space.
c
c     NOTE: the file 'nlsdim.inc' must be included before this file.
c
c*********************************************************************
c
	 complex*16 stvo,stvd
	 common	/stvcom/  stvo(mxdim),stvd(mxdim,5)
