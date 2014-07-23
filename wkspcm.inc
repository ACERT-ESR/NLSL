c  VERSION 1.0  (NLSPMC version)   2/5/99
c*********************************************************************
c
c     Declarations of common block to hold work space
c
c     NOTE: the file 'nlsdim.inc' must be included before this file.
c
c*********************************************************************
c
	integer iwsp
	double precision rwsp
	complex*16 cwsp1,cwsp2,cwsp3,cwsp4,cwsp5,c2wsp1,c2wsp2
	common	/wkspcm/  c2wsp1(MXDIM,MXSTEP),c2wsp2(MXEGV,MXEGV),
     #          cwsp1(MXDIM),cwsp2(MXDIM),cwsp3(MXDIM),cwsp4(MXDIM),
     #          cwsp5(MXDIM),rwsp(MXDIM),iwsp(MXDIM)
