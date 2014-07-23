c  VERSION 1.0  (NLSPMC version)   2/5/99
c*********************************************************************
c
c     Declarations of common block to hold 2D spectra
c
c     NOTE: the file 'nlsdim.inc' must be included before this file.
c  Add MXSITE for multi components.  Note, multiple spectra are 
c  indicated by MXPT=10*size of one spectrum.
c
c*********************************************************************
c
	double precision spectr,totspec
	complex*16 xspectr,cspectr,ctotspec
	common	/tdspec/  spectr(MXPT,MXSITE),
     #		xspectr(MXPT,MXSITE),totspec(MXPT),
     #		cspectr(MXPT,MXSITE),ctotspec(MXPT)
