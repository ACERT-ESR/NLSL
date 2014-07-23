c  VERSION 1.0  (NLSPMC version)   2/5/99
c*********************************************************************
c                    =========================
c                      include file LIMITS
c                    =========================
c
c  Defines array dimensioning parameters for NLS programs and other
c  numerical values
c
c LINELG  Length of line read
c WORDLG  Length of words
c NAMELG  Length of filenames
c NFPRM   Number of floating point parameters for simulation
c NIPRM   Number of integer point parameters for simulation
c MXSPEC  Maximum number of experimental spectra to fit
c MXSITE  Maximum number of components of each spectra
c MXVAR   Maximum number of simultanious variables for least squares
c MXORT   Maximum MOMD orientations
c MXDIM   Maximum dimension of matrix
c MMXDIM  Max combined size of basis sets = MXDIM*MXSITE
c MXSTEP  Maximum number CG steps ??
c MXEL    Maximum matrix elements
c MXEGV   Maximum number of eigenvalues
c MXFILE  Maximum number of data files
c MXNPT   Maximum number of points along one of the axis
c         NOTE: THIS SHOULD BE A POWER OF 2
c MXSPT   Maximum number of points in an individual spectrum
c              (can be mxnpt*mxnpt maximum)
c MXPT    Maximum total number of data points
c              (can be MXSPT*MXSPEC*MXSITE maximum) - store all spec.
c MXSPH   Maximum number of tensors that may be specified in spherical
c         form

c
c*********************************************************************
c
      integer LINELG,WORDLG,NAMELG,NFPRM,NIPRM,MXSPEC,MXSITE,MXVAR,
     #	MXORT,MXDIM,MMXDIM,MXSTEP,NVPRM,
     #           MXEL,MXEGV,MXFILE,
     #           MXNPT,MXSPT,MXPT,
     #           MXSPH,NALIAS,NSYMTR,NSYMBL
c
c later remove namelg and set all to 40.  see call to gettkn
c added 5 site weights as variable parameters.  Don't vary them
c all! 
c mxdim was 800, 
c mxsite was 3
c mxort was 20
c mxegv was 1024 ywc
c MAXDIM was 1700
      parameter (LINELG=80,WORDLG=40,NAMELG=40,NFPRM=52,NIPRM=25,
     #	MXSPEC=10,MXSITE=2,MXVAR=15,MXORT=16,MXDIM=20500,
     #	MMXDIM=MXDIM*MXSITE,MXSTEP=1500,
     #  NVPRM=41,MXEL=1000000,MXEGV=256,MXFILE=20,
     #           MXNPT=128,MXSPT=16384,MXPT=163840,
     #           MXSPH=3,NALIAS=9,NSYMTR=3,NSYMBL=5)
c
c*********************************************************************
c
