c NLSL Version 1.9.0 beta 2/4/15
c*********************************************************************
c                    =========================
c                          module NLSDIM
c                    =========================
c
c  Defines array dimensioning parameters for NLS programs
c
c MXSTEP  Maximum number of CG/Lanczos steps
c MXDIM   Maximum dimension of matrix
c MXEL    Maximum number of off-diagonal elements in matrix
c MXINP   Maximum number of points in an input file
c MXSPC   Maximum number of spectra that can be fit at once
c MXSITE  Maximum number of sites that can be fit for a given spectrum
c MXSPT   Maximum number of points in an individual spectrum
c         NOTE: THIS SHOULD BE A POWER OF 2
c MXPT    Maximum total number of data points
c MXVAR   Maximum number of parameters that may be varied
c MXJCOL  Maximum number of columns in Jacobian matrix for L-M least
c         squares (should be MXVAR+MXSITE)
c MXTV    Maximum number of variables in "transformed" covariance
c         matrix
c MXCMT   Maximum number of comment lines saved from a datafile
c MXTDG   Maximum total length of tridiagonal matrix element arrays
c         (upper limit needed is MXSTEP*MXSPC*MXSITE)
c MXTDM   Maximum number of blocks in tridiagonal matrix allocation
c         scheme (should be MXSITE*MXSPC) 
c MXSPH   Maximum number of tensors quantities in fprm array
c MXMTS   Dimension of mts array to save truncation indices/flags
c MXFILE  Maximum number of script files that may be open
c NSYMTR  Number of different symmetry types currently defined
c NSYMBL  Number of symbols currently defined
c*********************************************************************
c

      module nlsdim
      implicit none
c
      integer, parameter :: MXSTEP=2000,
     #           MXDIM=45000,MXDIM1=MXDIM+1,
     #           MXEL=50000000,
     #           MXSPC=4,
     #           MXCMT=16,
     #           MXSPT=512,
     #           MXPT=MXSPT*MXSPC,
     #           MXVAR=10,
     #           MXINP=4096,
     #           MXSITE=3,
     #           MXSPH=4,
     #           MXFILE=4,
     #           MXTDG=MXSTEP*MXSITE*MXSPC,
     #           MXTDM=MXSPC*MXSITE,
     #           MXJCOL=MXVAR+MXSITE,
     #           MXTV=MXJCOL+4*MXSITE,
     #           MXMTS=13
c
      integer, parameter :: NFPRM=43,
     #           NVPRM=35,
     #           NIPRM=24,
     #           NALIAS=12,
     #           NSYMTR=3,
     #           NSYMBL=5
c
c*********************************************************************
c
      end module nlsdim
