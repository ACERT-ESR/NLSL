! NLSL Version 1.9.0 beta 2/4/15
!*********************************************************************
!                    =========================
!                          module NLSDIM
!                    =========================
!
!  Defines array dimensioning parameters for NLS programs
!
! MXSTEP  Maximum number of CG/Lanczos steps
! MXDIM   Maximum dimension of matrix
! MXEL    Maximum number of off-diagonal elements in matrix
! MXINP   Maximum number of points in an input file
! MXSPC   Maximum number of spectra that can be fit at once
! MXSITE  Maximum number of sites that can be fit for a given spectrum
! MXSPT   Maximum number of points in an individual spectrum
!         NOTE: THIS SHOULD BE A POWER OF 2
! MXPT    Maximum total number of data points
! MXVAR   Maximum number of parameters that may be varied
! MXJCOL  Maximum number of columns in Jacobian matrix for L-M least
!         squares (should be MXVAR+MXSITE+MXSPC)
! MXTV    Maximum number of variables in "transformed" covariance
!         matrix
! MXCMT   Maximum number of comment lines saved from a datafile
! MXTDG   Maximum total length of tridiagonal matrix element arrays
!         (upper limit needed is MXSTEP*MXSPC*MXSITE)
! MXTDM   Maximum number of blocks in tridiagonal matrix allocation
!         scheme (should be MXSITE*MXSPC) 
! MXSPH   Maximum number of tensors quantities in fprm array
! MXMTS   Dimension of mts array to save truncation indices/flags
! MXFILE  Maximum number of script files that may be open
! NSYMTR  Number of different symmetry types currently defined
! NSYMBL  Number of symbols currently defined
!*********************************************************************
!

      module nlsdim
      implicit none
!
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
     #           MXJCOL=MXVAR+MXSITE+MXSPC,
     #           MXTV=MXJCOL+4*MXSITE,
     #           MXMTS=13
!
      integer, parameter :: NFPRM=43,
     #           NFLMPR=4,
     #           NILMPR=4,
     #           NVPRM=35,
     #           NIPRM=24,
     #           NALIAS=12,
     #           NSYMTR=3,
     #           NSYMBL=5
!
!*********************************************************************
!
      end module nlsdim
