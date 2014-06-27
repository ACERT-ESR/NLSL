c NLSL Version 1.3.2 2/27/94    
c*********************************************************************
c
c       Declarations of common block to hold matrix for EPRLL
c
c       These declarations include the additional jzmat and kzmat
c       arrays used with the updated versions of MATRLL and SCMVM
c       subroutines.
c
c       Notes:
c               1) The dimensions of arrays used declared here are
c                  determined by parameters declared in the include 
c                  file nlsdim.inc, therefore the inclusion of the 
c                  file must follow the inclusion of nlsdim.inc.
c
c               2) The diagonal elements of Z are stored in the complex 
c                  array zdiag in common /eprmat/. The upper diagonal of the 
c                  real and imaginary parts of Z are stored separately by row 
c                  in the zmat array (also in common /eprmat/), skipping
c                  zero elements. Imaginary elements are packed starting
c                  at the beginning of zmat, and real elements starting 
c                  at the upper limit of zmat. 
c                  The companion arrays jzmat and kzmat respectively
c                  give the the location of the first imaginary and real 
c                  elements of each row of Z within zmat, and izmat gives
c                  the column index corresponding to each element in zmat.
c
c         zdiag    : diagonal elements of Z
c         zmat     : real and imaginary elements of Z upper diagonal
c         jzmat(i) : index of 1st imaginary element of row i in zmat
c         kzmat(i) : index of 1st real element of row i in zmat
c         izmat(i) : Z column index of ith element of zmat
c
c       written by DJS 11-SEP-87
c
c*********************************************************************
c
         integer izmat,jzmat,kzmat
         double precision zmat,zdiag
         common /scmat/ zmat(MXEL),zdiag(2,MXDIM),
     #          izmat(MXEL),jzmat(MXDIM+1),kzmat(MXDIM+1)
c
