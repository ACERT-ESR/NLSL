c NLSL Version 1.3.2
c**********************************************************************
c
c               (S)parse (C)omplex (M)atrix (V)ector (M)ultiply
c               ===============================================
c
c       This subroutine will do a sparse matrix-vector multiplication
c       of the form y=Z*x where Z is sparse and complex. 
c
c       NB: This routine utilizes the new matrix storage convention
c       implemented by DJS in which only the upper half of the
c       matrix is stored. It should be used only with versions
c       of the MATRLL routine that store the matrix according to
c       this convention (described below). 
c
c       The diagonal elements of Z are stored in the complex 
c       array zdiag in common /eprmat/. The upper diagonal of the 
c       real and imaginary parts of Z are stored separately by row 
c       in the zmat array (also in common /eprmat/), skipping
c       zero elements. Imaginary elements are packed starting
c       at the beginning of zmat, and real elements starting 
c       at the upper limit of zmat. 
c       The companion arrays jzmat and kzmat respectively
c       give the the location of the first imaginary and real 
c       elements of each row of Z within zmat, and izmat gives
c       the column index corresponding to each element in zmat.
c
c       The matrix-vector multiply will not be performed if the 
c       user-halt flags defined in stdio.inc have been set 
c
c       Calling parameters:
c       ------------------
c
c        vectors :
c
c         in /eprmat/
c
c         zdiag    : diagonal elements of Z
c         zmat     : real and imaginary elements of Z upper diagonal
c         jzmat(i) : index of 1st imaginary element of row i in zmat
c         kzmat(i) : index of 1st real element of row i in zmat
c         izmat(i) : Z column index of ith element of zmat
c
c         arguments:
c
c         x        : input vector for matrix-vector multiplication
c         y        : resultant vector 
c
c       scalars :
c
c             ndim  : number of rows in matrix
c
c       Local Variables:
c       ---------------
c
c             accr  : real part of dot product of a row of the 
c                     matrix with the input vector
c             acci  : imaginary part of dot product of a row of the 
c                     matrix with the input vector
c
c       Notes:
c       -----
c               The routine does not use complex double precision
c               arithmetic explicitly for two reasons.  First is
c               that, unfortunately, not all F77 compilers support
c               complex double precision arithmetic.  Second, most of
c               floating point operations involve the multiplication
c               of a purely real or purely imaginary number by a
c               complex number.  This is more efficiently done by
c               treating the complex number an ordered pair or real
c               numbers, since there is no unnecessary multiplications
c               by zero performed.
c
c
c       Includes:
c               nlsdim.inc
c               rndoff.inc
c               eprmat.inc
c
c       Uses: 
c
c       written by DJS 3-OCT-86
c
c**********************************************************************
c
      subroutine scmvm(x,y,ndim)
c
      use nlsdim
      use rnddbl
      use eprmat
c
      integer ndim
      double precision x,y    
      dimension x(2,MXDIM),y(2,MXDIM)
c
      integer j,k,m,n,n1
      double precision accr,acci
c
c######################################################################
c
c----------------------------------------------------------------------
c     do diagonal elements first 
c----------------------------------------------------------------------
c
      do n=1,ndim 
        y(1,n)=zdiag(1,n)*x(1,n)-zdiag(2,n)*x(2,n)
        y(2,n)=zdiag(1,n)*x(2,n)+zdiag(2,n)*x(1,n)
      end do
c
c----------------------------------------------------------------------
c     loop over rows (columns) of matrix for off-diagonal elements
c----------------------------------------------------------------------
c
      do n=1,ndim
        n1=n+1
c
        accr=0.0D0
        acci=0.0D0
c
c       imaginary matrix elements
c
        if (jzmat(n) .ne. jzmat(n1) ) then
           do j=jzmat(n),jzmat(n1)-1
              m=izmat(j)
              acci=acci+zmat(j)*x(1,m)
              y(2,m)=y(2,m)+zmat(j)*x(1,n) 
              accr=accr-zmat(j)*x(2,m)
              y(1,m)=y(1,m)-zmat(j)*x(2,n)
           end do
        endif
c
c       real matrix elements
c
        if (kzmat(n) .ne. kzmat(n1)) then
           do k=kzmat(n),kzmat(n1)-1
              j = mxel-k+1
              m=izmat(j)
              accr=accr+zmat(j)*x(1,m)
              y(1,m)=y(1,m)+zmat(j)*x(1,n)          
              acci=acci+zmat(j)*x(2,m)
              y(2,m)=y(2,m)+zmat(j)*x(2,n)
           end do
        endif
c
        y(1,n)=y(1,n)+accr
c*djs        if (abs(y(1,n)).lt.rndoff) y(1,n)=0.0D0
        y(2,n)=y(2,n)+acci
c*djs        if (abs(y(2,n)).lt.rndoff) y(2,n)=0.0D0
c
      end do
c
      return
      end
