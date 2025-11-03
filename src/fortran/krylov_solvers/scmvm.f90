! NLSL Version 1.3.2
!**********************************************************************
!
!               (S)parse (C)omplex (M)atrix (V)ector (M)ultiply
!               ===============================================
!
!       This subroutine will do a sparse matrix-vector multiplication
!       of the form y=Z*x where Z is sparse and complex. 
!
!       NB: This routine utilizes the new matrix storage convention
!       implemented by DJS in which only the upper half of the
!       matrix is stored. It should be used only with versions
!       of the MATRLL routine that store the matrix according to
!       this convention (described below). 
!
!       The diagonal elements of Z are stored in the complex 
!       array zdiag in common /eprmat/. The upper diagonal of the 
!       real and imaginary parts of Z are stored separately by row 
!       in the zmat array (also in common /eprmat/), skipping
!       zero elements. Imaginary elements are packed starting
!       at the beginning of zmat, and real elements starting 
!       at the upper limit of zmat. 
!       The companion arrays jzmat and kzmat respectively
!       give the the location of the first imaginary and real 
!       elements of each row of Z within zmat, and izmat gives
!       the column index corresponding to each element in zmat.
!
!       The matrix-vector multiply will not be performed if the 
!       user-halt flags defined in stdio.inc have been set 
!
!       Calling parameters:
!       ------------------
!
!        vectors :
!
!         in /eprmat/
!
!         zdiag    : diagonal elements of Z
!         zmat     : real and imaginary elements of Z upper diagonal
!         jzmat(i) : index of 1st imaginary element of row i in zmat
!         kzmat(i) : index of 1st real element of row i in zmat
!         izmat(i) : Z column index of ith element of zmat
!
!         arguments:
!
!         x        : input vector for matrix-vector multiplication
!         y        : resultant vector 
!
!       scalars :
!
!             ndim  : number of rows in matrix
!
!       Local Variables:
!       ---------------
!
!             accr  : real part of dot product of a row of the 
!                     matrix with the input vector
!             acci  : imaginary part of dot product of a row of the 
!                     matrix with the input vector
!
!       Notes:
!       -----
!               The routine does not use complex double precision
!               arithmetic explicitly for two reasons.  First is
!               that, unfortunately, not all F77 compilers support
!               complex double precision arithmetic.  Second, most of
!               floating point operations involve the multiplication
!               of a purely real or purely imaginary number by a
!               complex number.  This is more efficiently done by
!               treating the complex number an ordered pair or real
!               numbers, since there is no unnecessary multiplications
!               by zero performed.
!
!
!       Includes:
!               nlsdim.inc
!               rndoff.inc
!               eprmat.inc
!
!       Uses: 
!
!       written by DJS 3-OCT-86
!
!**********************************************************************
!
      subroutine scmvm(x,y,ndim)
!
      use nlsdim
      use rnddbl
      use eprmat
!
      integer ndim
      double precision x,y    
      dimension x(2,MXDIM),y(2,MXDIM)
!
      integer j,k,m,n,n1
      double precision accr,acci
!
!######################################################################
!
!----------------------------------------------------------------------
!     do diagonal elements first 
!----------------------------------------------------------------------
!
      do n=1,ndim 
        y(1,n)=zdiag(1,n)*x(1,n)-zdiag(2,n)*x(2,n)
        y(2,n)=zdiag(1,n)*x(2,n)+zdiag(2,n)*x(1,n)
      end do
!
!----------------------------------------------------------------------
!     loop over rows (columns) of matrix for off-diagonal elements
!----------------------------------------------------------------------
!
      do n=1,ndim
        n1=n+1
!
        accr=0.0D0
        acci=0.0D0
!
!       imaginary matrix elements
!
        if (jzmat(n) .ne. jzmat(n1) ) then
           do j=jzmat(n),jzmat(n1)-1
              m=izmat(j)
              acci=acci+zmat(j)*x(1,m)
              y(2,m)=y(2,m)+zmat(j)*x(1,n) 
              accr=accr-zmat(j)*x(2,m)
              y(1,m)=y(1,m)-zmat(j)*x(2,n)
           end do
        endif
!
!       real matrix elements
!
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
!
        y(1,n)=y(1,n)+accr
!*djs        if (abs(y(1,n)).lt.rndoff) y(1,n)=0.0D0
        y(2,n)=y(2,n)+acci
!*djs        if (abs(y(2,n)).lt.rndoff) y(2,n)=0.0D0
!
      end do
!
      return
      end
