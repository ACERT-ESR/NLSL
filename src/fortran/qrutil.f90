c  NLSL Version 1.3.2
c  This file contains subroutines QTBVEC and RSOLVE
c----------------------------------------------------------------------
c                    =========================
c                         subroutine QTBVEC   
c                    =========================
c
c    This subroutine is one of the steps necessary to complete the
c    solution of a set of linear equations A*x=b using the QR decomposition
c    of the matrix A provided by subroutine QRFAC form the MINPACK library:
c             A*P = Q*R
c    where R is upper triangular, Q is an orthogonal matrix, and P
c    is a permutation matrix that accounts for any column pivoting that
c    occurred during the solution of A. Given the vector b, this routine
c    calculates the vector Q(transpose)*b that may be subsequently be
c    used to solve the original equation for x by solving the triangular
c    matrix equation
c
c    Input parameters:
c
c      M:     Number of rows in Q
c
c      N:     Number of columns in Q (and dimension of R,RDIAG,and QTB)
c
c      Q:     The strict upper triangle of Q is assumed to contain
c             the strict upper triangle of R and the lower triangle of Q
c             contains a factored form of Q, (as returned by QRFAC). 
c
c      LDQ:   Leading dimension of the Q array
c
c      RDIAG: Diagonal elements of R (also returned by QRFAC)
c
c      B:     Contains M-vector b.
c
c      QTB:   On output, the first N elements contain the N-vector
c             Q(transpose)*b product. 
c----------------------------------------------------------------------
      subroutine qtbvec(m,n,q,ldq,qraux,b,qtb)
      implicit none
      integer ldq,n,m
      double precision q(ldq,n),qraux(n),b(m),qtb(m)
c
      integer i,j
      double precision qtemp,sum,temp,ZERO
      parameter (ZERO=0.0d0)
c
      do i=1,m
         qtb(i)=b(i)
      end do
c
      do j=1,n
         qtemp=q(j,j)
         if (qraux(j).ne.ZERO) then
            q(j,j)=qraux(j)
            sum=ZERO
            do i=j,m
               sum=sum+q(i,j)*qtb(i)
            end do
            temp=-sum/q(j,j)
            do i=j,m
               qtb(i)=qtb(i)+q(i,j)*temp
            end do
         end if
c     
         q(j,j)=qtemp
      end do
c     
      return
      end


c----------------------------------------------------------------------
c                      =========================
c                          subroutine RSOLVE
c                      =========================
c
c    This subroutine is one of the steps necessary to complete the
c    solution of a set of linear equations A*x=b using the QR decomposition
c    of the matrix A provided by subroutine QRFAC form the MINPACK library:
c
c             A*P = Q*R
c
c    where R is upper triangular, Q is an orthogonal matrix, and P
c    is a permutation matrix that accounts for any column pivoting that
c    occurred during the solution of A. Given the vectors x and Qtb
c    [Q(transpose)*b], this routine solves the diagonal equation R*x=Qt*b
c    used to solve the original equation for x by solving the triangular
c    matrix equation
c
c    Input parameters:
c
c      M:     Number of rows in Q
c
c      N:     Number of columns in Q (and dimension of R,RDIAG,and QTB)
c
c      Q:     The strict upper triangle of Q is assumed to contain
c             the strict upper triangle of R and the lower triangle of Q
c             contains a factored form of Q, (as returned by QRFAC). 
c
c      LDQ:   Leading dimension of the Q array
c
c      RDIAG: Diagonal elements of R (also returned by QRFAC)
c
c      QTB:   On input, the first N elements contain the N-vector
c             Q(transpose)*b product. 
c
c      X:     On output, an N-vector containing solution of A*x=b
c
c      RSD:   On input, the M-vector b for which the solution/residuals
c             are desired
c
c             On output, an M-vector containing the residuals of the
c             least-squares problem, b-A*x.
c
c----------------------------------------------------------------------
      subroutine rsolve( m,n,q,ldq,qraux,qtb,x,rcalc,rsd ) 
      implicit none
      integer m,n,ldq
      double precision q(ldq,n),qraux(n),qtb(n),x(m),rsd(m)
      logical rcalc
c
      integer i,j,ju,jj
      double precision qtemp,sum,temp
c
      double precision ZERO
      parameter(ZERO=0.0d0)
c
c######################################################################
c
      do i=1,n
         x(i)=qtb(i)
      end do
c
c  --- Solve R*x = Q(transpose)*b
c
      do jj=1,n
         j=n-jj+1
         if (q(j,j).eq.ZERO) go to 4
         x(j)=x(j)/q(j,j)
         if (j.gt.1) then
            temp=-x(j)
            do i=1,j-1
               x(i)=temp*q(i,j)+x(i)
            end do
         end if
         end do
c
c  --- If required, calculate residual vector b-A*x 
c
 4    if (rcalc) then
         do i=1,m
            if (i.le.n) then
               rsd(i)=ZERO
            else
               rsd(i)=qtb(i)
            end if
         end do
c
         ju=min0(n,m-1)
         do jj=1,ju
            j=ju-jj+1
            if (qraux(j).ne.ZERO) then
               qtemp=q(j,j)
               q(j,j)=qraux(j)
c
               sum=ZERO
               do i=j,m
                  sum=sum+q(i,j)*rsd(j)
               end do
               temp=-sum/q(j,j)
c
               do i=j,m
                  rsd(i)=temp*q(i,j)+rsd(i)
               end do
               q(j,j)=qtemp
            end if
         end do
      end if
c
      return
      end
