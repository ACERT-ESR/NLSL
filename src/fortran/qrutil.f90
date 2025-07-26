!  NLSL Version 1.3.2
!  This file contains subroutines QTBVEC and RSOLVE
!----------------------------------------------------------------------
!                    =========================
!                         subroutine QTBVEC   
!                    =========================
!
!    This subroutine is one of the steps necessary to complete the
!    solution of a set of linear equations A*x=b using the QR decomposition
!    of the matrix A provided by subroutine QRFAC form the MINPACK library:
!             A*P = Q*R
!    where R is upper triangular, Q is an orthogonal matrix, and P
!    is a permutation matrix that accounts for any column pivoting that
!    occurred during the solution of A. Given the vector b, this routine
!    calculates the vector Q(transpose)*b that may be subsequently be
!    used to solve the original equation for x by solving the triangular
!    matrix equation
!
!    Input parameters:
!
!      M:     Number of rows in Q
!
!      N:     Number of columns in Q (and dimension of R,RDIAG,and QTB)
!
!      Q:     The strict upper triangle of Q is assumed to contain
!             the strict upper triangle of R and the lower triangle of Q
!             contains a factored form of Q, (as returned by QRFAC). 
!
!      LDQ:   Leading dimension of the Q array
!
!      RDIAG: Diagonal elements of R (also returned by QRFAC)
!
!      B:     Contains M-vector b.
!
!      QTB:   On output, the first N elements contain the N-vector
!             Q(transpose)*b product. 
!----------------------------------------------------------------------
      subroutine qtbvec(m,n,q,ldq,qraux,b,qtb)
      implicit none
      integer ldq,n,m
      double precision q(ldq,n),qraux(n),b(m),qtb(m)
!
      integer i,j
      double precision qtemp,sum,temp,ZERO
      parameter (ZERO=0.0d0)
!
      do i=1,m
         qtb(i)=b(i)
      end do
!
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
!     
         q(j,j)=qtemp
      end do
!     
      return
      end


!----------------------------------------------------------------------
!                      =========================
!                          subroutine RSOLVE
!                      =========================
!
!    This subroutine is one of the steps necessary to complete the
!    solution of a set of linear equations A*x=b using the QR decomposition
!    of the matrix A provided by subroutine QRFAC form the MINPACK library:
!
!             A*P = Q*R
!
!    where R is upper triangular, Q is an orthogonal matrix, and P
!    is a permutation matrix that accounts for any column pivoting that
!    occurred during the solution of A. Given the vectors x and Qtb
!    [Q(transpose)*b], this routine solves the diagonal equation R*x=Qt*b
!    used to solve the original equation for x by solving the triangular
!    matrix equation
!
!    Input parameters:
!
!      M:     Number of rows in Q
!
!      N:     Number of columns in Q (and dimension of R,RDIAG,and QTB)
!
!      Q:     The strict upper triangle of Q is assumed to contain
!             the strict upper triangle of R and the lower triangle of Q
!             contains a factored form of Q, (as returned by QRFAC). 
!
!      LDQ:   Leading dimension of the Q array
!
!      RDIAG: Diagonal elements of R (also returned by QRFAC)
!
!      QTB:   On input, the first N elements contain the N-vector
!             Q(transpose)*b product. 
!
!      X:     On output, an N-vector containing solution of A*x=b
!
!      RSD:   On input, the M-vector b for which the solution/residuals
!             are desired
!
!             On output, an M-vector containing the residuals of the
!             least-squares problem, b-A*x.
!
!----------------------------------------------------------------------
      subroutine rsolve( m,n,q,ldq,qraux,qtb,x,rcalc,rsd ) 
      implicit none
      integer m,n,ldq
      double precision q(ldq,n),qraux(n),qtb(n),x(m),rsd(m)
      logical rcalc
!
      integer i,j,ju,jj
      double precision qtemp,sum,temp
!
      double precision ZERO
      parameter(ZERO=0.0d0)
!
!######################################################################
!
      do i=1,n
         x(i)=qtb(i)
      end do
!
!  --- Solve R*x = Q(transpose)*b
!
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
!
!  --- If required, calculate residual vector b-A*x 
!
 4    if (rcalc) then
         do i=1,m
            if (i.le.n) then
               rsd(i)=ZERO
            else
               rsd(i)=qtb(i)
            end if
         end do
!
         ju=min0(n,m-1)
         do jj=1,ju
            j=ju-jj+1
            if (qraux(j).ne.ZERO) then
               qtemp=q(j,j)
               q(j,j)=qraux(j)
!
               sum=ZERO
               do i=j,m
                  sum=sum+q(i,j)*rsd(j)
               end do
               temp=-sum/q(j,j)
!
               do i=j,m
                  rsd(i)=temp*q(i,j)+rsd(i)
               end do
               q(j,j)=qtemp
            end if
         end do
      end if
!
      return
      end
