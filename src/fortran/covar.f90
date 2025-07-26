        subroutine covar(n,r,ldr,ipvt,tol,wa)
        integer n,ldr
        integer ipvt(n)
        double precision tol
        double precision r(ldr,n),wa(n)
!       **********
!
!       subroutine covar
!
!       Given an m by n matrix A, the problem is to determine
!       the covariance matrix corresponding to A, defined as
!
!                                T
!                       inverse(A *A)
!
!       This subroutine completes the solution of the problem
!       if it is provided with the necessary information from the
!       QR factorization, with column pivoting,  of A.  That is, if
!       A*P = Q*R, where P is a permutation matrix, Q has orthogonal
!       columns, and R is an upper triangular matrix with diagonal
!       elements of nonincreasing magnitude, then covar expects the
!       full upper triangle of R and the permutation matrix P.
!       the covariance matrix is then computed as
!
!                                  T     T
!                       P*inverse(R *R)*P
!
!       if A is nearly rank deficient, it may be desirable to compute
!       the covariance matrix corresponding to the linearly independent
!       columns of A.  To define the numerical rank of A, covar uses
!       the tolerance tol.  If l is the largest integer such that
!
!                       abs(R(l,l)) .gt. tol*abs(R(1,1)) ,
!
!       then covar computes the covariance matrix corresponding to
!       the first l columns of r.  For k greater than l, column
!       and row ipvt(k) of the covariance matrix are set to zero.
!
!       The subroutine statement is
!
!               subroutine covar(n,r,ldr,ipvt,tol,wa)
!
!
!       where
!
!               n is a positive integer input variable set to the order of r.
!
!               R is an n by n array.  On input, the full upper triangle must
!                 contain the full upper triangle of the matrix R.  On 
!                 output r contains the square symmetric covariance 
!                 matrix.
!
!               ldr is a positive integer input variable not less than n
!                 which specifies the leading dimension of the array R.
!
!               ipvt is an integer input array of length n which defines the
!                 permutation matrix P such that A*P = Q*R. Column j of P
!                 is column ipvt(j) of the identity matrix.
!
!               tol is a nonnegative input variable used to define the
!                 numerical rank of A in the manner described above.
!
!               wa is a work array of length n
!
!       Subprograms called
!
!               Fortran-supplied ... dabs
!
!     Argonne National Laboratory. MINPACK project. March 1980.
!     Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More
!
!       ***********
        integer i,ii,j,jj,k,km1,l
        logical sing
        double precision one,temp,tolr,zero
        data one,zero /1.0d0,0.0d0/
!
!----------------------------------------------------------------------
!       Form the inverse of R in the full upper triangle of R.
!----------------------------------------------------------------------
        tolr = tol*dabs(r(1,1))
        l = 0
        do 40 k = 1, n
           if (dabs(r(k,k)) .le. tolr) goto 50
           r(k,k) = one/r(k,k)
           km1 = k - 1
           if (km1 .ge. 1) then
              do 20 j = 1, km1
                 temp = r(k,k)*r(j,k)
                 r(j,k) = zero
                 do 10 i = 1, j
                    r(i,k) = r(i,k) - temp*r(i,j)
 10              continue
 20           continue
           end if
           l = k
 40     continue
 50     continue
!
!----------------------------------------------------------------------
!       Form the full upper triangle of the inverse of (R transpose)*R
!       in the full upper triangle of R
!----------------------------------------------------------------------
        if (l .ge. 1) then
           do 100 k = 1, l
              km1 = k - 1
              if (km1 .ge. 1) then
                 do 70 j = 1, km1
                    temp = r(j,k)
                    do 60 i = 1, j
                       r(i,j) = r(i,j) + temp*r(i,k)
 60                 continue
 70              continue
              end if
              temp = r(k,k)
              do 90 i = 1, k
                 r(i,k) = temp*r(i,k)
 90           continue
 100       continue
!
        end if
!
!----------------------------------------------------------------------
!       Form the full lower triangle of the covariance matrix
!       in the strict lower triangle of R and in wa
!----------------------------------------------------------------------
        do 130 j = 1, n
           jj = ipvt(j)
           sing = j .gt. l
           do 120 i = 1, j
              if (sing) r(i,j) = zero
              ii = ipvt(i)
              if (ii .gt. jj) r(ii,jj) = r(i,j)
              if (ii .lt. jj) r(jj,ii) = r(i,j)
 120       continue
           wa(jj) = r(j,j)
 130    continue
!
!----------------------------------------------------------------------
!       Symmetrize the covariance matrix in R.
!----------------------------------------------------------------------
        do 150 j = 1, n
           do 140 i = 1, j
              r(i,j) = r(j,i)
 140       continue
           r(j,j) = wa(j)
 150    continue
!
        return
        end
