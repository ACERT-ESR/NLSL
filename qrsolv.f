c NLSPMC Version 1.0 2/5/99
      subroutine qrsolv(n,r,ldr,ipvt,diag,qtb,x,sdiag,wa)
      integer n,ldr
      integer ipvt(n)
      real*8 r(ldr,n),diag(n),qtb(n),x(n),sdiag(n),wa(n)
c     **********
c
c     subroutine qrsolv
c
c     Given an m by n matrix A, an n by n diagonal matrix D,
c     and an m-vector B, the problem is to determine an X which
c     solves the system
c
c           A*X = B ,     D*X = 0 ,
c
c     in the least squares sense.
c
c     This subroutine completes the solution of the problem
c     if it is provided with the necessary information from the
c     QR factorization, with column pivoting, of A. That is, if
c     A*P = Q*R, where P is a permutationf matrix, Q has orthogonal
c     columns, and R is an upper triangular matrix with diagonal
c     elements of nonincreasing magnitude, then qrsolv expects
c     the full upper triangle of R, the permutation matrix P,
c     and the first n components of (Q transpose)*B. The system
c     A*X = B, D*X = 0, is then equivalent to
c
c                  t       t
c           R*Z = Q *B ,  P *D*P*Z = 0 ,
c
c     where X = P*Z. If this system does not have full rank,
c     then a least squares solution is obtained. On output qrsolv
c     also provides an upper triangular matrix S such that
c
c            t   t               t
C           P *(A *A + D*D)*P = S *S .
c
c     S is computed within qrsolv and may be of separate interest.
c
c     The subroutine statement is
c
c       subroutine qrsolv(n,r,ldr,ipvt,diag,qtb,x,sdiag,wa)
c
c     where
c
c       n is a positive integer input variable set to the order of r.
c
c       r is an n by n array. On input the full upper triangle
c         must contain the full upper triangle of the matrix r.
c         On output the full upper triangle is unaltered, and the
c         strict lower triangle contains the strict upper triangle
c         (transposed) of the upper triangular matrix S
c
c       ldr is a positive integer input variable not less than n
c         which specifies the leading dimension of the array r.
c
c       ipvt is an integer input array of length n which defines the
c         permutation matrix p such that A*P = Q*R. Column j of P
c         is column ipvt(j) of the identity matrix.
c
c       diag is an input array of length n which must contain the
c         diagonal elements of the matrix D.
c
c       qtb is an input array of length n which must contain the first
c         n elements of the vector (Q transpose)*B.
c
c       x is an output array of length n which contains the least
c         squares solution of the system A*X = B, D*X = 0.
c
c       sdiag is an output array of length n which contains the
c         diagonal elements of the upper triangular matrix S.
c
c       wa is a work array of length n.
c
c     Subprograms called
c
c       Fortran-supplied ... dabs,dsqrt
c
c     Argonne National Laboratory. MINPACK project. March 1980.
c     Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More
c
c     **********
      integer i,j,jp1,k,kp1,l,nsing
      real*8 cos,cotan,p5,p25,qtbpj,sin,sum,tan,temp,zero
      data p5,p25,zero /5.0d-1,2.5d-1,0.0d0/
c
c----------------------------------------------------------------------
c     Copy R and (Q transpose)*B to preserve input and initialize S.
c     In particular, save the diagonal elements of R in X.
c----------------------------------------------------------------------
      do 20 j = 1, n
         do 10 i = j, n
            r(i,j) = r(j,i)
   10       continue
         x(j) = r(j,j)
         wa(j) = qtb(j)
   20    continue
c
c----------------------------------------------------------------------
c     Eliminate the diagonal matrix D using a Givens rotation.
c----------------------------------------------------------------------
      do 100 j = 1, n
c
c        Prepare the row of D to be eliminated, locating the
c        diagonal element using P from the QR factorization.
c
         l = ipvt(j)
         if (diag(l) .ne. zero) then
            do 30 k = j, n
               sdiag(k) = zero
 30         continue
            sdiag(j) = diag(l)
c
c        The transformations to eliminate the row of D
c        modify only a single element of (Q transpose)*B
c        beyond the first n, which is initially zero.
c
            qtbpj = zero
            do 80 k = j, n
c
c           Determine a Givens rotation which eliminates the
c           appropriate element in the current row of D
c
            if (sdiag(k) .ne. zero) then
               if (dabs(r(k,k)) .ge. dabs(sdiag(k))) then
                  tan = sdiag(k)/r(k,k)
                  cos = p5/dsqrt(p25+p25*tan**2)
                  sin = cos*tan
               else
                  cotan = r(k,k)/sdiag(k)
                  sin = p5/dsqrt(p25+p25*cotan**2)
                  cos = sin*cotan
               end if
c
c           Compute the modified diagonal element of R and
c           the modified element of ((Q transpose)*B,0).
c
               r(k,k) = cos*r(k,k) + sin*sdiag(k)
               temp = cos*wa(k) + sin*qtbpj
               qtbpj = -sin*wa(k) + cos*qtbpj
               wa(k) = temp
c
c           Accumulate the tranformation in the row of S.
c
               kp1 = k + 1
               if (kp1 .le. n) then
                  do 60 i = kp1, n
                     temp = cos*r(i,k) + sin*sdiag(i)
                     sdiag(i) = -sin*r(i,k) + cos*sdiag(i)
                     r(i,k) = temp
 60               continue
               end if
            end if
 80      continue
      end if
c
c        Store the diagonal element of S and restore
c        the corresponding diagonal element of r.
c
         sdiag(j) = r(j,j)
         r(j,j) = x(j)
  100    continue
c
c     Solve the triangular system for z. If the system is
c     singular, then obtain a least squares solution.
c
      nsing = n
      do 110 j = 1, n
         if (sdiag(j) .eq. zero .and. nsing .eq. n) nsing = j - 1
         if (nsing .lt. n) wa(j) = zero
  110    continue
c
      if (nsing .ge. 1) then
         do 140 k = 1, nsing
            j = nsing - k + 1
            sum = zero
            jp1 = j + 1
            if (jp1 .le. nsing) then
               do 120 i = jp1, nsing
                  sum = sum + r(i,j)*wa(i)
 120           continue
            end if
            wa(j) = (wa(j) - sum)/sdiag(j)
 140     continue
      end if
c
c     Permute the components of Z back to components of X.
c
      do 160 j = 1, n
         l = ipvt(j)
         x(l) = wa(j)
 160  continue
c
      return
      end
