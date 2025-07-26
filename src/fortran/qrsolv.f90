      subroutine qrsolv(n,r,ldr,ipvt,diag,qtb,x,sdiag,wa)
      integer n,ldr
      integer ipvt(n)
      double precision r(ldr,n),diag(n),qtb(n),x(n),sdiag(n),wa(n)
!     **********
!
!     subroutine qrsolv
!
!     Given an m by n matrix A, an n by n diagonal matrix D,
!     and an m-vector B, the problem is to determine an X which
!     solves the system
!
!           A*X = B ,     D*X = 0 ,
!
!     in the least squares sense.
!
!     This subroutine completes the solution of the problem
!     if it is provided with the necessary information from the
!     QR factorization, with column pivoting, of A. That is, if
!     A*P = Q*R, where P is a permutationf matrix, Q has orthogonal
!     columns, and R is an upper triangular matrix with diagonal
!     elements of nonincreasing magnitude, then qrsolv expects
!     the full upper triangle of R, the permutation matrix P,
!     and the first n components of (Q transpose)*B. The system
!     A*X = B, D*X = 0, is then equivalent to
!
!                  t       t
!           R*Z = Q *B ,  P *D*P*Z = 0 ,
!
!     where X = P*Z. If this system does not have full rank,
!     then a least squares solution is obtained. On output qrsolv
!     also provides an upper triangular matrix S such that
!
!            t   t               t
C           P *(A *A + D*D)*P = S *S .
!
!     S is computed within qrsolv and may be of separate interest.
!
!     The subroutine statement is
!
!       subroutine qrsolv(n,r,ldr,ipvt,diag,qtb,x,sdiag,wa)
!
!     where
!
!       n is a positive integer input variable set to the order of r.
!
!       r is an n by n array. On input the full upper triangle
!         must contain the full upper triangle of the matrix r.
!         On output the full upper triangle is unaltered, and the
!         strict lower triangle contains the strict upper triangle
!         (transposed) of the upper triangular matrix S
!
!       ldr is a positive integer input variable not less than n
!         which specifies the leading dimension of the array r.
!
!       ipvt is an integer input array of length n which defines the
!         permutation matrix p such that A*P = Q*R. Column j of P
!         is column ipvt(j) of the identity matrix.
!
!       diag is an input array of length n which must contain the
!         diagonal elements of the matrix D.
!
!       qtb is an input array of length n which must contain the first
!         n elements of the vector (Q transpose)*B.
!
!       x is an output array of length n which contains the least
!         squares solution of the system A*X = B, D*X = 0.
!
!       sdiag is an output array of length n which contains the
!         diagonal elements of the upper triangular matrix S.
!
!       wa is a work array of length n.
!
!     Subprograms called
!
!       Fortran-supplied ... dabs,dsqrt
!
!     Argonne National Laboratory. MINPACK project. March 1980.
!     Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More
!
!     **********
      integer i,j,jp1,k,kp1,l,nsing
      double precision cos,cotan,p5,p25,qtbpj,sin,sum,tan,temp,zero
      data p5,p25,zero /5.0d-1,2.5d-1,0.0d0/
!
!----------------------------------------------------------------------
!     Copy R and (Q transpose)*B to preserve input and initialize S.
!     In particular, save the diagonal elements of R in X.
!----------------------------------------------------------------------
      do 20 j = 1, n
         do 10 i = j, n
            r(i,j) = r(j,i)
   10       continue
         x(j) = r(j,j)
         wa(j) = qtb(j)
   20    continue
!
!----------------------------------------------------------------------
!     Eliminate the diagonal matrix D using a Givens rotation.
!----------------------------------------------------------------------
      do 100 j = 1, n
!
!        Prepare the row of D to be eliminated, locating the
!        diagonal element using P from the QR factorization.
!
         l = ipvt(j)
         if (diag(l) .ne. zero) then
            do 30 k = j, n
               sdiag(k) = zero
 30         continue
            sdiag(j) = diag(l)
!
!        The transformations to eliminate the row of D
!        modify only a single element of (Q transpose)*B
!        beyond the first n, which is initially zero.
!
            qtbpj = zero
            do 80 k = j, n
!
!           Determine a Givens rotation which eliminates the
!           appropriate element in the current row of D
!
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
!
!           Compute the modified diagonal element of R and
!           the modified element of ((Q transpose)*B,0).
!
               r(k,k) = cos*r(k,k) + sin*sdiag(k)
               temp = cos*wa(k) + sin*qtbpj
               qtbpj = -sin*wa(k) + cos*qtbpj
               wa(k) = temp
!
!           Accumulate the tranformation in the row of S.
!
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
!
!        Store the diagonal element of S and restore
!        the corresponding diagonal element of r.
!
         sdiag(j) = r(j,j)
         r(j,j) = x(j)
  100    continue
!
!     Solve the triangular system for z. If the system is
!     singular, then obtain a least squares solution.
!
      nsing = n
      do 110 j = 1, n
         if (sdiag(j) .eq. zero .and. nsing .eq. n) nsing = j - 1
         if (nsing .lt. n) wa(j) = zero
  110    continue
!
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
!
!     Permute the components of Z back to components of X.
!
      do 160 j = 1, n
         l = ipvt(j)
         x(l) = wa(j)
 160  continue
!
      return
      end
