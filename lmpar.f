c NLSPMC Version 1.0 2/5/99
      subroutine lmpar(n,r,ldr,ipvt,diag,qtf,delta,par,x,sdiag,wa1,
     *                 wa2,gnvec,gradf)
      integer n,ldr
      integer ipvt(n)
      real*8 delta,par
      real*8 r(ldr,n),diag(n),qtf(n),x(n),sdiag(n),wa1(n),
     *                 wa2(n),gnvec(n), gradf(n)
c     **********
c
c     Subroutine lmpar
c
c     Given an m by n matrix a, an n by n nonsingular diagonal
c     matrix d, an m-vector b, and a positive number delta,
c     the problem is to determine a value for the parameter
c     par such that if x solves the system
c
c         A*x = f ,     sqrt(par)*D*x = 0 ,
c
c     in the least squares sense, and dxnorm is the Euclidean
c     norm of D*x, then either par is zero and
c
c           (dxnorm-delta) .le. 0.1*delta ,
c
c     or par is positive and
c
c           abs(dxnorm-delta) .le. 0.1*delta .
c
c     This subroutine completes the solution of the problem
c     if it is provided with the necessary information from the
c     QR factorization, with column pivoting, of A. That is, if
c     A*P = Q*R, where P is a permutation matrix, Q has orthogonal
c     columns, and R is an upper triangular matrix with diagonal
c     elements of nonincreasing magnitude, then lmpar expects
c     the full upper triangle of R, the permutation matrix P,
c     and the first n components of (Q transpose)*b. On output
c     lmpar also provides an upper triangular matrix S such that
c
c            t   t                   t
c           P *(A *A + par*D*D)*P = S *S .
c
c     S is employed within lmpar and may be of separate interest.
c
c     Only a few iterations are generally needed for convergence
c     of the algorithm. If, however, the limit of 10 iterations
c     is reached, then the output par will contain the best
c     value obtained so far.
c
c     The subroutine statement is
c
c       subroutine lmpar(n,r,ldr,ipvt,diag,qtf,delta,par,x,sdiag,
c                        wa1,wa2)
c
c     where
c
c       n is a positive integer input variable set to the order of R.
c
c       R is an n by n array. On input the full upper triangle
c         must contain the full upper triangle of the matrix R.
c         On output the full upper triangle is unaltered, and the
c         strict lower triangle contains the strict upper triangle
c         (transposed) of the upper triangular matrix S.
c
c       ldr is a positive integer input variable not less than n
c         which specifies the leading dimension of the array r.
c
c       ipvt is an integer input array of length n which defines the
c         permutation matrix P such that A*P = Q*R. Column j of P
c         is column ipvt(j) of the identity matrix.
c
c       diag is an input array of length n which must contain the
c         diagonal elements of the matrix D.
c
c       qtf is an input array of length n which must contain the first
c         n elements of the vector (Q transpose)*f (f is the vector
c         containing the function values to be minimized).
c
c       delta is a positive input variable which specifies an upper
c         bound on the Euclidean norm of D*x.
c
c       par is a nonnegative variable. On input par contains an
c         initial estimate of the Levenberg-Marquardt parameter.
c         On output par contains the final estimate.
c
c       x is an output array of length n which contains the least
c         squares solution of the system A*x = b, sqrt(par)*D*x = 0,
c         for the output par.
c
c       sdiag is an output array of length n which contains the
c         diagonal elements of the upper triangular matrix S.
c
c       gnvec is an output array of length n which contains the 
c         Gauss-Newton vector corresponding to the input matrix R.
c
c       gradf is an output array of length n which contains the
c         gradient (steepest descent) vector
c
c       wa1 and wa2 are work arrays of length n.
c
c     Subprograms called
c
c       MINPACK-supplied ... dpmpar,enorm,qrsolv
c
c       Fortran-supplied ... dabs,dmax1,dmin1,dsqrt
c
c     Argonne National Laboratory. MINPACK project. March 1980.
c     Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More
c
c     **********
      integer i,iter,j,jm1,jp1,k,l,nsing
      real*8 dxnorm,dwarf,fp,gnorm,parc,parl,paru,p1,p001,
     *                 sum,temp,zero
      real*8 dpmpar,enorm
      data p1,p001,zero /1.0d-1,1.0d-3,0.0d0/
c
c     dwarf is the smallest positive magnitude.
c
      dwarf = dpmpar(2)
c
c----------------------------------------------------------------------
c     Compute and store in x the Gauss-Newton direction. If the 
c     Jacobian is rank-deficient, obtain a least squares solution.
c     Store the scaled Gauss-Newton direction in gnvec.
c
c     The Gauss-Newton direction is the solution of the problem given
c     above for par=0. The solution to the above least-squares problem
c     in this case is equivalent to finding the vector
c
c           +        -1 T         +
c     x = -J f  = - R  Q f (i.e. J  is the pseudoinverse of the Jacobian)
c
c     See exercise 24 in chapter 6 of Dennis and Schnabel.
c----------------------------------------------------------------------
      nsing = n
      do 10 j = 1, n
         wa1(j) = qtf(j)
         if (r(j,j) .eq. zero .and. nsing .eq. n) nsing = j - 1
         if (nsing .lt. n) wa1(j) = zero
   10    continue
c
      if (nsing .ge. 1) then
         do 40 k = 1, nsing
            j = nsing - k + 1
            wa1(j) = wa1(j)/r(j,j)
            temp = wa1(j)
            jm1 = j - 1
         if (jm1 .ge. 1) then
            do 20 i = 1, jm1
               wa1(i) = wa1(i) - r(i,j)*temp
 20         continue
         end if
 40   continue
      end if
c
      do 60 j = 1, n
         l = ipvt(j)
         x(l) = wa1(j)
         gnvec(l)=x(l)
   60    continue
c
c----------------------------------------------------------------------
c     Initialize the iteration counter. If the (scaled) Gauss-Newton 
c     step size is within 10% of the trust-region bound, terminate the 
c     search (par will be returned as zero).
c----------------------------------------------------------------------
      iter = 0
      do 70 j = 1, n
         wa2(j) = diag(j)*x(j)
   70    continue
      dxnorm = enorm(n,wa2)
      fp = dxnorm - delta
      if (fp .le. p1*delta) go to 220
c
c----------------------------------------------------------------------
c     If the Jacobian is not rank deficient, the Gauss-Newton
c     step provides a lower bound, parl, for the zero of
c     the function. Otherwise set this bound to zero.
c----------------------------------------------------------------------
      parl = zero
      if (nsing .ge. n) then
         do 80 j = 1, n
            l = ipvt(j)
            wa1(j) = diag(l)*(wa2(l)/dxnorm)
 80      continue
         do 110 j = 1, n
            sum = zero
            jm1 = j - 1
            if (jm1 .ge. 1) then
               do 90 i = 1, jm1
                  sum = sum + r(i,j)*wa1(i)
 90            continue
            end if
            wa1(j) = (wa1(j) - sum)/r(j,j)
 110     continue
         temp = enorm(n,wa1)
         parl = ((fp/delta)/temp)/temp
      end if
c
c----------------------------------------------------------------------
c     Calculate the steepest descent direction of the residuals at
c     point where fvec has been evaluated. 
c
c                T       T T
c     Grad(f) = J f = P R Q f
c
c     This provides an upper bound, paru, for the zero of the function.
c----------------------------------------------------------------------
      do 140 j = 1, n
         sum = zero
         do 130 i = 1, j
            sum = sum + r(i,j)*qtf(i)
  130       continue
         l = ipvt(j)
         wa1(j) = sum/diag(l)
         gradf(j) = wa1(j)
  140    continue
      gnorm = enorm(n,wa1)
      paru = gnorm/delta
      if (paru .eq. zero) paru = dwarf/dmin1(delta,p1)
c
c----------------------------------------------------------------------
c     If the input par lies outside of the interval (parl,paru),
c     set par to the closer endpoint.
c----------------------------------------------------------------------
      par = dmax1(par,parl)
      par = dmin1(par,paru)
      if (par .eq. zero) par = gnorm/dxnorm
c
c----------------------------------------------------------------------
c     *** Beginning of an iteration.
c----------------------------------------------------------------------
  150 continue
         iter = iter + 1
c
c----------------------------------------------------------------------
c        Evaluate the function at the current value of par.
c----------------------------------------------------------------------
         if (par .eq. zero) par = dmax1(dwarf,p001*paru)
         temp = dsqrt(par)
         do 160 j = 1, n
            wa1(j) = temp*diag(j)
  160       continue
         call qrsolv(n,r,ldr,ipvt,wa1,qtf,x,sdiag,wa2)
         do 170 j = 1, n
            wa2(j) = diag(j)*x(j)
  170       continue
         dxnorm = enorm(n,wa2)
         temp = fp
         fp = dxnorm - delta
c
c----------------------------------------------------------------------
c        If the function is small enough, accept the current value
c        of par. Also test for the exceptional cases where parl
c        is zero or the number of iterations has reached 10.
c----------------------------------------------------------------------
         if (dabs(fp) .le. p1*delta
     *       .or. parl .eq. zero .and. fp .le. temp
     *            .and. temp .lt. zero .or. iter .eq. 10) go to 220
c
c----------------------------------------------------------------------
c        Compute the Newton correction.
c----------------------------------------------------------------------
         do 180 j = 1, n
            l = ipvt(j)
            wa1(j) = diag(l)*(wa2(l)/dxnorm)
  180       continue
         do 210 j = 1, n
            wa1(j) = wa1(j)/sdiag(j)
            temp = wa1(j)
            jp1 = j + 1
            if (n .ge. jp1) then
               do 190 i = jp1, n
                  wa1(i) = wa1(i) - r(i,j)*temp
 190           continue
            end if
  210       continue
         temp = enorm(n,wa1)
         parc = ((fp/delta)/temp)/temp
c
c----------------------------------------------------------------------
c        Depending on the sign of the function, update parl or paru.
c----------------------------------------------------------------------
         if (fp .gt. zero) parl = dmax1(parl,par)
         if (fp .lt. zero) paru = dmin1(paru,par)
c
c----------------------------------------------------------------------
c        Compute an improved estimate for par.
c----------------------------------------------------------------------
         par = dmax1(parl,par+parc)
c
c----------------------------------------------------------------------
c        *** End of an iteration.
c----------------------------------------------------------------------
         go to 150
  220 continue
c
c----------------------------------------------------------------------
c     Termination.
c----------------------------------------------------------------------
      if (iter .eq. 0) par = zero
c
      return
      end
