      subroutine qrfac(m,n,a,lda,pivot,ipvt,lipvt,rdiag,acnorm,wa)
      use rnddbl
      integer m,n,lda,lipvt
      integer ipvt(n)
      logical pivot
      double precision a(lda,n),rdiag(n),acnorm(n),wa(n)
c     **********
c
c     subroutine qrfac
c
c     This subroutine uses Householder transformations with column
c     pivoting (optional) to compute a QR factorization of the
c     m by n matrix A. That is, qrfac determines an orthogonal
c     matrix Q, a permutation matrix P, and an upper trapezoidal
c     matrix R with diagonal elements of nonincreasing magnitude,
c     such that A*P = Q*R. The Householder transformation for
c     column k, k = 1,2,...,min(m,n), is of the form
c
c                           T
c           i - (1/u(k))*u*u
c
c     where u has zeros in the first k-1 positions. The form of
c     this transformation and the method of pivoting first
c     appeared in the corresponding LINPACK subroutine.
c
c     The subroutine statement is
c
c       subroutine qrfac(m,n,a,lda,pivot,ipvt,lipvt,rdiag,acnorm,wa)
c
c     where
c
c       m is a positive integer input variable set to the number
c         of rows of a.
c
c       n is a positive integer input variable set to the number
c         of columns of a.
c
c       a is an m by n array. On input a contains the matrix for
c         which the qr factorization is to be computed. On output
c         the strict upper trapezoidal part of A contains the strict
c         upper trapezoidal part of R, and the lower trapezoidal
c         part of a contains a factored form of Q (the non-trivial
c         elements of the u vectors described above).
c
c       lda is a positive integer input variable not less than m
c         which specifies the leading dimension of the array a.
c
c       pivot is a logical input variable. If pivot is set true,
c         then column pivoting is enforced. If pivot is set false,
c         then no column pivoting is done.
c
c       ipvt is an integer output array of length n. ipvt
c         defines the permutation matrix P such that A*P = Q*R.
c         Column j of P is column ipvt(j) of the identity matrix.
c         If pivot is false, ipvt is not referenced.
c
c       lipvt is a positive integer input variable. If pivot is true
c         then lipvt specifies the number of leading columns in the 
c         Jacobian that may be pivoted. Columns lipvt+1 through n are
c         fixed at the right hand side of the matrix.
c
c       rdiag is an output array of length n which contains the
c         diagonal elements of r.
c
c       acnorm is an output array of length n which contains the
c         norms of the corresponding columns of the input matrix a.
c         if this information is not needed, then acnorm can coincide
c         with rdiag.
c
c       wa is a work array of length n. if pivot is false, then wa
c         can coincide with rdiag.
c
c     Subprograms called
c
c       MINPACK-supplied ... enorm
c
c       Fortran-supplied ... dmax1,dsqrt,min0
c
c     Argonne National Laboratory. MINPACK project. March 1980.
c     Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More
c
c     **********
      integer i,j,jp1,k,kmax,minmn
      double precision ajnorm,epsmch,one,p05,sum,temp,zero
      PARAMETER (EPSMCH=RNDOFF)
      double precision enorm
      data one,p05,zero /1.0d0,5.0d-2,0.0d0/
c
c     epsmch is the machine precision.
c
      if (lipvt .gt. n) lipvt = n
c
c----------------------------------------------------------------------
c     Compute the initial column norms and initialize several arrays.
c----------------------------------------------------------------------
      do j = 1, n
         acnorm(j) = enorm(m,a(1,j))
         rdiag(j) = acnorm(j)
         wa(j) = rdiag(j)
         if (pivot) ipvt(j) = j
      end do
c
c----------------------------------------------------------------------
c     Reduce A to R with Householder transformations.
c----------------------------------------------------------------------
      minmn = min0(m,n)
      do j = 1, minmn
         if (pivot) then
c
c----------------------------------------------------------------------
c        Bring the column of largest norm into the pivot position.
c----------------------------------------------------------------------
          kmax = j
          do k = j, lipvt
             if (rdiag(k) .gt. rdiag(kmax)) kmax = k
          end do
c
          if (j.ne.kmax) then
             do i = 1, m
                temp = a(i,j)
                a(i,j) = a(i,kmax)
                a(i,kmax) = temp
             end do
c
             rdiag(kmax) = rdiag(j)
             wa(kmax) = wa(j)
             k = ipvt(j)
             ipvt(j) = ipvt(kmax)
             ipvt(kmax) = k
          end if
       end if
c
c----------------------------------------------------------------------
c        Compute the Householder transformation to reduce the
c        j-th column of a to a multiple of the j-th unit vector.
c----------------------------------------------------------------------
       ajnorm = enorm(m-j+1,a(j,j))
c
       if (ajnorm .ne. zero) then
          if (a(j,j) .lt. zero) ajnorm = -ajnorm
          do i = j, m
             a(i,j) = a(i,j)/ajnorm
          end do
c
          a(j,j) = a(j,j) + one
c
c----------------------------------------------------------------------
c        Apply the transformation to the remaining columns
c        and update the norms.
c----------------------------------------------------------------------
          jp1 = j + 1
          if (jp1.le.n) then
             do k = jp1, n
                sum = zero
                do i = j, m
                   sum = sum + a(i,j)*a(i,k)
                end do
c
                temp = sum/a(j,j)
               do i = j, m
                  a(i,k) = a(i,k) - temp*a(i,j)
               end do
c
               if (pivot .and. rdiag(k) .ne. zero) then
                  temp = a(j,k)/rdiag(k)
                  rdiag(k) = rdiag(k)*dsqrt(dmax1(zero,one-temp**2))
                  if (p05*(rdiag(k)/wa(k))**2 .le. epsmch) then
                     rdiag(k) = enorm(m-j,a(jp1,k))
                     wa(k) = rdiag(k)
                  end if
               end if
            end do
c
         end if
        end if
c
        rdiag(j) = -ajnorm
      end do
c
      return
      end
