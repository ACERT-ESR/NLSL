      subroutine qrfac(m,n,a,lda,pivot,ipvt,lipvt,rdiag,acnorm,wa)
      use rnddbl
      integer m,n,lda,lipvt
      integer ipvt(n)
      logical pivot
      double precision a(lda,n),rdiag(n),acnorm(n),wa(n)
!     **********
!
!     subroutine qrfac
!
!     This subroutine uses Householder transformations with column
!     pivoting (optional) to compute a QR factorization of the
!     m by n matrix A. That is, qrfac determines an orthogonal
!     matrix Q, a permutation matrix P, and an upper trapezoidal
!     matrix R with diagonal elements of nonincreasing magnitude,
!     such that A*P = Q*R. The Householder transformation for
!     column k, k = 1,2,...,min(m,n), is of the form
!
!                           T
!           i - (1/u(k))*u*u
!
!     where u has zeros in the first k-1 positions. The form of
!     this transformation and the method of pivoting first
!     appeared in the corresponding LINPACK subroutine.
!
!     The subroutine statement is
!
!       subroutine qrfac(m,n,a,lda,pivot,ipvt,lipvt,rdiag,acnorm,wa)
!
!     where
!
!       m is a positive integer input variable set to the number
!         of rows of a.
!
!       n is a positive integer input variable set to the number
!         of columns of a.
!
!       a is an m by n array. On input a contains the matrix for
!         which the qr factorization is to be computed. On output
!         the strict upper trapezoidal part of A contains the strict
!         upper trapezoidal part of R, and the lower trapezoidal
!         part of a contains a factored form of Q (the non-trivial
!         elements of the u vectors described above).
!
!       lda is a positive integer input variable not less than m
!         which specifies the leading dimension of the array a.
!
!       pivot is a logical input variable. If pivot is set true,
!         then column pivoting is enforced. If pivot is set false,
!         then no column pivoting is done.
!
!       ipvt is an integer output array of length n. ipvt
!         defines the permutation matrix P such that A*P = Q*R.
!         Column j of P is column ipvt(j) of the identity matrix.
!         If pivot is false, ipvt is not referenced.
!
!       lipvt is a positive integer input variable. If pivot is true
!         then lipvt specifies the number of leading columns in the 
!         Jacobian that may be pivoted. Columns lipvt+1 through n are
!         fixed at the right hand side of the matrix.
!
!       rdiag is an output array of length n which contains the
!         diagonal elements of r.
!
!       acnorm is an output array of length n which contains the
!         norms of the corresponding columns of the input matrix a.
!         if this information is not needed, then acnorm can coincide
!         with rdiag.
!
!       wa is a work array of length n. if pivot is false, then wa
!         can coincide with rdiag.
!
!     Subprograms called
!
!       MINPACK-supplied ... enorm
!
!       Fortran-supplied ... dmax1,dsqrt,min0
!
!     Argonne National Laboratory. MINPACK project. March 1980.
!     Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More
!
!     **********
      integer i,j,jp1,k,kmax,minmn
      double precision ajnorm,epsmch,one,p05,sum,temp,zero
      PARAMETER (EPSMCH=RNDOFF)
      double precision enorm
      data one,p05,zero /1.0d0,5.0d-2,0.0d0/
!
!     epsmch is the machine precision.
!
      if (lipvt .gt. n) lipvt = n
!
!----------------------------------------------------------------------
!     Compute the initial column norms and initialize several arrays.
!----------------------------------------------------------------------
      do j = 1, n
         acnorm(j) = enorm(m,a(1,j))
         rdiag(j) = acnorm(j)
         wa(j) = rdiag(j)
         if (pivot) ipvt(j) = j
      end do
!
!----------------------------------------------------------------------
!     Reduce A to R with Householder transformations.
!----------------------------------------------------------------------
      minmn = min0(m,n)
      do j = 1, minmn
         if (pivot) then
!
!----------------------------------------------------------------------
!        Bring the column of largest norm into the pivot position.
!----------------------------------------------------------------------
          kmax = j
          do k = j, lipvt
             if (rdiag(k) .gt. rdiag(kmax)) kmax = k
          end do
!
          if (j.ne.kmax) then
             do i = 1, m
                temp = a(i,j)
                a(i,j) = a(i,kmax)
                a(i,kmax) = temp
             end do
!
             rdiag(kmax) = rdiag(j)
             wa(kmax) = wa(j)
             k = ipvt(j)
             ipvt(j) = ipvt(kmax)
             ipvt(kmax) = k
          end if
       end if
!
!----------------------------------------------------------------------
!        Compute the Householder transformation to reduce the
!        j-th column of a to a multiple of the j-th unit vector.
!----------------------------------------------------------------------
       ajnorm = enorm(m-j+1,a(j,j))
!
       if (ajnorm .ne. zero) then
          if (a(j,j) .lt. zero) ajnorm = -ajnorm
          do i = j, m
             a(i,j) = a(i,j)/ajnorm
          end do
!
          a(j,j) = a(j,j) + one
!
!----------------------------------------------------------------------
!        Apply the transformation to the remaining columns
!        and update the norms.
!----------------------------------------------------------------------
          jp1 = j + 1
          if (jp1.le.n) then
             do k = jp1, n
                sum = zero
                do i = j, m
                   sum = sum + a(i,j)*a(i,k)
                end do
!
                temp = sum/a(j,j)
               do i = j, m
                  a(i,k) = a(i,k) - temp*a(i,j)
               end do
!
               if (pivot .and. rdiag(k) .ne. zero) then
                  temp = a(j,k)/rdiag(k)
                  rdiag(k) = rdiag(k)*dsqrt(dmax1(zero,one-temp**2))
                  if (p05*(rdiag(k)/wa(k))**2 .le. epsmch) then
                     rdiag(k) = enorm(m-j,a(jp1,k))
                     wa(k) = rdiag(k)
                  end if
               end if
            end do
!
         end if
        end if
!
        rdiag(j) = -ajnorm
      end do
!
      return
      end
