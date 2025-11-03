      subroutine dchex(r,ldr,p,k,l,z,ldz,nz,c,s,job)
      integer ldr,p,k,l,ldz,nz,job
      double precision r(ldr,*),z(ldz,*),s(*)
      double precision c(*)
!
!     dchex updates the cholesky factorization
!
!                   a = trans(r)*r
!
!     of a positive definite matrix a of order p under diagonal
!     permutations of the form
!
!                   trans(e)*a*e
!
!     where e is a permutation matrix.  specifically, given
!     an upper triangular matrix r and a permutation matrix
!     e (which is specified by k, l, and job), dchex determines
!     a orthogonal matrix u such that
!
!                           u*r*e = rr,
!
!     where rr is upper triangular.  at the users option, the
!     transformation u will be multiplied into the array z.
!     if a = trans(x)*x, so that r is the triangular part of the
!     qr factorization of x, then rr is the triangular part of the
!     qr factorization of x*e, i.e. x with its columns permuted.
!     for a less terse description of what dchex does and how
!     it may be applied, see the linpack guide.
!
!     the matrix q is determined as the product u(l-k)*...*u(1)
!     of plane rotations of the form
!
!                           (    c(i)       s(i) )
!                           (                    ) ,
!                           (    -s(i)      c(i) )
!
!     where c(i) is double precision, the rows these rotations operate
!     on are described below.
!
!     there are two types of permutations, which are determined
!     by the value of job.
!
!     1. right circular shift (job = 1).
!
!         the columns are rearranged in the following order.
!
!                1,...,k-1,l,k,k+1,...,l-1,l+1,...,p.
!
!         u is the product of l-k rotations u(i), where u(i)
!         acts in the (l-i,l-i+1)-plane.
!
!     2. left circular shift (job = 2).
!         the columns are rearranged in the following order
!
!                1,...,k-1,k+1,k+2,...,l,k,l+1,...,p.
!
!         u is the product of l-k rotations u(i), where u(i)
!         acts in the (k+i-1,k+i)-plane.
!
!     on entry
!
!         r      double precision(ldr,p), where ldr.ge.p.
!                r contains the upper triangular factor
!                that is to be updated.  elements of r
!                below the diagonal are not referenced.
!
!         ldr    integer.
!                ldr is the leading dimension of the array r.
!
!         p      integer.
!                p is the order of the matrix r.
!
!         k      integer.
!                k is the first column to be permuted.
!
!         l      integer.
!                l is the last column to be permuted.
!                l must be strictly greater than k.
!
!         z      double precision(ldz,nz), where ldz.ge.p.
!                z is an array of nz p-vectors into which the
!                transformation u is multiplied.  z is
!                not referenced if nz = 0.
!
!         ldz    integer.
!                ldz is the leading dimension of the array z.
!
!         nz     integer.
!                nz is the number of columns of the matrix z.
!
!         job    integer.
!                job determines the type of permutation.
!                       job = 1  right circular shift.
!                       job = 2  left circular shift.
!
!     on return
!
!         r      contains the updated factor.
!
!         z      contains the updated matrix z.
!
!         c      double precision(p).
!                c contains the cosines of the transforming rotations.
!
!         s      double precision(p).
!                s contains the sines of the transforming rotations.
!
!     linpack. this version dated 08/14/78 .
!     g.w. stewart, university of maryland, argonne national lab.
!
!     dchex uses the following functions and subroutines.
!
!     blas drotg
!     fortran min0
!
      integer i,ii,il,iu,j,jj,km1,kp1,lmk,lm1
      double precision t
!
!     initialize
!
      km1 = k - 1
      kp1 = k + 1
      lmk = l - k
      lm1 = l - 1
!
!     perform the appropriate task.
!
      go to (10,130), job
!
!     right circular shift.
!
   10 continue
!
!        reorder the columns.
!
         do 20 i = 1, l
            ii = l - i + 1
            s(i) = r(ii,l)
   20    continue
         do 40 jj = k, lm1
            j = lm1 - jj + k
            do 30 i = 1, j
               r(i,j+1) = r(i,j)
   30       continue
            r(j+1,j+1) = 0.0d0
   40    continue
         if (k .eq. 1) go to 60
            do 50 i = 1, km1
               ii = l - i + 1
               r(i,k) = s(ii)
   50       continue
   60    continue
!
!        calculate the rotations.
!
         t = s(1)
         do 70 i = 1, lmk
            call drotg(s(i+1),t,c(i),s(i))
            t = s(i+1)
   70    continue
         r(k,k) = t
         do 90 j = kp1, p
            il = max0(1,l-j+1)
            do 80 ii = il, lmk
               i = l - ii
               t = c(ii)*r(i,j) + s(ii)*r(i+1,j)
               r(i+1,j) = c(ii)*r(i+1,j) - s(ii)*r(i,j)
               r(i,j) = t
   80       continue
   90    continue
!
!        if required, apply the transformations to z.
!
         if (nz .lt. 1) go to 120
         do 110 j = 1, nz
            do 100 ii = 1, lmk
               i = l - ii
               t = c(ii)*z(i,j) + s(ii)*z(i+1,j)
               z(i+1,j) = c(ii)*z(i+1,j) - s(ii)*z(i,j)
               z(i,j) = t
  100       continue
  110    continue
  120    continue
      go to 260
!
!     left circular shift
!
  130 continue
!
!        reorder the columns
!
         do 140 i = 1, k
            ii = lmk + i
            s(ii) = r(i,k)
  140    continue
         do 160 j = k, lm1
            do 150 i = 1, j
               r(i,j) = r(i,j+1)
  150       continue
            jj = j - km1
            s(jj) = r(j+1,j+1)
  160    continue
         do 170 i = 1, k
            ii = lmk + i
            r(i,l) = s(ii)
  170    continue
         do 180 i = kp1, l
            r(i,l) = 0.0d0
  180    continue
!
!        reduction loop.
!
         do 220 j = k, p
            if (j .eq. k) go to 200
!
!              apply the rotations.
!
               iu = min0(j-1,l-1)
               do 190 i = k, iu
                  ii = i - k + 1
                  t = c(ii)*r(i,j) + s(ii)*r(i+1,j)
                  r(i+1,j) = c(ii)*r(i+1,j) - s(ii)*r(i,j)
                  r(i,j) = t
  190          continue
  200       continue
            if (j .ge. l) go to 210
               jj = j - k + 1
               t = s(jj)
               call drotg(r(j,j),t,c(jj),s(jj))
  210       continue
  220    continue
!
!        apply the rotations to z.
!
         if (nz .lt. 1) go to 250
         do 240 j = 1, nz
            do 230 i = k, lm1
               ii = i - km1
               t = c(ii)*z(i,j) + s(ii)*z(i+1,j)
               z(i+1,j) = c(ii)*z(i+1,j) - s(ii)*z(i,j)
               z(i,j) = t
  230       continue
  240    continue
  250    continue
  260 continue
      return
      end
