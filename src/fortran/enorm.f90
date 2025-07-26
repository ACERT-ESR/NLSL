      function enorm(n,x)
      implicit none
      integer n
      double precision x(n),enorm
!     **********
!
!     function enorm
!
!     Given an n-vector x, this function calculates the
!     Euclidean norm of x.
!
!     The Euclidean norm is computed by accumulating the sum of
!     squares in three different sums. The sums of squares for the
!     small and large components are scaled so that no overflows
!     occur. Non-destructive underflows are permitted. Underflows
!     and overflows do not occur in the computation of the unscaled
!     sum of squares for the intermediate components.
!     The definitions of small, intermediate and large components
!     depend on two constants, rdwarf and rgiant. The main
!     restrictions on these constants are that rdwarf**2 not
!     underflow and rgiant**2 not overflow. The constants
!     given here are suitable for every known computer.
!
!     the function statement is
!
!       double precision function enorm(n,x)
!
!     where
!
!       n is a positive integer input variable.
!
!       x is an input array of length n.
!
!     Subprograms called
!
!       Fortran-supplied ... dabs,dsqrt
!
!     Argonne National Laboratory. MINPACK project. March 1980.
!     Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More
!
!     **********
      integer i
      double precision agiant,floatn,one,rdwarf,rgiant,s1,s2,s3,xabs,
     *                 x1max,x3max,zero
      data one,zero,rdwarf,rgiant /1.0d0,0.0d0,3.834d-20,1.304d19/
      s1 = zero
      s2 = zero
      s3 = zero
      x1max = zero
      x3max = zero
      enorm = 0.0d0
      floatn = n
      agiant = rgiant/floatn
      do 90 i = 1, n
         xabs = dabs(x(i))

         if (xabs .ge. agiant) then
!
!     ----- Sum for large components --------
!
            if (xabs .gt. x1max) then
               s1 = one + s1*(x1max/xabs)**2
               x1max = xabs
            else
               s1 = s1 + (xabs/x1max)**2
            end if

            else if (xabs .le. rdwarf) then
!
!     -----  Sum for small components -------
!
               if (xabs .gt. x3max) then
                  s3 = one + s3*(x3max/xabs)**2
                  x3max = xabs
               else
                  if (xabs .ne. zero) s3 = s3 + (xabs/x3max)**2
               end if
!
            else
!
!     ----- Sum for intermediate components ------
!
            s2 = s2 + xabs**2
         end if
 90   continue
!
!     ----- Calculation of norm ------
!
      if (s1 .ne. zero) then
         enorm = x1max*dsqrt(s1+(s2/x1max)/x1max)
      else
         if (s2 .ne. zero) then
            if (s2 .ge. x3max)
     *         enorm = dsqrt(s2*(one+(x3max/s2)*(x3max*s3)))
            if (s2 .lt. x3max)
     *         enorm = dsqrt(x3max*((s2/x3max)+(x3max*s3)))
         else
            enorm = x3max*dsqrt(s3)
         end if
      end if
      return
      end
