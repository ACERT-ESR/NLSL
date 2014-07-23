c NLSPMC Version 1.0 2/5/99
      function dpmpar(i)
      implicit none
      double precision dpmpar
      integer i
c     **********
c
c     function dpmpar
c
c     this function provides double precision machine parameters
c
c     the function statement is
c
c       double precision function dpmpar(i)
c
c     where
c
c       i is an integer input variable set to 1, 2, or 3 which
c         selects the desired machine parameter. if the machine has
c         t base b digits and its smallest and largest exponents are
c         emin and emax, respectively, then these parameters are
c
c         dpmpar(1) = b**(1 - t), the machine precision,
c
c         dpmpar(2) = b**(emin - 1), the smallest magnitude,
c
c         dpmpar(3) = b**emax*(1 - b**(-t)), the largest magnitude.
c
c     Argonne National Laboratory. MINPACK project. March 1980.
c     Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More
c
c     *******
c     Adapted for use on IBM 3090  David Budil May 1991
c
c        t = 56, emin=-128, emax=127
c     **********
c
      real*8 dmach(3)
      integer*2 imach(12)
      equivalence(imach,dmach)
c
        data imach / 13328, 0, 0, 0,
     1               16, 0, 0, 0,
     2               32767, -1, -1, -1 /
c
      dpmpar = dmach(i)
      return
      end
