! Version 1.4 10/10/94
!*********************************************************************
!
!       This complex double precision function returns the dot
!       product of the two input vectors without complex 
!       conjugation.
!
!       Calling Parameters:
!       ------------------
!
!       vectors:
!               x    : input vector
!               y    : input vector
!
!       scalars:
!               ndim : dimension of the vectors x and y
!
!
!       Local Variables:
!       ---------------
!               ir   : row counter
!               accr : real part of dot product
!               acci : imaginary part of dot product
!
!       Notes:
!       -----
!               It is assumed that the input vectors are stored in a
!               double precision complex array or as a real double 
!               precision array dimensioned as 2 X mxdim in the 
!               calling routine.
!
!
!       Includes:
!               nlsdim.inc
!               rndoff.inc
!
!       Uses:
!
!       written by DJS 26-AUG-87
!
!*********************************************************************
!
      function zdotu(x,y,ndim)
!
      use nlsdim
      use rnddbl
!
      implicit none
!
      double complex zdotu
!
      integer ndim
      double precision x,y
      dimension x(2,MXDIM),y(2,MXDIM)
!
      integer ir
      double precision accx,accy,accr,acci,scale
!
      double precision ZERO
      parameter (ZERO=0.0D0)
!
!######################################################################
!
      accx=ZERO
      accy=ZERO
      accr=ZERO
      acci=ZERO
!
      do ir=1,ndim
         accx=accx+x(1,ir)*x(1,ir)+x(2,ir)*x(2,ir)
         accy=accy+y(1,ir)*y(1,ir)+y(2,ir)*y(2,ir) 
         accr=accr+x(1,ir)*y(1,ir)-x(2,ir)*y(2,ir)
         acci=acci+x(1,ir)*y(2,ir)+x(2,ir)*y(1,ir)
      end do
!
      scale=sqrt(accx*accy)
!
      if (dabs(accr)/scale.lt.RNDOFF) accr=ZERO
      if (dabs(acci)/scale.lt.RNDOFF) acci=ZERO
!
      zdotu=dcmplx(accr,acci)
!
      return
      end
