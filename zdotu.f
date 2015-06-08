c     Version 1.5  5/2/94
c*********************************************************************
c
c                    =========================
c                         function  ZDOTU
c                    =========================
c
c       This complex double precision function returns the dot
c       product of the two input vectors without complex 
c       conjugation.
c
c       Calling Parameters:
c       ------------------
c
c       vectors:
c               x    : input vector
c               y    : input vector
c
c       scalars:
c               ndim : dimension of the vectors x and y
c
c
c       Local Variables:
c       ---------------
c               ir   : row counter
c               accr : real part of dot product
c               acci : imaginary part of dot product
c
c       Notes:
c       -----
c               It is assumed that the input vectors are stored in a
c               double precision complex array or as a real double 
c               precision array dimensioned as 2 X mxdim in the 
c               calling routine.
c
c
c       Includes:
c               stddim.inc
c               rndoff.inc
c
c       Uses:
c
c       written by DJS 26-AUG-87
c
c*********************************************************************
c
      function zdotu(x,y,ndim)
c
      complex*16 zdotu
c
      include 'stddim.inc'
      include 'rndoff.inc'
c
      integer ndim
      double precision x,y
      dimension x(2,MXDIM),y(2,MXDIM)
c
      integer ir
      double precision accx,accy,accr,acci,scale
c
c######################################################################
c
      accx=0.0D0
      accy=0.0D0
      accr=0.0D0
      acci=0.0D0
c
      do 10 ir=1,ndim
        accx=accx+x(1,ir)*x(1,ir)+x(2,ir)*x(2,ir)
        accy=accy+y(1,ir)*y(1,ir)+y(2,ir)*y(2,ir) 
        accr=accr+x(1,ir)*y(1,ir)-x(2,ir)*y(2,ir)
        acci=acci+x(1,ir)*y(2,ir)+x(2,ir)*y(1,ir)
 10   continue
c
      scale=sqrt(accx*accy)
c
      if (dabs(accr)/scale.lt.rndoff) accr=0.0D0
      if (dabs(acci)/scale.lt.rndoff) acci=0.0D0
c
      zdotu=dcmplx(accr,acci)
c
      return
      end
