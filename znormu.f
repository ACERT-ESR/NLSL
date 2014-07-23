c  VERSION 1.0  (NLSPMC version)   2/5/99
c**********************************************************************
c
c                                ZNORMU
c                                ======
c
c       ZNORMU is a complex double precision function for calculating 
c       the Euclidean pseudo-norm of a complex vector, i.e. the square 
c       root of the sum of the squares of the components, NOT the 
c       square root of the sum of the moduli of the components.
c
c       Calling Parameters:
c       ------------------
c
c       vectors:
c
c               v       : a NDIM dimensional complex vector whose
c                         Euclidean pseudo-norm is desired.
c
c       scalars:
c
c               ndim    : number of elements in the vector v.
c
c
c       Local Variables:
c       ---------------
c       
c       vectors:
c
c       scalars:
c
c               ir      : index for loop over components of v
c               amp     : amplitude of pseudo-norm of v.
c               phase   : phase of pseudo-norm of v.
c               tempr   : temporarily stores the real part of the 
c                         pseudo-norm when checking for roundoff.
c               tempi   : temporarily stores the real part of the
c                         pseudo-norm when checking for roundoff.
c               accum   : holds the partial sums in loop over the 
c                         in the loop over the components of v.
c
c       Notes:
c       ------
c               1)  The size of the vector v is determined by the 
c                   parameter MXDIM in the include file stddim.inc.
c                   In the calling routine, the vector v should be 
c                   dimensioned as complex*16 v(mxdim), or 
c                   alternatively as real*8 v(2,mxdim).  
c                   
c               2)  The square root of the final value of accum is 
c                   not calculated using the usual Fortran square root
c                   function.  The methods used here allows one to 
c                   take the square root of a number on the negative
c                   real axis.
c
c               3)  The standard include files are used to specify
c                   the dimensions of arrays, the unit roundoff error
c                   of the machine, etc.
c
c       written by DJS 5-12-87          
c
c       Includes:
c               nlsdim.inc
c               rndoff.inc
c
c       Uses:
c
c**********************************************************************
c
      function znormu(v,ndim)
c
      complex*16 znormu
c
      include 'limits.inc'
      include 'rndoff.inc'
c
      integer ndim
      double precision v
      dimension v(2,mxdim)
c
      integer ir
      double precision amp,phase,accv,accr,acci
c
c######################################################################
c
      accv=0.0D0
      accr=0.0D0
      acci=0.0D0
c
      do 10 ir=1,ndim
        accv=accv+v(1,ir)*v(1,ir)+v(2,ir)*v(2,ir)
        accr=accr+v(1,ir)*v(1,ir)-v(2,ir)*v(2,ir)
        acci=acci+2.0D0*v(1,ir)*v(2,ir)
 10   continue
c
      amp=sqrt(sqrt(accr*accr+acci*acci))
      phase=0.5D0*datan2(acci,accr)
c
      if (amp/sqrt(accv).lt.rndoff) amp=0.0D0
c
      accr=amp*cos(phase)
      acci=amp*sin(phase)
c
      znormu=dcmplx(accr,acci)
c
      return
      end
