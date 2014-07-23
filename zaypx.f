c  NLSPMC VERSION 1.0   2/5/99
c**********************************************************************
c
c             complex double precision vector scale and add
c             ---------------------------------------------
c
c        This subroutine will update a complex double precision
c        vector Y by the formula,
c
c                         Y=x+a*Y ,
c
c        where a is a complex double precision constant and X is a
c        complex double presion vector.
c
c       Includes:
c               nlsdim.inc
c               rndoff.inc
c               
c       Uses:
c
c**********************************************************************
c
      subroutine zaypx(x,y,scale,ndim)
c
      include 'limits.inc'
      include 'rndoff.inc'
c
      integer ndim
      complex*16 scale
c
      complex*16 x,y
      dimension x(mxdim),y(mxdim)
c
      integer iel
      complex*16 czero
      parameter (czero=(0.0D0,0.0D0))
c
c######################################################################
c
      do 10 iel=1,ndim
        y(iel)=x(iel)+scale*y(iel)
c*djs        if (abs(y(iel)).lt.rndoff) y(iel)=czero
 10   continue
c
      return
      end
