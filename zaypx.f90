c  Version 1.4 10/10/94 
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
c
c**********************************************************************
c
      subroutine zaypx(x,y,scale,ndim)
c
      use nlsdim
c*djs      use rnddbl
c
      integer ndim
      double complex scale
c
      double complex x,y
      dimension x(MXDIM),y(MXDIM)
c
      integer iel
c
c*djs      double complex CZERO
c*djs      parameter (CZERO=(0.0D0,0.0D0))
c
c######################################################################
c
      do iel=1,ndim
        y(iel)=x(iel)+scale*y(iel)
c*djs        if (abs(y(iel)).lt.RNDOFF) y(iel)=CZERO
      end do
c
      return
      end
