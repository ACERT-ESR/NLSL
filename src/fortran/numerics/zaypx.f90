!  Version 1.4 10/10/94 
!**********************************************************************
!
!             complex double precision vector scale and add
!             ---------------------------------------------
!
!        This subroutine will update a complex double precision
!        vector Y by the formula,
!
!                         Y=x+a*Y ,
!
!        where a is a complex double precision constant and X is a
!        complex double presion vector.
!
!
!**********************************************************************
!
      subroutine zaypx(x,y,scale,ndim)
!
      use nlsdim
!*djs      use rnddbl
!
      integer ndim
      double complex scale
!
      double complex x,y
      dimension x(MXDIM),y(MXDIM)
!
      integer iel
!
!*djs      double complex CZERO
!*djs      parameter (CZERO=(0.0D0,0.0D0))
!
!######################################################################
!
      do iel=1,ndim
        y(iel)=x(iel)+scale*y(iel)
!*djs        if (abs(y(iel)).lt.RNDOFF) y(iel)=CZERO
      end do
!
      return
      end
