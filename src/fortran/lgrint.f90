! NLSL Version 1.5 7/1/95
!----------------------------------------------------------------------
!                    =========================
!                       function LGRINT
!                    =========================
!
! Given an array of five points corresponding to a function fn
! tabulated at equally spaced abscissae x(-2),x(-1),x(0),x(1),x(2),
! such that there is a maximum in fn at one of the internal point
! of the array, this routine finds a value p such that fn(x(0)+p)
! is the extremum, where p is given as a multiple of the spacing
! between the x-values.
!
! The minimization is based on the 5-point Lagrange formula for
! polynomial interpolation using equally-spaced x-values [formula
! 25.2.15 in Abramowitz and Stegun, Handbook of Mathematical Functions,
! Dover, 2nd edition (1978).]
!
! The extremum is found by setting the first derivative of formula
! 25.2.15 w.r.t. p equal to zero and finding the root of the 
! resulting polynomial using a method based on the polynomial
! root-finding algorithm given in Numerical Recipes by Press. et al.
!
! The starting point for the root search is assumed to be x(0) (p=0).
!
! Inputs:
!   fn   5-vector with function values
!
! Output:
!   err  is set to error estimate given by Abramowitz and Stegun
!   Function return is value of p corresponding to the function
!   maximum
!
!----------------------------------------------------------------------
!
      function lgrint(fn)
      implicit none
      double complex lgrint
      double precision fn(5)
!
      integer iter,j,m
      double precision ctmp,eps,yy,zz
      double complex a(4)
      double complex x,dx,x1,b,d,f,g,h,sq,gp,gm,g2,xx,ww
!
      integer MAXIT
      double precision TINY
      double complex CZERO
      parameter (CZERO=(0.,0.),TINY=1.e-15,MAXIT=100)
!
! ---------------------------------------------------------------
! Define coefficients for polynomial that is the first-derivative
! of the 5-point Lagrange polynomial with respect to p
! ---------------------------------------------------------------
      ctmp = ((fn(1)-fn(5))-8.0d0*(fn(2)-fn(4)))/12.0d0
      a(1) = dcmplx(ctmp,0.0d0)
      ctmp = ( -(fn(1)+fn(5))+16.0d0*(fn(2)+fn(4))-30.0d0*fn(3))/12.0d0
      a(2) = dcmplx(ctmp,0.0d0)
      ctmp =(-2.0d0*(fn(1)-fn(5))+6.0d0*(fn(2)-fn(4)))/12.0d0
      a(3) = dcmplx(ctmp,0.0d0)
      ctmp = ( 2.0*(fn(1)+fn(5))-8.0d0*(fn(2)+fn(4)))/12.0d0+fn(3)
      a(4) = dcmplx(ctmp,0.0d0)
!
      x=CZERO
      eps=1.0d-6
      m=3
!
!   --------------------------------------------------------------------
!   Start of polynomial root-finding 
!   Given the degree m and the m+1 complex coefficients a of the
!   fitting polynomial, and given eps, the desired fractional accuracy,
!   and given a complex value x as the seed for the search, the
!   following routine improves x by Laguerre's method until it converges
!   to a root of the polynomial.
!   --------------------------------------------------------------------
      do iter=1,MAXIT
         b=a(m+1)
         d=CZERO
         f=CZERO
         do j=m,1,-1
            f=x*f+d
            d=x*d+b
            b=x*b+a(j)
         end do
!
         if(abs(b).le.TINY) then
            dx=CZERO
         else if(abs(d).le.TINY.and.abs(f).le.TINY)then
            dx=cmplx(abs(b/a(m+1))**(1.d0/m),0.d0)
         else
            g=d/b
            g2=g*g
            h=g2-2.*f/b
            xx=(m-1)*(m*h-g2)
!           yy=abs(real(xx))
!           zz=abs(aimag(xx))
            yy=abs(DBLE(xx))
            zz=abs(DIMAG(xx))
            if(yy.lt.TINY.and.zz.lt.TINY) then
               sq=CZERO
            else if (yy.ge.zz) then
               ww=(1.0/yy)*xx
               sq=sqrt(yy)*sqrt(ww)
            else
               ww=(1.0/zz)*xx
               sq=sqrt(zz)*sqrt(ww)
            endif
            gp=g+sq
            gm=g-sq
            if(abs(gp).lt.abs(gm)) gp=gm
            dx=m/gp
         endif
         x1=x-dx
         if(x.eq.x1)then
            lgrint=x
            return
         endif
!
         x=x1
         if(abs(dx).le.eps*abs(x))then
            lgrint=x
            return
         endif
      end do
      print *, 'LGRINT: too many iterations'
      print *, '[execution paused, press enter to continue]'
      read (*,*)
      return
      end
