! Version 1.5.1 beta 2/3/96
!----------------------------------------------------------------------
!     This file contains various functions for making statistical
!     inferences from the residuals and curvature matrix returned
!     by the least-squares procedure. (Subroutines marked with *)
!
!       chi       Returns chi-squared for given confidence, degrees freedom
!       pchi      Chi-squared probability function
!       residx    Returns residual index of fit vs. data
!       corrl     Returns correlation coefficient of fit vs. data
!       gammq     Incomplete gamma function
!     * gcf       Continued-fraction approximation to incomplete gamma fn
!     * gser      Series approximation of incomplete gamma function
!       gammln    Returns ln( gamma )
!     * zbrac     Bracketing routine for Brent root-finder
!       zbrent    Brent's method for finding root of function pchi
!
!     Auxiliary functions:   betai, betacf, alphat, alphaf, ts, fs
!     Auxiliary subroutine:  pausex
!
!----------------------------------------------------------------------
!                    =========================
!                         function CHI
!                    =========================
!     Finds delta-chi for a given confidence level p and number of
!     degrees of freedom, nu.
!
!     This is done by finding the value of the chi-squared probability
!     function, gammq( nu/2, chi^2/2 ) corresponding to the desired
!     confidence level for a given number of degrees of freedom. The
!     value is obtained by finding the root of function pchi (coded
!     below), which returns the confidence level p as a function of 
!     chi. The root is first bracketed using routine zbrac, and the
!     root found by Brent's method (function zbrent).
!----------------------------------------------------------------------
      function chi( p, nu )
      implicit none
      double precision p, chi
      integer nu
!
      double precision x1,x2
      logical success
!
      double precision halfnu, conf
      common /chicom/ halfnu, conf
!
      double precision pchi,zbrent
      external pchi,zbrent
!
      chi=0.0d0
      halfnu=0.5d0*nu
      conf=p
      x1=1.0d0
      x2=2.0d0
      call zbrac( pchi, x1, x2, success )
      if (success) chi=zbrent( pchi, x1, x2, 1.0d-6 )
      return
      end


!----------------------------------------------------------------------
!                    =========================
!                         function PCHI
!                    =========================
!     Returns the chi-squared probability function given chi-squared
!     and the number degrees of freedom, nu.
!----------------------------------------------------------------------
      function pchi(chi)
      implicit none
      double precision pchi,chi
!
      double precision halfnu, conf
      common /chicom/ halfnu, conf
!
      double precision gammq
      external gammq
!
      pchi = gammq( halfnu, 0.5d0*chi ) - 1.0d0 + conf
      return
      end
!----------------------------------------------------------------------
!                    =========================
!                       function RESIDX
!                    =========================
!     Returns the residual index between experimental data
!     and calculated fit function in common blocks /expdat/ and /lmcom/
!----------------------------------------------------------------------
      function residx()
!
      use nlsdim
      use expdat
      use lmcom
!
      implicit none
      integer i,j,k
      double precision osum,residx,rsum,yc,yo
!
      osum=0.0d0
      rsum=0.0d0
      do i=1,nspc
         do j=1,npts(i)
            k=ixsp(i)+j-1
            yo=data(k)
            yc=data(k)-rmsn(i)*fvec(k)
            rsum=rsum+abs(yo-yc)
            osum=osum+abs(yo)
         end do
      end do
      residx=rsum/osum
      return
      end

!----------------------------------------------------------------------
!                    =========================
!                       function CORRL
!                    =========================
!     Returns the correlation function between experimental data
!     and calculated fit function in common blocks /expdat/ and /lmcom/
!----------------------------------------------------------------------
      function corrl()
!
      use nlsdim
      use expdat
      use lmcom
!
      implicit none
      integer i,j,k
      double precision cbar,cdsum,c2sum,dbar,dc,dd,d2sum,corrl
!
      dbar=0.0d0
      cbar=0.0d0
      do i=1,nspc
         do j=1,npts(i)
            k=ixsp(i)+j-1
            dbar=dbar+data(k)
            cbar=cbar+(data(k)-rmsn(i)*fvec(k))
         end do
      end do
      dbar=dbar/float(ndatot)
      cbar=cbar/float(ndatot)
!
      cdsum=0.0d0
      d2sum=0.0d0
      c2sum=0.0d0
      do i=1,nspc
         do j=1,npts(i)
            k=ixsp(i)+j-1
            dd=(data(k)-dbar)
            dc=(data(k)-rmsn(i)*fvec(k)-cbar)
            cdsum=cdsum+dd*dc
            d2sum=d2sum+dd*dd
            c2sum=c2sum+dc*dc
         end do
      end do
      corrl=cdsum/sqrt(d2sum*c2sum)
!
      return
      end
!----------------------------------------------------------------------
!                        =========================
!                             function BETAI
!                        =========================
!     Incomplete beta function, used in calculation of t-distribution
!     and F-distribution used in statistical analysis of parameters
!     in NLS fits.
!
!     Function returns the function I_x(A,B)
!      
!     From Numerical Recipes by W.H. Press et al.
!----------------------------------------------------------------------
!
      function betai(a,b,x)
      implicit none
      double precision betai,a,b,x,bt
!
      double precision ONE,TWO,ZERO
      parameter(ONE=1.0d0,TWO=2.0d0,ZERO=0.0D0)
!
      double precision betacf,gammln
      external betacf,gammln
!
      if (x.lt.ZERO .or. x.gt.ONE) then
         call pausex('bad argument X in BETAI')
      end if
      if (x.eq.ZERO .or. x.eq.ONE) then
         bt=ZERO
      else
!
!        ----------------------------------------------------
!        Calculate factors in front of the continued fraction 
!        ----------------------------------------------------
         bt=dexp( gammln(a+b)-gammln(a)-gammln(b)
     #          +a*dlog(x)+b*dlog(ONE-x) )
      end if
!
!     --------------------------------
!     Use continued fraction directly
!     --------------------------------
      if (x.lt.(a+ONE)/(a+b+TWO)) then
         betai=bt*betacf(a,b,x)/a
         return
!
!     ----------------------------------------------------------------
!     Use continued fraction after making the symmetry transformation  
!     ----------------------------------------------------------------
      else
         betai=ONE-bt*betacf(b,a,ONE-x)/b
         return
      end if
      end


      function betacf(a,b,x)
      implicit none
      double precision betacf,a,b,x
!
      integer m
      double precision am,aold,ap,app,az,bp,bm,bpp,bz,d,em,qab,qap,
     #                 qam,tem
!
      integer ITMAX
      double precision ONE,TWO,ZERO,EPS
      parameter(ITMAX=100,EPS=3.0d-7,ONE=1.0d0,TWO=2.0d0,ZERO=0.0D0)
!
      am=ONE
      bm=ONE
      az=ONE
      qab=a+b
      qap=a+ONE
      qam=a-ONE
      bz=ONE-qab*x/qap
      do m=1,ITMAX
         EM=m
         tem=em+em
!
!        ----------------------------
!        Even step of the recurrence
!        ----------------------------
         d=em*(b-m)*x/((qam+tem)*(a+tem))
         ap=az+d*am
         bp=bz+d*bm
!
!        ----------------------------
!        Odd step of the recurrence
!        ----------------------------
         d=-(a+em)*(qab+em)*x/((a+tem)*(qap+tem))
         app=ap+d*az
         bpp=bp+d*bz
!
!        ----------------------------------------------------
!        Save old answer and renormalize to prevent overflow
!        ----------------------------------------------------
         aold=az
         am=ap/bpp
         bm=bp/bpp
         az=app/bpp
         bz=ONE
!
!        ----------------------
!        Check for convergence
!        ----------------------
         if (abs(az-aold).lt.EPS*abs(az)) goto 1
      end do
!
      write(*,1000) abs((az-aold)/az)
 1000 format('BETACF: ITMAX exceeded, final error was ',g11.3)
!
    1 betacf=az
      return
      end

      function alphat(t)
      implicit none
      double precision alphat,t
!
      double precision a,b,alpha
      common /bcom/ a,b,alpha
!
      double precision betai
      external betai
!
      alphat = alpha - 0.5*betai( a, b, a/(a+0.5D0*t*t) )
      return
      end


      function alphaf(f)
      implicit none
      double precision alphaf,f
!
      double precision a,b,alpha
      common /bcom/ a,b,alpha
!
      double precision betai
      external betai
!
      alphaf = betai( b, a, b/(b+a*f) ) - alpha
      return
      end

      function ts( al, nu )
      implicit none
      double precision ts, al
      integer nu
!
      double precision x1,x2
      logical success
!
      double precision a, b, alpha
      common /bcom/ a, b, alpha
!
      double precision alphat,zbrent
      external alphat,zbrent
!
      alpha=al
      a=0.5d0*nu
      b=0.5d0
      x1=0.0d0
      x2=1.0d0
      call zbrac( alphat, x1, x2, success )
      if (success) ts=zbrent( alphat, x1, x2, 1.0d-6 )
      return
      end


      function fs( al, nu1, nu2 )
      implicit none
      double precision fs, al
      integer nu1,nu2
!
      double precision x1,x2
      logical success
!
      double precision a, b, alpha
      common /bcom/ a, b, alpha
!
      double precision alphaf,zbrent
      external alphaf,zbrent
!
      alpha=al
      a=0.5D0*nu1
      b=0.5D0*nu2
      x1=0.0d0
      x2=1.0d0
      call zbrac( alphaf, x1, x2, success )
      if (success) fs=zbrent( alphaf, x1, x2, 1.0d-6 )
      return
      end

!----------------------------------------------------------------------
!                    =========================
!                        function GAMMQ
!                    =========================
!
!  Incomplete gamma function. From Numerical Recipes by Press, et al.
!
!----------------------------------------------------------------------
!
      FUNCTION gammq(a,x)
      IMPLICIT none
      DOUBLE PRECISION a,gammq,x
CU    USES gcf,gser
      DOUBLE PRECISION gammcf,gamser,gln
      if(x.lt.0..or.a.le.0.) call pausex('bad arguments in gammq')
      if(x.lt.a+1.)then
        call gser(gamser,a,x,gln)
        gammq=1.-gamser
      else
        call gcf(gammcf,a,x,gln)
        gammq=gammcf
      endif
      return
      END
!
!----------------------------------------------------------------------
!                    =========================
!                        function GCF
!                    =========================
!
!  Continued-fraction calculation of incomplete gamma function. 
!  From Numerical Recipes by Press, et al.
!
!----------------------------------------------------------------------
      SUBROUTINE gcf(gammcf,a,x,gln)
      IMPLICIT none
      INTEGER ITMAX
      DOUBLE PRECISION a,gammcf,gln,x,EPS,FPMIN
      PARAMETER (ITMAX=100,EPS=3.e-7,FPMIN=1.e-30)
C     USES gammln
      INTEGER i
      DOUBLE PRECISION an,b,c,d,del,h,gammln
      gln=gammln(a)
      b=x+1.-a
      c=1./FPMIN
      d=1./b
      h=d
      do 11 i=1,ITMAX
        an=-i*(i-a)
        b=b+2.
        d=an*d+b
        if(abs(d).lt.FPMIN)d=FPMIN
        c=b+an/c
        if(abs(c).lt.FPMIN)c=FPMIN
        d=1./d
        del=d*c
        h=h*del
        if(abs(del-1.).lt.EPS)goto 1
11    continue
      call pausex('a too large, ITMAX too small in gcf')
1     gammcf=exp(-x+a*log(x)-gln)*h
      return
      END

!----------------------------------------------------------------------
!                    =========================
!                        function GSER
!                    =========================
!
!  Series approximation of incomplete gamma function. 
!  From Numerical Recipes by Press, et al.
!
!----------------------------------------------------------------------
      SUBROUTINE gser(gamser,a,x,gln)
      IMPLICIT none
      INTEGER ITMAX
      DOUBLE PRECISION a,gamser,gln,x,EPS
      PARAMETER (ITMAX=100,EPS=3.e-7)
C     USES gammln
      INTEGER n
      DOUBLE PRECISION ap,del,sum,gammln
      gln=gammln(a)
      if(x.le.0.)then
        if(x.lt.0.) call pausex('x < 0 in gser')
        gamser=0.
        return
      endif
      ap=a
      sum=1./a
      del=sum
      do 11 n=1,ITMAX
        ap=ap+1.
        del=del*x/ap
        sum=sum+del
        if(abs(del).lt.abs(sum)*EPS)goto 1
11    continue
      call pausex('a too large, ITMAX too small in gser')
1     gamser=sum*exp(-x+a*log(x)-gln)
      return
      END
!
!----------------------------------------------------------------------
!                    =========================
!                        function GAMMLN
!                    =========================
!
!  Returns the value ln(gamma(xx) for xx>1, with full accuracy for xx>1
!  From Numerical Recipes by Press, et al.
!
!----------------------------------------------------------------------
      FUNCTION gammln(xx)
      IMPLICIT none
      DOUBLE PRECISION gammln,xx
      INTEGER j
      DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)
      SAVE cof,stp
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,
     *24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,
     *-.5395239384953d-5,2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.000000000190015d0
      do 11 j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y
11    continue
      gammln=tmp+log(stp*ser/x)
      return
      END

!
!######################################################################
!
      SUBROUTINE zbrac(func,x1,x2,succes)
      IMPLICIT none
      INTEGER NTRY
      DOUBLE PRECISION x1,x2,func,FACTOR
      EXTERNAL func
      PARAMETER (FACTOR=1.6,NTRY=50)
      INTEGER j
      DOUBLE PRECISION f1,f2
      LOGICAL succes
      if (x1.eq.x2)
     *  call pausex('you have to guess an initial range in zbrac')
      f1=func(x1)
      f2=func(x2)
      succes=.true.
      do 11 j=1,NTRY
        if(f1*f2.lt.0.)return
        if(abs(f1).lt.abs(f2))then
          x1=x1+FACTOR*(x1-x2)
          f1=func(x1)
        else
          x2=x2+FACTOR*(x2-x1)
          f2=func(x2)
        endif
11    continue
      succes=.false.
      return
      END
      FUNCTION zbrent(func,x1,x2,tol)
      IMPLICIT none
      INTEGER ITMAX
      DOUBLE PRECISION zbrent,tol,x1,x2,func,EPS
      EXTERNAL func
      PARAMETER (ITMAX=100,EPS=3.e-8)
      INTEGER iter
      DOUBLE PRECISION a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
      a=x1
      b=x2
      fa=func(a)
      fb=func(b)
      e=0.0d0
      if((fa.gt.0..and.fb.gt.0.).or.(fa.lt.0..and.fb.lt.0.))
     *  call pausex('root must be bracketed for zbrent')
      c=b
      fc=fb
      do 11 iter=1,ITMAX
        if((fb.gt.0..and.fc.gt.0.).or.(fb.lt.0..and.fc.lt.0.))then
          c=a
          fc=fa
          d=b-a
          e=d
        endif
        if(abs(fc).lt.abs(fb)) then
          a=b
          b=c
          c=a
          fa=fb
          fb=fc
          fc=fa
        endif
        tol1=2.*EPS*abs(b)+0.5*tol
        xm=.5*(c-b)
        if(abs(xm).le.tol1 .or. fb.eq.0.)then
          zbrent=b
          return
        endif
        if(abs(e).ge.tol1 .and. abs(fa).gt.abs(fb)) then
          s=fb/fa
          if(a.eq.c) then
            p=2.*xm*s
            q=1.-s
          else
            q=fa/fc
            r=fb/fc
            p=s*(2.*xm*q*(q-r)-(b-a)*(r-1.))
            q=(q-1.)*(r-1.)*(s-1.)
          endif
          if(p.gt.0.) q=-q
          p=abs(p)
          if(2.*p .lt. min(3.*xm*q-abs(tol1*q),abs(e*q))) then
            e=d
            d=p/q
          else
            d=xm
            e=d
          endif
        else
          d=xm
          e=d
        endif
        a=b
        fa=fb
        if(abs(d) .gt. tol1) then
          b=b+d
        else
          b=b+sign(tol1,xm)
        endif
        fb=func(b)
11    continue
      call pausex('zbrent exceeding maximum iterations')
      zbrent=b
      return
      END

      subroutine pausex(mesg)
      character*(*) mesg
      print *, mesg
      print *, '[execution paused, press enter to continue]'
      read (*,*)
      return
      end
