c  VERSION 1.0  (NLSPMC version)   2/5/99
c----------------------------------------------------------------------
c                     =========================
c                         function GETDAT
c                     =========================
c
c Inputs 2D experimental EPR spectra for nonlinear least squares
c fitting programs.
c
c      function getd2d( filen,infmt,npt1,npt2,lumsg,iexp,icomb,
c                       stept1,stept2,spndat,cwdat,
c                       x1a,x2a,tmp1,tmp2,ydat,y2drv )
c
c Returns 0 for a successful read.
c
c Subroutine arguments:
c
c   filen       Name of data file
c   infmt       Specifies the input format of the datafile
c               0: ASCII format for stap or conp (DISSPLA 3D PLOT)
c               1: binary format (IBM in eagle, IEEE elsewhere)
c   npt1,npt2   Number of data points in 2D spectrum
c   lumsg       Logical unit for output of informational messages.
c               Set to 0 if no output is desired
c   iexp        Experimental type
c   icomb       Combination type
c   stept1      Step size in t1 specified in series option
c   stept2      Step size in t2 specified in series option
c   x1a,x2a,tmp1,tmp2,ydat,y2drv  Work spaces for spline
c
c Output:
c   smax        Maximum of the input data
c   spndat      2D splined spectrum
c   cwdat       cw-equivalent spectrum extracted from 2D spectrum
c
c----------------------------------------------------------------------
c
      function getd2d( filen,infmt,npt1,npt2,lumsg,iexp,icomb,
     #                 spt1,spt2,smax,spndat,cwdat,x1a,x2a,tmp1,
     #                 tmp2,ydat,y2drv )
      implicit none
c
      include 'limits.inc'
      include 'stdio.inc'
c
      character*30 filen
      integer getd2d,infmt,npt1,npt2,lumsg,iexp,icomb,i,j,k,
     #        npt1i,npt2i,nend,ibuf
      double precision spndat(npt1,npt2),cwdat(npt2),spt1,spt2,
     #        inif1,res1,inif2,res2,inif1i,res1i,inif2i,res2i,f1,f2,
     #        smin,smax
c
      double precision x1a(mxegv),x2a(mxegv),tmp1(mxegv),tmp2(mxegv),
     #        tmp3(mxegv),ydat(mxegv,mxegv),y2drv(mxegv,mxegv)
      integer imx,jmx
      double precision mxdat
c
c------ for binary input format
c
c      real wksp1(npt1,npt2)
c
c######################################################################
c
c----------------------------------------------------------------------
c     Input datafile according to specified format:
c----------------------------------------------------------------------
      if (lumsg.ne.0) write(lumsg,1005) filen
c
c      if (infmt .eq. 1) then
c
c   .....format #1: binary single precision format : NOT ACTIVE OPTION
c
c         open (unit=ludisk,file=filen,status='old',
c     #         access='sequential',form='unformatted',err=98)
c         read (ludisk,err=98) npt2i,npt1i
c         read (ludisk,err=98) inif2i,res2i,inif1i,res1i
c         read (ludisk,err=98) wksp1 
c         close (unit=ludisk)
c
c     ** convert single precision into double precision
c
c         do 10 i=1,npt1i
c         do 10 j=1,npt2i
c            ydat(i,j)=dble(wksp1(i,j))
c 10      continue
c      else
c
c  ......format #0 (default): standard ASCII format
c
         open (unit=ludisk,file=filen,status='old',
     #         access='sequential',form='formatted',err=98)
         read (ludisk,*,err=98)
         read (ludisk,1110,err=98) npt2i
         read (ludisk,1110,err=98) npt1i
         read (ludisk,1120,err=98) inif2i,res2i,inif1i,res1i,smin,smax
         read (ludisk,1130,err=98) ((ydat(i,j),j=1,npt2i),i=1,npt1i)
         close (unit=ludisk)
c
c      end if
c
      if ((npt1i.gt.mxegv).or.(npt2i.gt.mxegv)) then
         getd2d=-2
         return
      else
         go to 100
      end if
c
c     ** Abort processing if there was an input error
 98   close (unit=ludisk)
      getd2d=-1
      return
c
c----------------------------------------------------------------------
c     Perform 2D natural bicubic spline interpolation to get desired
c     resolution & number of points
c     (Ref. Numerical Recipes : spline,splint,splie2,splin2 routines)
c----------------------------------------------------------------------
c
 100  inif1=-5.0d2*dble(npt1-1)/spt1/dble(npt1)
      res1=-2.0d0*inif1/dble(npt1-1)	! step size in f1
      inif1=-5.0d2/spt1
      inif2=-5.0d2*dble(npt2-1)/spt2/dble(npt2)
      res2=-2.0d0*inif2/dble(npt2-1)
      inif2=-5.0d2/spt2
c
      if ( ((inif1-inif1i).lt.-1.0d0).or.((inif2-inif2i).lt.-1.0d0) )
     #                                        then
         getd2d=-3
         return
      end if
c
      if ((npt1i.ne.npt1).or.(npt2i.ne.npt2).or.
     #    (abs(inif1-inif1i).gt.1.0d-3).or.
     #    (abs(inif2-inif2i).gt.1.0d-3).or.
     #    (abs(res1-res1i).gt.1.0d-3).or.
     #    (abs(res2-res2i).gt.1.0d-3) ) then	! do spline...
c
         do 110 j=1,npt1i
            x1a(j)=inif1i+dble(j-1)*res1i
 110     continue
         do 112 j=1,npt2i
            x2a(j)=inif2i+dble(j-1)*res2i
 112     continue
c
c                         ** construct 2nd derivative table
         do 120 j=1,npt1i
            do 122 k=1,npt2i
               tmp1(k)=ydat(j,k)
 122        continue
            call ncspln(x2a,tmp1,npt2i,tmp2)
            do 124 k=1,npt2i
               y2drv(j,k)=tmp2(k)
 124        continue
 120     continue
c
c                         ** loop over f2 values
         do 130 i=1,npt2
            f2=inif2+dble(i-1)*res2
c
            do 132 j=1,npt1i
               do 134 k=1,npt2i
                  tmp1(k)=ydat(j,k)
                  tmp2(k)=y2drv(j,k)
 134           continue
               call splnay(x2a,tmp1,tmp2,npt2i,f2,1.0d0,tmp3(j),1)
 132        continue
c
            call ncspln(x1a,tmp3,npt1i,tmp2)
            call splnay(x1a,tmp3,tmp2,npt1i,inif1,res1,
     #                  spndat(1,i),npt1)
 130     continue
c
      else
c                 ** copy the spectrum
         do 200 j=1,npt1
         do 200 k=1,npt2
            spndat(j,k)=ydat(j,k)
 200     continue
      end if
c
c----------------------------------------------------------------------
c     Extract cw-equivalent spectrum from 2D spetrum
c----------------------------------------------------------------------
c
      if ((iexp.eq.1).or.(iexp.eq.3)) then
c
c *** FID experiment ***
c
c     The following restriction is only for obtaining quick working
c     version of NLSPMC program.  One should modify those codes keeping
c     in mind that correct information should be passed for the
c     calculation of cw-spectrum.  (init2, spt2, npt2 are used in
c     current version.  See calling sequence for xshft in pfunnew.f for
c     detail.)
c         1) inif1=inif2
c         2) res1=res2
c         3) npt1=npt2
c
        if ((npt1.ne.npt2).or.(abs(res1-res2).gt.1.0d-3).or.
     #      (abs(inif1-inif2).gt.1.0d-3)) then
           getd2d=-4
           return
        end if
c
        if (icomb.eq.1) then
c                              * positive diagonal
            do 210 i=1,npt2
 210           cwdat(i)=spndat(i,i)
        else
c                              * negative diagonal
c            mxdat=0.0D0
c            do 211 i=1,npt1
c              do 211 j=1,npt2
c                if(spndat(j,i).gt.mxdat) then
c                  mxdat=spndat(j,i)
c                  jmx=j
c                  imx=i
c                end if
c211         continue
c
c modified to pick up expected max at (65,65) for 128x128 data file.
c
            cwdat(1)=0.0d0
            do 212 i=2,npt2
              cwdat(i)=spndat(npt2-i+2,i)
 212        continue
c        stop
        end if      
      else if ((iexp.eq.2).or.(iexp.eq.4).or.(iexp.eq.5)) then
c
c *** ECHO experiment ***
c
         j=0
 220     j=j+1
         f1=inif1+(j-1)*res1
         if (f1.lt.0.0d0) go to 220
         if (-(f1-res1).lt.f1) j=j-1
c
         do 222 i=1,npt2
 222        cwdat(i)=spndat(j,i)
      end if
c
      getd2d=0
      return
c
c### format statements ############################################
c
 1005 format(/5x,'Opening file ',a,/)
 1110 format(i10)
 1120 format(6e10.4)
 1130 format(8e10.4)
c
      end


c----------------------------------------------------------------------
c                       =======================
c                         subroutine NCSPLN
c                       =======================
c Perform a natural cubic spline fit. This is a modification of 
c subroutine SPLINE from Numerical Recipes. Given the tabulated function
c of N points in arrays X and Y, return an array Y2 of length N which 
c contains the second derivative of the interpolating function at the 
c corresponding X values. The major difference from the original SPLINE 
c routine is that the second derivatives at the both boundaries of 
c the function are assumed to be zero (hence "natural cubic spline").
c----------------------------------------------------------------------
      subroutine ncspln(x,y,n,y2)
      implicit none
      integer nmax
      parameter(nmax=512)
      double precision zero,one,two,six
      parameter(zero=0.0d0,one=1.0d0,two=2.0D0,six=6.0d0)
c
      integer i,k,n
      double precision x(n),y(n),y2(n),u(nmax),sig,p
c      
      y2(1)=zero
      u(1) =zero
      do 18 i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1) + two
        y2(i)=(sig-one)/p
        u(i)=(six*( (y(i+1)-y(i))/(x(i+1)-x(i))
     &               -(y(i)-y(i-1))/(x(i)-x(i-1))  )
     &              /(x(i+1)-x(i-1))
     &         - sig*u(i-1))/p
18    continue
      y2(n)=zero
      do 19 k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
19    continue
      return
      end


c----------------------------------------------------------------------
c                         ====================
c                          subroutine SPLNAY
c                         ====================
c
c Modification of Numerical Recipes spline interpolation routine (SPLINT)
c which interpolates an entire array at once. Given the original N 
c data points in XA,YA and the second derivative YA2 (calculated by SPLINE
c or NCSPLN) and the set of NS X values starting at X0 and spaced by DX, 
c this routine returns the interpolated function values in Y.
c The original binary search in subroutine SPLINT has been replaced, 
c since the interpolations start near the lower bound and are repeated 
c for successively higher, but closely spaced x values.
c----------------------------------------------------------------------
      subroutine splnay(xa,ya,y2a,n,x0,dx,y,ns)
      implicit none
      integer i,klo,khi,n,ns
      double precision xa(n),ya(n),y2a(n),y(ns),x,x0,dx,h,a,b
c
      if (x0 .lt. xa(1) .or. x0+dx*(ns-1) .gt. xa(n)) 
     #    pause 'SPLNAY: Bad interpolated X range.'
      x=x0-dx
      klo=1
      khi=2
      do 3 i=1, ns
        x=x+dx
2       if (x.gt.xa(khi) .and. khi.lt.n) then
          klo=khi
          khi=khi+1
        goto 2
        endif
        h=xa(khi)-xa(klo)
        if (h .eq. 0.0D0) pause 'SPLNAY: Bad XA input.'
        a=(xa(khi)-x)/h
        b=(x-xa(klo))/h
        y(i)=a*ya(klo) + b*ya(khi) +
     *        ((a**3-a)*y2a(klo) + (b**3-b)*y2a(khi))*(h**2)/6.0D0
3     continue
      return
      end
