! NLS Version 1.5.1 beta 11/25/95
!----------------------------------------------------------------------
!                     =========================
!                         function GETDAT
!                     =========================
!
! Inputs experimental EPR spectra for nonlinear least squares fitting
! programs.
!
!      function getdat( filen,infmt,modebc,nspln,ideriv,norm,lumsg,
!     #                 comment,ncmts,spndat,bi,delb,stddev,xraw,yraw,y2)
!
! Returns 0 for a successful read.
!
! Subroutine arguments:
!
!   filen    Name of data file
!   infmt    Specifies the input format of the datafile
!            0: PLOT/GENPLOT ASCII format with one X,Y pair per line
!               (lines with 'c', 'C', or '!' in first column skipped)
!            1: PC-GENPLOT binary format files
!   modebc   Specifies the type of baseline correction to be applied
!            to the data, as follows:
!            0: No baseline correction
!            n: Fit a least squares line to the first and last <n> points
!               and subtract it from the spectrum
!   nspln    Number of points desired in the final splined data
!            (0 if no spline is desired)
!   ideriv   Derivative flag (nonzero for first derivative data)
!   norm     Normalization flag: 0=do not normalize data, otherwise,
!            adjust so that integral=1 (2nd integral for 1st deriv data) 
!   lumsg    Logical unit for output of informational messages.
!            Set to 0 if no output is desired
!   xraw     Work arrays to contain raw x,y data. Should be dimensioned
!   yraw     to maximum expected number of points
!   y2       Work array to contain spline data
!
! Output:
!   comment  File comment
!   ncmts     Number of comment lines in file
!   spndat   Splined data after baseline correction
!   nspln    Number of data points if input {nspln,bi,delb} incomplete
!   bi       Low field of input data
!   delb     Delta-field per point
!
! Functions coded in this file:
!
!    getasc    get ASCII 2-column data
!    sglint    Return single integral of data array
!    dblint    Return double integral of data array
!    linebc    Baseline-correct of data array
!    ncspln    Determine spline coefficients from array
!    splnay    Spline interpolation of array
!    normlz    Normalize array to unit integral
!
! Calls functions (not coded in this file)  
!    getbin    (in genio.c)
!
!----------------------------------------------------------------------
!
      function getdat( filen,infmt,modebc,nspln,ideriv,norm,lumsg,
     #                 comment,ncmts,spndat,bi,delb,stddev,xraw,yraw,
     #                 y2)
!
!   --- The following is needed for the definition of "MXCMT" and "MXINP"
!
      use nlsdim
!
      implicit none
      integer getdat,ideriv,infmt,lumsg,modebc,norm,nspln,ncmts
      double precision bi,delb,stddev,spndat(nspln),xraw(MXINP),
     #                 yraw(MXINP),y2(MXINP)
!
      character*30 filen
      character*80 comment(MXCMT)
!
      integer i,iret,npts
      double precision anorm,c0,c1
!
      integer getasc,getbin

      double precision sglint,dblint
      external getasc,getbin,sglint,dblint
!
!######################################################################
!
!----------------------------------------------------------------------
!       Input datafile according to specified format:
!----------------------------------------------------------------------
      if (lumsg.ne.0) write(lumsg,1005) filen
!
!.......format #1: PC-GENPLOT binary data (input routine written in C)
!
      if (infmt .eq. 1) then
        iret=getbin(filen,xraw,yraw,npts,comment,ncmts,MXINP)
      else
!
!......format #0 (default): standard ASCII format
!
        iret=getasc(filen,xraw,yraw,npts,comment,ncmts,MXCMT,MXINP)
      end if
!
!----------------------------------------------------------------------
!     Abort processing if there was an input error
!----------------------------------------------------------------------
      if (iret .ne. 0) then
         getdat=1
         return
      else
!                                     Too many points
         if (npts.gt.MXINP) then
            if (lumsg.ne.0) write (lumsg,1013) MXINP,npts-MXINP
            npts=MXINP
         end if
!                                     Report # points read
         if (lumsg.ne.0) then
            if (infmt.eq.1) then
               write(lumsg,1006) npts,xraw(1),xraw(npts)
            else
               write(lumsg,1007) npts,xraw(1),xraw(npts)
            end if
         end if
!
!                                     Check field range
!
         if ((xraw(1).ge.xraw(npts)) .or. 0.5d0*(xraw(1)+xraw(npts))
     #        .lt. abs(xraw(npts)-xraw(1)) ) then
            if (lumsg.ne.0) write(lumsg,1012)
         end if
!
      end if
!
!----------------------------------------------------------------------
!   Baseline correction.
!
!      Fit a line to first and last MODEBC points of the curve and 
!      subtract it from the data.
!----------------------------------------------------------------------
      stddev = 0.0d0
      call linebc( xraw,yraw,npts,modebc,c0,c1,stddev )
      if (modebc.ne.0) then
         if (lumsg.ne.0) write (lumsg,1008) modebc,c0,c1,stddev
      else
         if (lumsg.ne.0) write (lumsg,1009) stddev
      endif
!
!----------------------------------------------------------------------
! If specified, perform a natural cubic spline on the data in order 
! to interpolate the desired number of points to be used for fitting.
!----------------------------------------------------------------------
      if (nspln .gt. 0) then
         bi=xraw(1)
         delb=(xraw(npts)-xraw(1))/(nspln-1)
         call ncspln(xraw,yraw,npts,y2 )
         call splnay(xraw,yraw,y2,npts,bi,delb,spndat,nspln)
         if (lumsg.ne.0) write(lumsg,1010) nspln
      else
!
!     No spline: copy raw data directly
!
         bi=xraw(1)
         delb=(xraw(npts)-xraw(1))/(npts-1)
         nspln=npts
         do i=1,nspln
            spndat(i)=yraw(i)
         end do
         if (lumsg.ne.0) write (lumsg,1011) 
      end if
!
!----------------------------------------------------------------------
!     Normalize data to adjust integral (or double integral for
!     first-derivative data) to unity
!----------------------------------------------------------------------
      if (norm.ne.0) then
         call normlz( spndat,nspln,delb,ideriv,lumsg,anorm )
         if (anorm.ne.0.0d0) stddev=stddev/anorm
      end if
!
      getdat=0
      return
!
!### format statements ############################################
!
 1005 format(//'Opening file ',a/)
 1006 format(i6,' raw data points read in PC-GENPLOT binary format'/
     #       '    Field range:',f10.2,' to ',f10.2)
 1007 format(i6,' raw data points read in ASCII format'/
     #       '    Field range:',f10.2,' to ',f10.2/)
 1008 format('  Data baseline-corrected using linear fit to ',i3,
     #' points at each end'/'  Intercept=',g14.7,'  Slope=',g14.7,
     #' Noise: ',g14.7)
 1009 format('  Spectral noise: ',g14.7)
 1010 format('  Data splined to ',i4,' points')
 1011 format('  Input data not splined')
 1012 format(' *** Questionable field range ***')
 1013 format(' *** Maximum of ',i5,' input pts reached: remaining',i5,
     #     ' pts ignored ***')
      end


!-----------------------------------------------------------------------
!                       =========================
!                          function GETASC
!                       =========================
! Input (X,Y) data in ASCII format with one (X,Y) pair per line
! Each line is read in as a character variable; if the first character 
! of the line is a 'c' or a '!', it is considered to be a comment and
! the line is stored in the comment array. Otherwise, the (X,Y) pair
! is input from the line using a FORTRAN internal read.
!
!  Arguments:
!    filen:   Name of file to be read
!    x:       Double precision array to receive x data
!    y:       Double precision array to receive y data
!    comment: Character array to receive file comments
!    ncmt:    Number of comments read from file
!    maxcmt:  Maximum number of comments that may be returned
!    maxinp:  Maximum number of points that may be returned
!-----------------------------------------------------------------------
      function getasc(filen,x,y,n,comment,ncmt,maxcmt,maxinp)
!
      use stdio
!
      implicit none
      integer getasc,n,ncmt,maxcmt,maxinp,lth
      character numstr*30
!     double precision x(n),y(n),xin,yin // n is "output only"
      double precision x(*),y(*),xin,yin
!
      character*30 filen
      character*80 comment(maxcmt),line
!
      logical ftoken
      external ftoken
!
!######################################################################
!
      ncmt=0
      open(ludisk,file=filen,status='old', 
     #      access='sequential',form='formatted',err=99 )

    7 read (ludisk,1007,end=10,err=99) line
      if (  line(1:1) .eq. 'c' 
     # .or. line(1:1) .eq. 'C'
     # .or. line(1:1) .eq. '!' ) then
         if (ncmt .lt. maxcmt) then
            ncmt=ncmt+1
            comment(ncmt)=line
         endif
!
         goto 7
      end if
!
      call gettkn(line,numstr,lth)
      if (.not.ftoken(numstr,lth,x(1)) ) go to 99
      call gettkn(line,numstr,lth)
      if (.not.ftoken(numstr,lth,y(1)) ) go to 99
 1007 format(a)
!
!---------------------------------------------
! loop to read input lines, skipping comments
!---------------------------------------------
      n=1
 8    read (ludisk,1007,end=10,err=99) line
        if ( line(1:1) .eq. 'c' 
     #  .or. line(1:1) .eq. 'C'
     #  .or. line(1:1) .eq. '!' ) then
        if (ncmt .lt. maxcmt) then
          ncmt=ncmt+1
          comment(ncmt)=line
        endif
        goto 8
      endif
!
      call gettkn(line,numstr,lth)
      if (.not.ftoken(numstr,lth,xin) ) go to 99
      call gettkn(line,numstr,lth)
      if (.not.ftoken(numstr,lth,yin) ) go to 99
      n=n+1
      if (n.le.maxinp) then
         x(n)=xin
         y(n)=yin
      end if
      goto 8
!
!----------------------------------------------------------------------
!     Normal return
!----------------------------------------------------------------------
 10   close( ludisk )
      if (n.gt.0) getasc=0
      return

!----------------------------------------------------------------------
!     Error opening or reading file
!----------------------------------------------------------------------
 99   close( ludisk )
      n=0
      getasc=1
      return
      end


!----------------------------------------------------------------------
!                      =========================
!                          function SGLINT
!                      =========================
! 
! Calculate the integral of the function tabulated in ARRY using the
! trapezoidal rule. It is assumed that the ordinates are equally spaced.
! The return value is normalized to a spacing DX=1.
! N is number of points in ARRY.
!----------------------------------------------------------------------
      function sglint( arry, n )
      implicit none
      integer i,n
      double precision arry(n), sglint

      sglint=0.0D0
      do 10 i=1,n
        sglint=sglint+arry(i)
 10   continue
      sglint=sglint-0.5D0*(arry(1)+arry(n))
      return
      end

!----------------------------------------------------------------------
!                      =========================
!                          function DBLINT
!                      =========================
! 
! Calculate the double integral of the function tabulated in ARRY using 
! the trapezoidal rule. It is assumed that the ordinates are equally spaced.
! The return value is normalized to a spacing DX=1.
! N is number of points in ARRY.
!----------------------------------------------------------------------
      function dblint( arry, n )
      implicit none
      integer i,n
      double precision arry(n),sglint,dblint

      sglint=0.0D0
      dblint=0.0D0
      do 10 i=1,n
        sglint=sglint+arry(i)
        dblint=dblint+sglint
 10   continue
      dblint=dblint-0.5*(arry(1)+sglint)
      return
      end


!----------------------------------------------------------------------
!                       =======================
!                         subroutine LINEBC
!                       =======================
!
! Given a set of N points in arrays X and Y, fit a least-squares line to 
! the first and last NEND points of the spectrum and subtract it from the
! Y array. Returns slope and intercept of the baseline and the RMS deviation
! of the data from the baseline at the endpoints (a measure of spectral
! noise) 
!
!       Inputs
!          x      array of x-values
!          y      array of y-values (replaced on output)
!          n      number of points
!          modebc number of endpoints to use for baseline correction
!                 if modebc.eq.0, baseline is calculated but not
!
!       Outputs
!          c0     Zeroth order (intercept) coefficient for baseline
!          c1     First order (slope) coefficient for baseline
!          stddev Std deviation of data from baseline (spectral RMS noise)
!
!----------------------------------------------------------------------
      subroutine linebc(x,y,n,modebc,c0,c1,stddev)
      implicit none
      integer n,modebc
      double precision x(n),y(n),c0,c1,stddev
      double precision d,sn,sx,sx2,sxy,sy
      integer i,k,nend
!
      double precision ZERO,TWO
      parameter(ZERO=0.0D0,TWO=2.0D0)
!
      nend=min0(modebc,n/3)
      if (nend.le.0) nend=max0(10,n/20)
!
      sx =ZERO
      sy =ZERO
      sx2=ZERO
      sxy=ZERO
      do i=1,nend
         k=n-i+1
        sx =sx + x(i)+x(k)
        sy =sy + y(i)+y(k)
        sx2=sx2 + x(i)**2 + x(k)**2
        sxy=sxy + x(i)*y(i) + x(k)*y(k)
      end do
!
!     --------------------------------------------------------------
!     Calculate slope and intercept of baseline and spectral noise
!     --------------------------------------------------------------
!
      sn=2*nend
      d =sn*sx2 - sx*sx
      c0=(sx2*sy - sx*sxy) / d
      c1=(sn*sxy-sx*sy) / d
!
!     --------------------------------
!     Calculate variance of linear fit
!     --------------------------------
!
      stddev=ZERO
      do i=1,nend
         k=n-i+1
         stddev=stddev+((y(i)+y(k))-TWO*c0-c1*(x(i)+x(k)))**2
      end do
      stddev=dsqrt(stddev/dfloat(2*nend-1))
!
!     ------------------------------
!     Subtract baseline from data
!     ------------------------------
!
      if (modebc.gt.0) then
         do i=1,n
            y(i)=y(i)-(c0+c1*x(i))
         end do
      end if
      return
      end

!----------------------------------------------------------------------
!                       =======================
!                         subroutine NCSPLN
!                       =======================
! Perform a natural cubic spline fit. This is a modification of 
! subroutine SPLINE from Numerical Recipes. Given the tabulated function
! of N points in arrays X and Y, return an array Y2 of length N which 
! contains the second derivative of the interpolating function at the 
! corresponding X values. The major difference from the original SPLINE 
! routine is that the second derivatives at the both boundaries of 
! the function are assumed to be zero (hence "natural cubic spline").
!----------------------------------------------------------------------
      subroutine ncspln(x,y,n,y2)
      implicit none
      integer NMAX
      parameter(NMAX=4096)
      double precision ZERO,ONE,TWO,SIX
      parameter(ZERO=0.0d0,ONE=1.0d0,TWO=2.0D0,SIX=6.0d0)
!
      integer i,k,n
      double precision x(n),y(n),y2(n),u(NMAX),sig,p
!      
      y2(1)=ZERO
      u(1) =ZERO
      do 18 i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1) + TWO
        y2(i)=(sig-ONE)/p
        u(i)=(SIX*( (y(i+1)-y(i))/(x(i+1)-x(i))
     &               -(y(i)-y(i-1))/(x(i)-x(i-1))  )
     &              /(x(i+1)-x(i-1))
     &         - sig*u(i-1))/p
18    continue
      y2(n)=ZERO
      do 19 k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
19    continue
      return
      end

!----------------------------------------------------------------------
!                         ====================
!                          subroutine SPLNAY
!                         ====================
!
! Modification of Numerical Recipes spline interpolation routine (SPLINT)
! which interpolates an entire array at once. Given the original N 
! data points in XA,YA and the second derivative YA2 (calculated by SPLINE
! or NCSPLN) and the set of NS X values starting at X0 and spaced by DX, 
! this routine returns the interpolated function values in Y.
! The original binary search in subroutine SPLINT has been replaced, 
! since the interpolations start near the lower bound and are repeated 
! for successively higher, but closely spaced x values.
!----------------------------------------------------------------------
      subroutine splnay(xa,ya,y2a,n,x0,dx,y,ns)
      implicit none
      integer i,klo,khi,n,ns
      double precision xa(n),ya(n),y2a(n),y(ns),x,x0,dx,h,a,b
!
      if (x0 .lt. xa(1) .or. x0+dx*(ns-1) .gt. xa(n)) then
        print *, 'SPLNAY: Bad interpolated X range.'
        print *, '[execution paused, press enter to continue]'
        read (*,*)
      end if
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
        if (h .eq. 0.0D0) then
          print *, 'SPLNAY: Bad XA input.'
          print *, '[execution paused, press enter to continue]'
          read (*,*)
        end if
        a=(xa(khi)-x)/h
        b=(x-xa(klo))/h
        y(i)=a*ya(klo) + b*ya(khi) +
     *        ((a**3-a)*y2a(klo) + (b**3-b)*y2a(khi))*(h**2)/6.0D0
3     continue
      return
      end

!----------------------------------------------------------------------
!                    =========================
!                       subroutine NORMLZ
!                    =========================
!
!     Normalizes an array to unit integral, given the x-spacing dx
!
!     Inputs
!       arry    Double precision y array
!       n       Length of array
!       dx      Assumed x-value increment (equally spaced array assumed)
!       ideriv  Derivative flag (0=0th, otherwise 1st)
!       lumsg   Logical unit to receive output messages (0 for no msgs)
!
!     Outputs
!       arry    Normalized array
!
!     Uses
!       sglint  Single integration of a data array
!       dblint  Double integration of a data array
!
!----------------------------------------------------------------------
!
      subroutine normlz(arry,n,dx,ideriv,lumsg,anorm)
      implicit none
      integer n,ideriv,lumsg
      double precision arry(n),dx,anorm
!
      double precision base,sntgrl
      integer i
!
      double precision sglint,dblint
      external sglint,dblint
!
      if (ideriv.eq.0) then
!
!     ---------------------------------------------------------------------
!     Find normalization constant by single integration of 0th-deriv data
!     ---------------------------------------------------------------------
         anorm=dx*sglint(arry,n)
!
      else
!
!     -----------------------------------------------------------------
!      If first-derivative data are to be normalized, first subtract
!      a constant so that the integral of the curve is zero 
!     -----------------------------------------------------------------
         sntgrl=sglint(arry,n)
         base=sntgrl/dfloat(n)
         do i=1,n
            arry(i)=arry(i)-base
         end do
         if (lumsg.ne.0) write(lumsg,1008) base
!
!        ------------------------------------------------------------------- 
!        Find normalization constant by double integration of 1st-deriv data
!        ------------------------------------------------------------------- 
         anorm=dx*dx*dblint(arry,n)
      end if
!
      if (lumsg.ne.0) write (lumsg,1012) anorm
      if (anorm.ne.0.0D0) then
         do i=1,n
            arry(i)=arry(i)/anorm
         end do
         
      end if
      return
!
 1008 format('  Baseline adjusted by subtracting ',g14.7,
     #        ' before normalization')
 1012 format('  Data normalized: integral = ',g14.7)
      end
