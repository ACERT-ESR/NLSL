c     Version 1.5   5/2/94
c*********************************************************************
c
c               MODIFIED BESSEL FUNCTION OF THE FIRST
c               KIND OF ORDER ONE AND REAL ARGUMENT
c               -------------------------------------
c
c       This double precision function subroutine calculates the
c       modified Bessel function of the first kind of order one
c       and real argument by either the Taylor series expansion
c       for small arguments or the first term of the asymptotic
c       series for sufficiently large arguments.
c
c       written by DJS 10-SEP-87
c
c       Includes:
c               rndoff.inc
c               pidef.inc
c
c       Uses:
c
c*********************************************************************
c
      double precision function bessi1(z)
c
      include 'rndoff.inc'
      include 'pidef.inc'
c
      double precision z
c
      integer i,j
      double precision x,y,smax,temp1,temp2,temp3,phase,sum
c
      integer nmax
      parameter (nmax=40)
c
      double precision series
      dimension series(nmax)
c
      double precision cutoff
      parameter (cutoff=20.0D0)
c
c#####################################################################
c
      if (z.gt.0.0D0) then
        phase=1.0D0
        y=z
      else
        phase=-1.0D0
        y=-z
      end if
c
c----------------------------------------------------------------------
c     set answer to zero if argument is too small, otherwise
c----------------------------------------------------------------------
      if (y.lt.rndoff) then
        bessi1=0.0D0
c
c----------------------------------------------------------------------
c     Use Taylor series expansion for small to moderate arguments or
c----------------------------------------------------------------------
      else if (y.le.cutoff) then
        x=y*y*0.25D0
        temp1=1.0D0
        smax=1.0D0
        i=1
 10     temp1=(temp1/dble(i))*(x/dble(i+1))
          if (i.gt.nmax) then
            write(*,1000)
            stop
          end if
          series(i)=temp1
          i=i+1
          if (temp1.gt.smax) smax=temp1
          if (temp1/smax.gt.rndoff) go to 10
          bessi1=0.0D0
          do 20 j=i-1,1,-1
            bessi1=bessi1+series(j)
 20         continue
          bessi1=phase*y*0.5D0*(bessi1+1.0D0)
c
c----------------------------------------------------------------------
c     asymptotic expansion for large arguments
c----------------------------------------------------------------------
      else
        x=0.125D0/y
        sum=3.0D0*x
        temp3=sum
        smax=1.0D0
        i=2
 30     temp1=dble(2*i-1)
        temp1=temp1*temp1-4.0D0
        temp2=(x*temp1)/dble(i)
        if (temp2.gt.1.0D0) go to 40
        temp3=temp3*temp2
        if (temp3.gt.smax) smax=temp3
        if (temp3/smax.lt.rndoff) go to 40
        sum=sum+temp3
        i=i+1
        go to 30
 40     bessi1=dexp(y)*(1.0D0-sum)/dsqrt(y*(pi+pi))
      end if
c
      return
c
c----------------------------------------------------------------------
 1000 format('bessi0: Taylor series did not converge')
      end
