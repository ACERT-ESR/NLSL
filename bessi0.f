c      Version 1.5  5/2/94
c*********************************************************************
c
c               MODIFIED BESSEL FUNCTION OF THE FIRST
c               KIND OF ORDER ZERO AND REAL ARGUMENT
c               -------------------------------------
c
c       This double precision function subroutine calculates the
c       modified Bessel function of the first kind of order zero
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
      double precision function bessi0(z)
c
      include 'rndoff.inc'
      include 'pidef.inc'
c
      double precision z
c
      integer i,j
      double precision x,y,smax,temp1,temp2,temp3,sum
c
      integer nmax
      parameter (nmax=40)
c
      double precision tser
      dimension tser(nmax)
c
      double precision cutoff
      parameter (cutoff=20.0D0)
c
c######################################################################
c
      y=abs(z)
c
c------------------------------------------------------------
c     Set function value to unity if argument is too small
c------------------------------------------------------------
      if (y.lt.rndoff) then
        bessi0=1.0D0
c
c-------------------------------------------------------------
c     Taylor series expansion for small to moderate arguments
c-------------------------------------------------------------
      else if (y.le.cutoff) then
        x=y*y*0.25D0
        temp1=1.0D0
        smax=1.0D0
        i=1
 10     temp1=(temp1/dble(i))*(x/dble(i))
          if (i.gt.nmax) then
            write(*,1000)
            stop
          end if
          tser(i)=temp1
          i=i+1
          if (temp1.gt.smax) smax=temp1
          if (temp1/smax.gt.rndoff) go to 10
c
        bessi0=0.0D0
        do 20 j=i-1,1,-1
          bessi0=bessi0+tser(j)
 20       continue
        bessi0=bessi0+1.0D0
c
c----------------------------------------------
c     Asymptotic expansion for large arguments
c----------------------------------------------
      else
        x=0.125D0/y
        sum=0.0D0
        temp3=1.0D0
        smax=1.0D0
        i=1
 30     temp1=dble(2*i-1)
          temp2=(x*temp1)*(temp1/dble(i))
          if (temp2.gt.1.0D0) go to 40
          temp3=temp3*temp2
          if (temp3.gt.smax) smax=temp3
          if (temp3/smax.lt.rndoff) go to 40
            sum=sum+temp3
            i=i+1
            go to 30
 40     bessi0=dexp(y)*((sum+1.0D0)/dsqrt(y*(pi+pi)))
      end if
c
      return
c
 1000 format('bessi0: Taylor series did not converge')
      end
