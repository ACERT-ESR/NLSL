c  Version 1.4 10/10/94
c*********************************************************************
c
c               MODIFIED BESSEL FUNCTIONS OF THE FIRST
c               KIND OF INTEGER ORDER AND REAL ARGUMENT
c               ---------------------------------------
c
c       This double precision function subroutine calculates the
c       value of the modified Bessel function of the first kind of
c       integer order and strictly real argument via a backward
c       recurrence scheme taken from "Numerical Recipes", W. H. Press,
c       B. P. Flannery, S. A. Teulosky and W. T. Vetterling, 1st ed.,
c       Cambridge Univ. Press, 1986.
c
c       The differs from the standalone version in the addition of
c       the parameter ierr. It is used to report non-convergence of
c       Bessel function evaluations instead of halting when such
c       errors occur.   
c       
c       written by DJS 10-SEP-87
c       Bug in small-argument Taylor series expansion fixed by DEB OCT-92
c
c
c       Includes:
c               rndoff.inc
c
c       Uses:
c               bessi0.f
c               bessi1.f
c
c*********************************************************************
c
      function bessi(n,z,ierr)
c
      use rnddbl
c
      implicit none
      integer n,ierr
      double precision bessi,z
c
      double precision bessi0,bessi1
      external bessi0,bessi1
c
      integer iacc
      double precision BIGNO,BIGNI,ONE
      parameter (iacc=40,BIGNO=1.0D10,BIGNI=1.0D-10,ONE=1.0D0)
c
      integer i,m,mmax
      double precision x,phase,twobyx,bi,bip,bim
c
      intrinsic abs
c
c#####################################################################
c
      ierr=0
c
c---------------------------------------------------------------------
c       get proper phase factor if argument is negative with the
c       following rules
c                                         n
c       I (z) = I (z)  and   I (-z) = (-1) I (z)
c        n       -n           n             n
c---------------------------------------------------------------------
c
      m=abs(n)
      x=abs(z)
c
      if ((z.lt.0.0D0).and.(mod(m,2).eq.1)) then
        phase=-ONE
      else
        phase=ONE
      end if
c
c---------------------------------------------------------------------
c       return proper values if argument is zero
c---------------------------------------------------------------------
c       
      if (x.lt.RNDOFF) then
        if (m.eq.0) then
          bessi=ONE
        else
          bessi=0.0D0
        end if
        return
      end if
c
c---------------------------------------------------------------------
c       call bessi0 if n=0, bessi1 if n=1, or go through
c       downward recurrence if n>1.
c---------------------------------------------------------------------
c
      if (m.eq.0) then
        bessi=phase*bessi0(x,ierr)
        if (ierr.ne.0) return
      else if (m.eq.1) then
        bessi=phase*bessi1(x,ierr)
        if (ierr.ne.0) return
      else
        bessi=0.0D0
        twobyx=2.0D0/x
        bip=0.0D0
        bi=ONE
        mmax=2*((m+int(sqrt(dble(iacc*m)))))
        do i=mmax,1,-1
          bim=bip+dble(i)*twobyx*bi
          bip=bi
          bi=bim
          if (abs(bi).gt.BIGNO) then
            bessi=bessi*BIGNI
            bi=bi*BIGNI
            bip=bip*BIGNI
          end if
          if (i.eq.m) bessi=bip
        end do
        bessi=phase*bessi*bessi0(x,ierr)/bi
        if (ierr.ne.0) return
      end if
c
      return
      end

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
      double precision function bessi0(z,ierr)
c
      use rnddbl
      use pidef
c
      implicit none
      double precision z
      integer ierr
c
      integer i,j
      double precision x,y,smax,temp1,temp2,temp3,sum
c
      integer NMAX
      parameter (NMAX=40)
c
      double precision tser
      dimension tser(NMAX)
c
      double precision CUTOFF,ONE
      parameter (CUTOFF=20.0D0,ONE=1.0D0)
c
c######################################################################
c
      bessi0=0.0D0
      y=abs(z)
c
c------------------------------------------------------------
c     Set function value to unity if argument is too small
c------------------------------------------------------------
      if (y.lt.RNDOFF) then
        bessi0=ONE
c
c-------------------------------------------------------------
c     Taylor series expansion for small to moderate arguments
c-------------------------------------------------------------
      else if (y.le.CUTOFF) then
        x=y*y*0.25D0
        temp1=ONE
        smax=ONE
        i=1
 10     temp1=(temp1/dble(i))*(x/dble(i))
          if (i.gt.NMAX) then
            ierr=-1
            return
          end if
          tser(i)=temp1
          i=i+1
          if (temp1.gt.smax) smax=temp1
          if (temp1/smax.gt.RNDOFF) go to 10
c
        bessi0=0.0D0
        do j=i-1,1,-1
           bessi0=bessi0+tser(j)
        end do
        bessi0=bessi0+ONE
c
c----------------------------------------------
c     Asymptotic expansion for large arguments
c----------------------------------------------
      else
        x=0.125D0/y
        sum=0.0D0
        temp3=ONE
        smax=ONE
        i=1
 30     temp1=dble(2*i-1)
          temp2=(x*temp1)*(temp1/dble(i))
          if (temp2.gt.ONE) go to 40
          temp3=temp3*temp2
          if (temp3.gt.smax) smax=temp3
          if (temp3/smax.lt.RNDOFF) go to 40
            sum=sum+temp3
            i=i+1
            go to 30
 40     bessi0=dexp(y)*((sum+ONE)/dsqrt(y*(PI+PI)))
      end if
c
      return
      end

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
      double precision function bessi1(z,ierr)
c
      use rnddbl
      use pidef
c
      implicit none
      double precision z
      integer ierr
c
      integer i,j
      double precision x,y,smax,temp1,temp2,temp3,phase,sum
c
      integer NMAX
      parameter (NMAX=40)
c
      double precision series
      dimension series(NMAX)
c
      double precision CUTOFF,ONE
      parameter (CUTOFF=20.0D0,ONE=1.0D0)
c
c#####################################################################
c
      bessi1=0.0D0
      if (z.gt.0.0D0) then
        phase=ONE
        y=z
      else
        phase=-ONE
        y=-z
      end if
c
c----------------------------------------------------------------------
c     set answer to zero if argument is too small, otherwise
c----------------------------------------------------------------------
      if (y.lt.RNDOFF) then
        bessi1=0.0D0
c
c----------------------------------------------------------------------
c     Use Taylor series expansion for small to moderate arguments or
c----------------------------------------------------------------------
      else if (y.le.CUTOFF) then
         x=y*y*0.25D0
         temp1=ONE
         smax=ONE
         i=1
 10      temp1=(temp1/dble(i))*(x/dble(i+1))
         if (i.gt.NMAX) then
            ierr=-1
            return
         end if
         series(i)=temp1
         i=i+1
         if (temp1.gt.smax) smax=temp1
         if (temp1/smax.gt.RNDOFF) go to 10
         bessi1=0.0D0
         do j=i-1,1,-1
            bessi1=bessi1+series(j)
         end do
         bessi1=phase*y*0.5D0*(bessi1+ONE)
c     
c----------------------------------------------------------------------
c     asymptotic expansion for large arguments
c----------------------------------------------------------------------
      else
         x=0.125D0/y
         sum=3.0D0*x
         temp3=sum
         smax=ONE
         i=2
 30      temp1=dble(2*i-1)
         temp1=temp1*temp1-4.0D0
         temp2=(x*temp1)/dble(i)
         if (temp2.gt.ONE) go to 40
         temp3=temp3*temp2
         if (temp3.gt.smax) smax=temp3
         if (temp3/smax.lt.RNDOFF) go to 40
         sum=sum+temp3
         i=i+1
         go to 30
 40      bessi1=dexp(y)*(ONE-sum)/dsqrt(y*(PI+PI))
      end if
c
      return
      end
