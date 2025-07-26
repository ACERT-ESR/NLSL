!  Version 1.4 10/10/94
!*********************************************************************
!
!               MODIFIED BESSEL FUNCTIONS OF THE FIRST
!               KIND OF INTEGER ORDER AND REAL ARGUMENT
!               ---------------------------------------
!
!       This double precision function subroutine calculates the
!       value of the modified Bessel function of the first kind of
!       integer order and strictly real argument via a backward
!       recurrence scheme taken from "Numerical Recipes", W. H. Press,
!       B. P. Flannery, S. A. Teulosky and W. T. Vetterling, 1st ed.,
!       Cambridge Univ. Press, 1986.
!
!       The differs from the standalone version in the addition of
!       the parameter ierr. It is used to report non-convergence of
!       Bessel function evaluations instead of halting when such
!       errors occur.   
!       
!       written by DJS 10-SEP-87
!       Bug in small-argument Taylor series expansion fixed by DEB OCT-92
!
!
!       Includes:
!               rndoff.inc
!
!       Uses:
!               bessi0.f
!               bessi1.f
!
!*********************************************************************
!
      function bessi(n,z,ierr)
!
      use rnddbl
!
      implicit none
      integer n,ierr
      double precision bessi,z
!
      double precision bessi0,bessi1
      external bessi0,bessi1
!
      integer iacc
      double precision BIGNO,BIGNI,ONE
      parameter (iacc=40,BIGNO=1.0D10,BIGNI=1.0D-10,ONE=1.0D0)
!
      integer i,m,mmax
      double precision x,phase,twobyx,bi,bip,bim
!
      intrinsic abs
!
!#####################################################################
!
      ierr=0
!
!---------------------------------------------------------------------
!       get proper phase factor if argument is negative with the
!       following rules
!                                         n
!       I (z) = I (z)  and   I (-z) = (-1) I (z)
!        n       -n           n             n
!---------------------------------------------------------------------
!
      m=abs(n)
      x=abs(z)
!
      if ((z.lt.0.0D0).and.(mod(m,2).eq.1)) then
        phase=-ONE
      else
        phase=ONE
      end if
!
!---------------------------------------------------------------------
!       return proper values if argument is zero
!---------------------------------------------------------------------
!       
      if (x.lt.RNDOFF) then
        if (m.eq.0) then
          bessi=ONE
        else
          bessi=0.0D0
        end if
        return
      end if
!
!---------------------------------------------------------------------
!       call bessi0 if n=0, bessi1 if n=1, or go through
!       downward recurrence if n>1.
!---------------------------------------------------------------------
!
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
!
      return
      end

!*********************************************************************
!
!               MODIFIED BESSEL FUNCTION OF THE FIRST
!               KIND OF ORDER ZERO AND REAL ARGUMENT
!               -------------------------------------
!
!       This double precision function subroutine calculates the
!       modified Bessel function of the first kind of order zero
!       and real argument by either the Taylor series expansion
!       for small arguments or the first term of the asymptotic
!       series for sufficiently large arguments.
!
!       written by DJS 10-SEP-87
!
!       Includes:
!               rndoff.inc
!               pidef.inc
!
!       Uses:
!
!*********************************************************************
!
      double precision function bessi0(z,ierr)
!
      use rnddbl
      use pidef
!
      implicit none
      double precision z
      integer ierr
!
      integer i,j
      double precision x,y,smax,temp1,temp2,temp3,sum
!
      integer NMAX
      parameter (NMAX=40)
!
      double precision tser
      dimension tser(NMAX)
!
      double precision CUTOFF,ONE
      parameter (CUTOFF=20.0D0,ONE=1.0D0)
!
!######################################################################
!
      bessi0=0.0D0
      y=abs(z)
!
!------------------------------------------------------------
!     Set function value to unity if argument is too small
!------------------------------------------------------------
      if (y.lt.RNDOFF) then
        bessi0=ONE
!
!-------------------------------------------------------------
!     Taylor series expansion for small to moderate arguments
!-------------------------------------------------------------
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
!
        bessi0=0.0D0
        do j=i-1,1,-1
           bessi0=bessi0+tser(j)
        end do
        bessi0=bessi0+ONE
!
!----------------------------------------------
!     Asymptotic expansion for large arguments
!----------------------------------------------
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
!
      return
      end

!*********************************************************************
!
!               MODIFIED BESSEL FUNCTION OF THE FIRST
!               KIND OF ORDER ONE AND REAL ARGUMENT
!               -------------------------------------
!
!       This double precision function subroutine calculates the
!       modified Bessel function of the first kind of order one
!       and real argument by either the Taylor series expansion
!       for small arguments or the first term of the asymptotic
!       series for sufficiently large arguments.
!
!       written by DJS 10-SEP-87
!
!       Includes:
!               rndoff.inc
!               pidef.inc
!
!       Uses:
!
!*********************************************************************
!
      double precision function bessi1(z,ierr)
!
      use rnddbl
      use pidef
!
      implicit none
      double precision z
      integer ierr
!
      integer i,j
      double precision x,y,smax,temp1,temp2,temp3,phase,sum
!
      integer NMAX
      parameter (NMAX=40)
!
      double precision series
      dimension series(NMAX)
!
      double precision CUTOFF,ONE
      parameter (CUTOFF=20.0D0,ONE=1.0D0)
!
!#####################################################################
!
      bessi1=0.0D0
      if (z.gt.0.0D0) then
        phase=ONE
        y=z
      else
        phase=-ONE
        y=-z
      end if
!
!----------------------------------------------------------------------
!     set answer to zero if argument is too small, otherwise
!----------------------------------------------------------------------
      if (y.lt.RNDOFF) then
        bessi1=0.0D0
!
!----------------------------------------------------------------------
!     Use Taylor series expansion for small to moderate arguments or
!----------------------------------------------------------------------
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
!     
!----------------------------------------------------------------------
!     asymptotic expansion for large arguments
!----------------------------------------------------------------------
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
!
      return
      end
