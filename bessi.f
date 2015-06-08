c     Version 1.5  5/2/94
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
c       written by DJS 10-SEP-87
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
      double precision function bessi(n,z)
c
      include 'rndoff.inc'
c
      integer n
      double precision z
c
      double precision bessi0,bessi1
      external bessi0,bessi1
c
      integer iacc
      double precision bigno,bigni
      parameter (iacc=40,bigno=1.0D10,bigni=1.0D-10)
c
      integer i,m,mmax
      double precision x,phase,twobyx,bi,bip,bim
c
c#####################################################################
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
        phase=-1.0D0
      else
        phase = 1.0D0
      end if
c
c---------------------------------------------------------------------
c       return proper values if argument is zero
c---------------------------------------------------------------------
c       
      if (x.lt.rndoff) then
        if (m.eq.0) then
          bessi=1.0D0
        else
          bessi=0.0D0
        end if
        go to 9999
      end if
c
c---------------------------------------------------------------------
c       call bessi0 if n=0, bessi1 if n=1, or go through
c       downward recurrence if n>1.
c---------------------------------------------------------------------
c
      if (m.eq.0) then
        bessi=phase*bessi0(x)
      else if (m.eq.1) then
        bessi=phase*bessi1(x)
      else
        bessi=0.0D0
        twobyx=2.0D0/x
        bip=0.0D0
        bi=1.0D0
        mmax=2*((m+int(sqrt(dble(iacc*m)))))
        do 10 i=mmax,1,-1
          bim=bip+dble(i)*twobyx*bi
          bip=bi
          bi=bim
          if (abs(bi).gt.bigno) then
            bessi=bessi*bigni
            bi=bi*bigni
            bip=bip*bigni
          end if
          if (i.eq.m) bessi=bip
 10     continue
        bessi=phase*bessi*bessi0(x)/bi
      end if
c
c---------------------------------------------------------------------
c       return to calling program
c---------------------------------------------------------------------
c
 9999 return
      end
