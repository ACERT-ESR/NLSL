c    Version 1.5  5/2/94
c*********************************************************************
c
c                   ASSOCIATED LEGENDRE FUNCTIONS
c                   -----------------------------
c
c       This double precision function subroutine calculates the
c       value of the associated Legendre function via an upward
c       recurrence scheme taken from "Numerical Recipes", W. H. Press,
c       B. P. Flannery, S. A. Teulosky and W. T. Vetterling,1st ed.,
c       Cambridge Univ. Press, 1986.
c
c       written by DJS 10-SEP-87
c
c       Includes:
c
c       Uses:
c
c*********************************************************************
c
      double precision function plgndr(l,m,z)
c
      integer l,m
      double precision z
c
      integer i
      double precision pmm,temp1,temp2,pmmp1,pmmp2
c
      intrinsic abs,sqrt
c
c#####################################################################
c
      if ((m.lt.0).or.(m.gt.l).or.(abs(z).gt.1.0D0)) then
        go to 9999
      end if
c
      pmm=1.0D0
c
c---------------------------------------------------------------------
c
c       calculate the starting point for the upward recurrence
c       on l using the formula :
c
c                m        m             2 m/2
c               P (z)=(-1) (2m-1)!! (1-z )
c                m
c
c---------------------------------------------------------------------
c
      if (m.gt.0) then
        temp1=sqrt((1.0D0-z)*(1.0D0+z))
        temp2=1.0D0
        do 10 i=1,m
          pmm=-pmm*temp1*temp2
 10       temp2=temp2+2.0D0
        end if
c
c---------------------------------------------------------------------
c       
c       do upward recursion on l using the formulae :
c
c                     m              m              m
c               (l-m)P  (z)= z(2l-1)P (z) - (l+m-1)P (z)
c                     l              l-1            l-2
c
c       or
c                m               m
c               P (z) = z(2m+1) P (z)           if l=m+1
c                m+1             m
c
c---------------------------------------------------------------------
c
        if (l.eq.m) then
          plgndr=pmm
        else
          pmmp1=z*(2*m+1)*pmm
          if (l.eq.m+1) then
            plgndr=pmmp1
          else
            do 20 i=m+2,l
              pmmp2=(z*(2*i-1)*pmmp1-(i+m-1)*pmm)/(i-m)
              pmm=pmmp1
 20           pmmp1=pmmp2
              plgndr=pmmp2
            end if
          end if
c
c---------------------------------------------------------------------
c     return to calling program
c---------------------------------------------------------------------
c
          return
c
c---------------------------------------------------------------------
c     exit from program if improper arguments detected
c---------------------------------------------------------------------
c
 9999     stop
          end
