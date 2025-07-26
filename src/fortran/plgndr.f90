! Version 1.4 10/10/94
!*********************************************************************
!
!                   ASSOCIATED LEGENDRE FUNCTIONS
!                   -----------------------------
!
!       This double precision function subroutine calculates the
!       value of the associated Legendre function via an upward
!       recurrence scheme taken from "Numerical Recipes", W. H. Press,
!       B. P. Flannery, S. A. Teulosky and W. T. Vetterling,1st ed.,
!       Cambridge Univ. Press, 1986.
!
!       written by DJS 10-SEP-87
!
!       Includes:
!
!       Uses:
!
!*********************************************************************
!
      function plgndr(l,m,z)
      implicit none
!
      integer l,m
      double precision plgndr,z
!
      integer i
      double precision pmm,temp1,temp2,pmmp1,pmmp2
!
      intrinsic abs,sqrt
!
!#####################################################################
!
      if ((m.lt.0).or.(m.gt.l).or.(abs(z).gt.1.0D0)) then
        go to 9999
      end if
!
      pmm=1.0D0
!
!---------------------------------------------------------------------
!
!       calculate the starting point for the upward recurrence
!       on l using the formula :
!
!                m        m             2 m/2
!               P (z)=(-1) (2m-1)!! (1-z )
!                m
!
!---------------------------------------------------------------------
!
      if (m.gt.0) then
        temp1=sqrt((1.0D0-z)*(1.0D0+z))
        temp2=1.0D0
        do i=1,m
          pmm=-pmm*temp1*temp2
          temp2=temp2+2.0D0
        end do
      end if
!
!---------------------------------------------------------------------
!       
!       do upward recursion on l using the formulae :
!
!                     m              m              m
!               (l-m)P  (z)= z(2l-1)P (z) - (l+m-1)P (z)
!                     l              l-1            l-2
!
!       or
!                m               m
!               P (z) = z(2m+1) P (z)           if l=m+1
!                m+1             m
!
!---------------------------------------------------------------------
!
        if (l.eq.m) then
          plgndr=pmm
        else
          pmmp1=z*(2*m+1)*pmm
          if (l.eq.m+1) then
            plgndr=pmmp1
          else
            do i=m+2,l
              pmmp2=(z*(2*i-1)*pmmp1-(i+m-1)*pmm)/(i-m)
              pmm=pmmp1
              pmmp1=pmmp2
              plgndr=pmmp2
            end do
          end if
          end if
!
!---------------------------------------------------------------------
!     return to calling program
!---------------------------------------------------------------------
!
          return
!
!---------------------------------------------------------------------
!     exit from program if improper arguments detected
!---------------------------------------------------------------------
!
 9999     write (*,1000) l,m,z
 1000     format('Improper arguments for PLGNDR: ',
     #    i3,',',i3,',',g13.6)
          stop
          end
