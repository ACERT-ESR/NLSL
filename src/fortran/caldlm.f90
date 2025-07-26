!       VERSION 1.0     12/17/88
!**********************************************************************
!
!                  2
!       Calculate d   (beta) as in "Angular Momentum" by Brink and
!                  k,m
!       Satchler, p.24 ,second ed., Clarendon Press, Oxford (1979).
!       Beta is specified in degrees.
!
!                            2
!       note : dlm(i+3,j+3)=d   (psi)
!                            i,j
!
!       written by DJS 11-SEP-87
!
!       Includes:
!               rndoff.inc
!               pidef.inc
!
!       Uses:
!
!**********************************************************************
!
      subroutine caldlm(d2km,beta)
!
      use pidef
      use rnddbl
!
      double precision beta,d2km
      dimension d2km(5,5)
!
      double precision dsq32,dsq38,c,sn,sn2,cs,cs2
!
!######################################################################
!
      dsq32=dsqrt(3.0D0/2.0D0)
      dsq38=dsqrt(3.0D0/8.0D0)
!
      if (dabs(beta).lt.rndoff) then
        sn=0.0D0
        sn2=0.0D0
        cs=1.0D0
        cs2=1.0D0
      else if (dabs(90.0D0-beta).lt.rndoff) then
        sn=1.0D0
        sn2=1.0D0
        cs=0.0D0
        cs2=0.0D0
      else
        c=beta*pi/180.0D0
        sn=dsin(c)
        sn2=sn*sn
        cs=dcos(c)
        cs2=cs*cs
      end if
!
      d2km(5,5)=0.25D0*(1.0D0+cs)*(1.0D0+cs)
      d2km(1,1)=d2km(5,5)
!
      d2km(5,4)=-0.5D0*sn*(1.0D0+cs)
      d2km(4,5)=-d2km(5,4)
      d2km(1,2)=-d2km(5,4)
      d2km(2,1)=d2km(5,4)
!
      d2km(5,3)=dsq38*sn2
      d2km(3,5)=d2km(5,3)
      d2km(1,3)=d2km(5,3)
      d2km(3,1)=d2km(5,3)
!
      d2km(5,2)=0.5D0*sn*(cs-1.0D0)
      d2km(4,1)=d2km(5,2)
      d2km(1,4)=-d2km(5,2)
      d2km(2,5)=-d2km(5,2)
!
      d2km(5,1)=(0.5d0*(1.0D0-cs))**2
      d2km(1,5)=d2km(5,1)
!
      d2km(4,4)=0.5D0*(2.0D0*cs-1.0D0)*(cs+1.0D0)
      d2km(2,2)=d2km(4,4)
!
      d2km(4,2)=0.5d0*(2.0D0*cs+1.0D0)*(1.0D0-cs)
      d2km(2,4)=d2km(4,2)
!
      d2km(4,3)=-dsq32*sn*cs
      d2km(3,2)=d2km(4,3)
      d2km(3,4)=-d2km(4,3)
      d2km(2,3)=-d2km(4,3)
!
      d2km(3,3)=0.5D0*(3.0D0*cs2-1.0D0)
!
!----------------------------------------------------------------------
!     return to calling program
!----------------------------------------------------------------------
!
      return
      end
