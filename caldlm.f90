c       VERSION 1.0     12/17/88
c**********************************************************************
c
c                  2
c       Calculate d   (beta) as in "Angular Momentum" by Brink and
c                  k,m
c       Satchler, p.24 ,second ed., Clarendon Press, Oxford (1979).
c       Beta is specified in degrees.
c
c                            2
c       note : dlm(i+3,j+3)=d   (psi)
c                            i,j
c
c       written by DJS 11-SEP-87
c
c       Includes:
c               rndoff.inc
c               pidef.inc
c
c       Uses:
c
c**********************************************************************
c
      subroutine caldlm(d2km,beta)
c
      use pidef
      use rnddbl
c
      double precision beta,d2km
      dimension d2km(5,5)
c
      double precision dsq32,dsq38,c,sn,sn2,cs,cs2
c
c######################################################################
c
      dsq32=dsqrt(3.0D0/2.0D0)
      dsq38=dsqrt(3.0D0/8.0D0)
c
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
c
      d2km(5,5)=0.25D0*(1.0D0+cs)*(1.0D0+cs)
      d2km(1,1)=d2km(5,5)
c
      d2km(5,4)=-0.5D0*sn*(1.0D0+cs)
      d2km(4,5)=-d2km(5,4)
      d2km(1,2)=-d2km(5,4)
      d2km(2,1)=d2km(5,4)
c
      d2km(5,3)=dsq38*sn2
      d2km(3,5)=d2km(5,3)
      d2km(1,3)=d2km(5,3)
      d2km(3,1)=d2km(5,3)
c
      d2km(5,2)=0.5D0*sn*(cs-1.0D0)
      d2km(4,1)=d2km(5,2)
      d2km(1,4)=-d2km(5,2)
      d2km(2,5)=-d2km(5,2)
c
      d2km(5,1)=(0.5d0*(1.0D0-cs))**2
      d2km(1,5)=d2km(5,1)
c
      d2km(4,4)=0.5D0*(2.0D0*cs-1.0D0)*(cs+1.0D0)
      d2km(2,2)=d2km(4,4)
c
      d2km(4,2)=0.5d0*(2.0D0*cs+1.0D0)*(1.0D0-cs)
      d2km(2,4)=d2km(4,2)
c
      d2km(4,3)=-dsq32*sn*cs
      d2km(3,2)=d2km(4,3)
      d2km(3,4)=-d2km(4,3)
      d2km(2,3)=-d2km(4,3)
c
      d2km(3,3)=0.5D0*(3.0D0*cs2-1.0D0)
c
c----------------------------------------------------------------------
c     return to calling program
c----------------------------------------------------------------------
c
      return
      end
