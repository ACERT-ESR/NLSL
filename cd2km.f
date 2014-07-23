c       VERSION 1.0     2/5/99
c**********************************************************************
c                    =========================
c                      subroutine CD2KM
c                    =========================
c
c                  2
c       Calculate D   (alpha,beta,gamma) as in "Angular Momentum" 
c                  k,m
c       by Brink and Satchler, p.24, 2nd ed., Clarendon Press, 
c       Oxford (1979).
c
c                                     2
c       note : d2km(1,i+3,j+3)=real{ D   (alpha,beta,gamma) }
c                                     i,j
c
c                                     2
c              d2km(2,i+3,j+3)=imag{ D   (alpha,beta,gamma) }
c                                     i,j
c
c       written by DEB using original routines of DJS and GM
c
c       Includes:
c               rndoff.inc
c               pidef.inc
c
c       Uses:
c
c**********************************************************************
c
      subroutine cd2km(d2km,alpha,beta,gamma)
c
      include 'rndoff.inc'
c
      double precision alpha,beta,gamma,d2km
      dimension d2km(2,5,5)
c
      double precision dsq32,dsq38,d,cb,sb,cb2,sb2,cd,sd
      integer i,j
c
c######################################################################
c
      dsq32=dsqrt(3.0D0/2.0D0)
      dsq38=dsqrt(3.0D0/8.0D0)
c
      call setcs(beta,cb,sb )
      cb2=cb*cb
      sb2=sb*sb
c
c----------------------------------------------------------------------
c     Set real parts of D2KM elements first
c----------------------------------------------------------------------
      d2km(1,5,5)=0.25D0*(1.0D0+cb)*(1.0D0+cb)
      d2km(1,1,1)=d2km(1,5,5)
c
      d2km(1,5,4)=-0.5D0*sb*(1.0D0+cb)
      d2km(1,4,5)=-d2km(1,5,4)
      d2km(1,1,2)=-d2km(1,5,4)
      d2km(1,2,1)=d2km(1,5,4)
c
      d2km(1,5,3)=dsq38*sb2
      d2km(1,3,5)=d2km(1,5,3)
      d2km(1,1,3)=d2km(1,5,3)
      d2km(1,3,1)=d2km(1,5,3)
c
      d2km(1,5,2)=0.5D0*sb*(cb-1.0D0)
      d2km(1,4,1)=d2km(1,5,2)
      d2km(1,1,4)=-d2km(1,5,2)
      d2km(1,2,5)=-d2km(1,5,2)
c
      d2km(1,5,1)=(0.5d0*(1.0D0-cb))**2
      d2km(1,1,5)=d2km(1,5,1)
c
      d2km(1,4,4)=0.5D0*(2.0D0*cb-1.0D0)*(cb+1.0D0)
      d2km(1,2,2)=d2km(1,4,4)
c
      d2km(1,4,2)=0.5d0*(2.0D0*cb+1.0D0)*(1.0D0-cb)
      d2km(1,2,4)=d2km(1,4,2)
c
      d2km(1,4,3)=-dsq32*sb*cb
      d2km(1,3,2)=d2km(1,4,3)
      d2km(1,3,4)=-d2km(1,4,3)
      d2km(1,2,3)=-d2km(1,4,3)
c
      d2km(1,3,3)=0.5D0*(3.0D0*cb2-1.0D0)
c
c----------------------------------------------------------------------
c    Set imaginary part of D2KM elements
c----------------------------------------------------------------------
      if (dabs(alpha).gt.rndoff .or. dabs(gamma).gt.rndoff) then
        do 20 i=1,5
          do 19 j=1,5
            d=(i-3)*alpha+(j-3)*gamma
            call setcs(d,cd,sd)
            d2km(2,i,j)=-d2km(1,i,j)*sd
            d2km(1,i,j)=d2km(1,i,j)*cd
   19       continue
   20     continue
      else
        do 30 i=1,5
          do 29 j=1,5
            d2km(2,i,j)=0.0D0
   29       continue
   30     continue
      end if
c
      return
      end


c----------------------------------------------------------------------
c                    =========================
c                       subroutine SETCS
c                    =========================
c   Sets sine and cosine of an angle (given in degrees) after checking
c   special cases
c----------------------------------------------------------------------
      subroutine setcs( ang, ca, sa )
      double precision ang,c,ca,sa
c
      include 'pidef.inc'
      include 'rndoff.inc'
c
c Change suggested by Liang made to nino's version 4/29/97:
c
      if (dabs(ang).lt.rndoff) then
         sa=0.0D0
         ca=1.0D0
      else if (dabs(180.0D0-ang).lt.rndoff) then
         sa=0.0D0
         ca=-1.0D0
      else if (dabs(90.0D0-ang).lt.rndoff) then
         sa=1.0D0
         ca=0.0D0
      else
        c=ang*pi/180.0D0
        sa=dsin(c)
        ca=dcos(c)
      end if
      return
      end
