c  VERSION 1.0  (NLSPMC version)   2/5/99
c**********************************************************************
c                         ====================
c                           subroutine ANXLK
c                         ====================
c
c       This subroutine calculates the xlk coefficients for
c       the potential-dependent part of the diffusion operator
c       for use in the SLE matrix calculation (MATRL, MATRF).
c       Array xlk (in /eprdat/) should be dimensioned to at least
c       (5,5), so that it can accommodate even L,K values from 0
c       to 8 (twice the highest possible L of a potential coefficient)
c
c       Notes:
c          The summation over the X(L,K)'s proceeds from 0 to 2*lptmx
c          X(L,K) is nonzero only for even L,K
c          X(L,K) = X(L,-K)  [similar to potential coefficients]
c
c          This code has been updated to include the possibility of a
c          nonaxial diffusion tensor (Rx .ne. Ry .ne. Rz). It was 
c          developed from subroutine CALXLK by D.J. Schneider which was 
c          based on the original subroutine matr written by G.Moro.
c
c
c       written by DEB 8-JUL-92
c
c       Includes:
c               nlsdim.inc
c               eprprm.inc
c               rndoff.inc
c
c       Uses:
c               w3j.f
c
c**********************************************************************
c
      subroutine anxlk(rp,rm,rz)
      implicit none
c
      include 'limits.inc'
      include 'simparm.inc'
      include 'rndoff.inc'
c
      double precision rp,rm,rz
c
      integer k,kx,kptf,k1f,k1i,k1,k1abs,k1x,k2,k2abs,k2f,k2i,k2x,
     #     l,lx,l1,l1x,l2,l2x
      double precision factl,wig1,rj2u,rjuju,term
c
      double precision w3j
      external w3j
c
c######################################################################
c
c
      do 20 lx=1,5
        do 10 kx=1,5
          xlk(lx,kx)=0.0D0
  10    continue
  20  continue
c
c----------------------------------------------------------------------
c     exit if no potential
c----------------------------------------------------------------------
c
      if(ipt.eq.0) return
c
c---------------------------------------------------------------------
c     calculate xlk coefficients
c----------------------------------------------------------------------
c
c --- loop over L
c
      do 100 l=0,lptmx*2,2
         lx=1+l/2
         factl=dble(2*l+1)
         kptf=min(l,2*kptmx+2)
c
c --- loop over K ---
c
         do 90 k=0,kptf,2
            kx=1+k/2
c     
c------------------------------
c         (J*R*J)U term
c------------------------------
c
            rj2u=cpot(lx,kx)*( rp*dble(l*(l+1)-k*k) + rz*dble(k*k) )
            if (k+2.le.l .and. k+2.le.kptmx) rj2u=rj2u+rm*cpot(lx,kx+1)
     #           * dsqrt(dble((l+k+1)*(l+k+2)*(l-k-1)*(l-k) ) )
            if (k-2.ge.0) rj2u=rj2u+rm*cpot(lx,kx-1)
     #           * dsqrt(dble((l-k+1)*(l-k+2)*(l+k-1)*(l+k) ) )
            if (k-2.lt.0) rj2u=rj2u+rm*cpot(lx,kx+1)
     #           * dsqrt(dble((l-k+1)*(l-k+2)*(l+k-1)*(l+k) ) )
            xlk(lx,kx)=-0.5d0*rj2u
c
c------------------------------
c      (JU)*R*(JU) term
c------------------------------
c
            rjuju = 0.0d0
c
c      --- loop over L1
c
            do 80 l1=0,lptmx,2
               l1x=1+l1/2
               k1f=min(l1,kptmx)
               k1i=-k1f
c
c        --- loop over L2
c
               do 70 l2=0,lptmx,2
                  l2x=1+l2/2
                  if(l1+l2.ge.l) then
                     wig1=w3j(l1,l,l2,0,0,0)
                  else
                     wig1=0.0D0
                     go to 70
                  end if
c
c         --- loop over K1
c
                  do 60 k1=k1i,k1f,2
                     k1abs=abs(k1)
                     k1x=1+k1abs/2
c
c        ---- loop over K2
c
                     k2i=max(k-k1-2,-kptmx,-l2)
                     k2f=min(k-k1+2,kptmx,l2)
c                
                     do 50 k2=k2i,k2f,2
                        k2abs=abs(k2)
                        k2x=1+k2abs/2
c                
c        ------- (J_+ U)(J_+ U) term
c                
                        if( (k2.eq.k-k1-2) .and. (abs(k1+1).le.l1)
     #                       .and.(abs(k2+1).le.l2) ) then
                          term=rm*w3j(l1,l,l2,k1+1,-k,k2+1)* dsqrt(
     #                    dble((l1-k1)*(l1+k1+1)*(l2-k2)*(l2+k2+1)) )
c
c       ------- (J_- U)(J_- U) term
c
                       else if ( (k2.eq.k-k1+2) .and. (abs(k1-1).le.l1)
     #                         .and. (abs(k2-1).le.l2) ) then
                         term=rm*w3j(l1,l,l2,k1-1,-k,k2-1)* dsqrt(
     #                   dble((l1+k1)*(l1-k1+1)*(l2+k2)*(l2-k2+1)) )
c                
c      ------- (J_+ U)(J_- U) term
c                
                else if (k2.eq.k-k1) then
                   if((abs(k1+1).le.l1).and.(abs(k2-1).le.l2)) then
                      term=rp*w3j(l1,l,l2,k1+1,-k,k2-1)*
     #                 dsqrt(dble((l1-k1)*(l1+k1+1)*(l2+k2)*(l2-k2+1)))
                   else
                      term=0.0D0
                   end if
c
c      ------ (Jz U)(Jz U) term
c
                   if (k2abs.le.l2) 
     #                  term=term+rz*dble(k1*k2)*w3j(l1,l,l2,k1,-k,k2)
c
                else
                   term=0.0d0
                end if
c
                rjuju=rjuju+cpot(l1x,k1x)*cpot(l2x,k2x)*wig1*term
c
 50          continue
 60       continue
 70    continue
 80   continue
      xlk(lx,kx)=xlk(lx,kx)-0.25D0*factl*rjuju
 90   continue
 100  continue
c
      do 120 lx=1,5
        do 110 kx=1,5
          if (abs(xlk(lx,kx)).lt.rndoff) xlk(lx,kx)=0.0D0
 110    continue
 120  continue
c
      return
      end
