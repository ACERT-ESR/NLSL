!  VERSION 1.3  (NLS version)   1/13/93
!**********************************************************************
!                         ====================
!                           subroutine ANXLK
!                         ====================
!
!       This subroutine calculates the xlk coefficients for
!       the potential-dependent part of the diffusion operator
!       for use in the SLE matrix calculation (MATRL, MATRF).
!       Array xlk (in /eprdat/) should be dimensioned to at least
!       (5,5), so that it can accommodate even L,K values from 0
!       to 8 (twice the highest possible L of a potential coefficient)
!
!       Notes:
!          The summation over the X(L,K)'s proceeds from 0 to 2*lptmx
!          X(L,K) is nonzero only for even L,K
!          X(L,K) = X(L,-K)  [similar to potential coefficients]
!
!          This code has been updated to include the possibility of a
!          nonaxial diffusion tensor (Rx .ne. Ry .ne. Rz). It was 
!          developed from subroutine CALXLK by D.J. Schneider which was 
!          based on the original subroutine matr written by G.Moro.
!
!
!       written by DEB 8-JUL-92
!
!       Includes:
!               rndoff.inc
!               eprprm.inc
!
!       Uses:
!               w3j.f
!
!**********************************************************************
!
      subroutine anxlk(rp,rm,rz)
!
      use rnddbl
      use eprprm
!
      implicit none
      double precision rp,rm,rz
!
      integer k,kx,kptf,k1f,k1i,k1,k1abs,k1x,k2,k2abs,k2f,k2i,k2x,
     #     l,lx,l1,l1x,l2,l2x
      double precision factl,wig1,rj2u,rjuju,term
!
      double precision w3j
      external w3j
!
!######################################################################
!
!
      do lx=1,5
        do kx=1,5
          xlk(lx,kx)=0.0D0
        end do
      end do
!
!----------------------------------------------------------------------
!     exit if no potential
!----------------------------------------------------------------------
!
      if(ipt.eq.0) return
!
!---------------------------------------------------------------------
!     calculate xlk coefficients
!----------------------------------------------------------------------
!
! --- loop over L
!
      do 100 l=0,lptmx*2,2
         lx=1+l/2
         factl=dble(2*l+1)
         kptf=min(l,2*kptmx+2)
!
! --- loop over K ---
!
         do 90 k=0,kptf,2
            kx=1+k/2
!     
!------------------------------
!         (J*R*J)U term
!------------------------------
!
            rj2u=cpot(lx,kx)*( rp*dble(l*(l+1)-k*k) + rz*dble(k*k) )
            if (k+2.le.l .and. k+2.le.kptmx) rj2u=rj2u+rm*cpot(lx,kx+1)
     #           * dsqrt(dble((l+k+1)*(l+k+2)*(l-k-1)*(l-k) ) )
            if (k-2.ge.0) rj2u=rj2u+rm*cpot(lx,kx-1)
     #           * dsqrt(dble((l-k+1)*(l-k+2)*(l+k-1)*(l+k) ) )
            if (k-2.lt.0) rj2u=rj2u+rm*cpot(lx,kx+1)
     #           * dsqrt(dble((l-k+1)*(l-k+2)*(l+k-1)*(l+k) ) )
            xlk(lx,kx)=-0.5d0*rj2u
!
!------------------------------
!      (JU)*R*(JU) term
!------------------------------
!
            rjuju = 0.0d0
!
!      --- loop over L1
!
            do 80 l1=0,lptmx,2
               l1x=1+l1/2
               k1f=min(l1,kptmx)
               k1i=-k1f
!
!        --- loop over L2
!
               do 70 l2=0,lptmx,2
                  l2x=1+l2/2
                  if(l1+l2.ge.l) then
                     wig1=w3j(l1,l,l2,0,0,0)
                  else
                     wig1=0.0D0
                     go to 70
                  end if
!
!         --- loop over K1
!
                  do 60 k1=k1i,k1f,2
                     k1abs=abs(k1)
                     k1x=1+k1abs/2
!
!        ---- loop over K2
!
                     k2i=max(k-k1-2,-kptmx,-l2)
                     k2f=min(k-k1+2,kptmx,l2)
!                
                     do 50 k2=k2i,k2f,2
                        k2abs=abs(k2)
                        k2x=1+k2abs/2
!                
!        ------- (J_+ U)(J_+ U) term
!                
                        if( (k2.eq.k-k1-2) .and. (abs(k1+1).le.l1)
     #                       .and.(abs(k2+1).le.l2) ) then
                          term=rm*w3j(l1,l,l2,k1+1,-k,k2+1)* dsqrt(
     #                    dble((l1-k1)*(l1+k1+1)*(l2-k2)*(l2+k2+1)) )
!
!       ------- (J_- U)(J_- U) term
!
                       else if ( (k2.eq.k-k1+2) .and. (abs(k1-1).le.l1)
     #                         .and. (abs(k2-1).le.l2) ) then
                         term=rm*w3j(l1,l,l2,k1-1,-k,k2-1)* dsqrt(
     #                   dble((l1+k1)*(l1-k1+1)*(l2+k2)*(l2-k2+1)) )
!                
!      ------- (J_+ U)(J_- U) term
!                
                else if (k2.eq.k-k1) then
                   if((abs(k1+1).le.l1).and.(abs(k2-1).le.l2)) then
                      term=rp*w3j(l1,l,l2,k1+1,-k,k2-1)*
     #                 dsqrt(dble((l1-k1)*(l1+k1+1)*(l2+k2)*(l2-k2+1)))
                   else
                      term=0.0D0
                   end if
!
!      ------ (Jz U)(Jz U) term
!
                   if (k2abs.le.l2) 
     #                  term=term+rz*dble(k1*k2)*w3j(l1,l,l2,k1,-k,k2)
!
                else
                   term=0.0d0
                end if
!
                rjuju=rjuju+cpot(l1x,k1x)*cpot(l2x,k2x)*wig1*term
!
 50          continue
 60       continue
 70    continue
 80   continue
      xlk(lx,kx)=xlk(lx,kx)-0.25D0*factl*rjuju
 90   continue
 100  continue
!
      do 120 lx=1,5
        do 110 kx=1,5
          if (abs(xlk(lx,kx)).lt.rndoff) xlk(lx,kx)=0.0D0
 110    continue
 120  continue
!
      return
      end
