c Version 1.5 5/2/94
c----------------------------------------------------------------------
c                    ======================
c                      subroutine FFT
c                    ======================
c
c     Discrete Fast Fourier Transform for a complex-valued array.
c     Adapted from Numerical Recipes by Press et al.
c----------------------------------------------------------------------
      subroutine fft(data,nn,isign)
      implicit none
      integer i,istep,j,m,mmax,n,nn,isign
      double precision wr,wi,wpr,wpi,wtemp,tempi,tempr,theta
      double precision data(2*nn)
c
      n=2*nn
      j=1
      do i=1,n,2
         if(j.gt.i)then
            tempr=data(j)
            tempi=data(j+1)
            data(j)=data(i)
            data(j+1)=data(i+1)
            data(i)=tempr
            data(i+1)=tempi
         endif
         m=n/2
1        if ((m.ge.2).and.(j.gt.m)) then
            j=j-m
            m=m/2
            go to 1
         endif
         j=j+m
      end do
c
      mmax=2
2     if (n.gt.mmax) then
         istep=2*mmax
         theta=6.28318530717959d0/(isign*mmax)
         wpr=-2.d0*dsin(0.5d0*theta)**2
         wpi=dsin(theta)
         wr=1.d0
         wi=0.d0
         do m=1,mmax,2
            do i=m,n,istep
               j=i+mmax
               tempr=sngl(wr)*data(j)-sngl(wi)*data(j+1)
               tempi=sngl(wr)*data(j+1)+sngl(wi)*data(j)
               data(j)=data(i)-tempr
               data(j+1)=data(i+1)-tempi
               data(i)=data(i)+tempr
               data(i+1)=data(i+1)+tempi
            end do
            wtemp=wr
            wr=wr*wpr-wi*wpi+wr
            wi=wi*wpr+wtemp*wpi+wi
         end do
         mmax=istep
      go to 2
      endif
      return
      end
