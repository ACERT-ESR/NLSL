c     Version 1.5  5/2/94
c----------------------------------------------------------------------
c                    =========================
c                       subroutine GCONVL
c                    =========================
c
c     Convolute given spectrum with a Gaussian lineshape having
c     a specified derivative peak-to-peak linewidth
c
c----------------------------------------------------------------------
      subroutine gconvl( spectr,wline,dfld,nfld )
      implicit none
c
      integer nfld
      double precision spectr(2,nfld),wline,dfld
c
      integer i,k,nft,no2
      double precision df,f,g,gib,gnorm
c
      include 'stddim.inc'
      include 'eprdat.inc'
      include 'pidef.inc'
c
      double precision tmpdat
      common /temp/ tmpdat(2,MXPT)
c
      double precision EIGHT,EPS,ONE,TWO,THIRD,ZERO
      parameter ( EIGHT=8.0D0,ONE=1.0D0,TWO=2.0D0,ZERO=0.0D0,
     #            EPS=1.0D-6,THIRD=0.3333333333333333D0 )
c
c######################################################################
c
      if (wline.lt.EPS) return
c
c      ----------------------------
c      Find next largest power of 2
c      ----------------------------
       nft=2
    1  nft=nft+nft
          if (nft.lt.nfld) goto 1    
c
c      -----------------------------------------------
c       Store calculation in tmpdat, zero-pad, and FFT
c      -----------------------------------------------
      do i=1,nfld
         tmpdat(1,i)=spectr(1,i)
         tmpdat(2,i)=spectr(2,i)
      end do
      do i=nfld+1,nft
         tmpdat(1,i)=ZERO
         tmpdat(2,i)=ZERO
      end do
c
      call fft( tmpdat,nft,+1 )
c
c     --------------------------------------------------
c     Convolute with Gaussian function by multiplying in
c     with a Gaussian in Fourier space.
c     --------------------------------------------------
      no2=nft/2
      df=TWO*PI/(dfloat(nft)*dfld)
      gnorm=ONE/dfloat(nft)
      tmpdat(1,1)=tmpdat(1,1)*gnorm
      tmpdat(2,1)=tmpdat(2,1)*gnorm
      f=df
      do i=2,no2
         k=nft-i+2
         g=gnorm*dexp( -(f*wline)**2/EIGHT )
         tmpdat(1,i)=tmpdat(1,i)*g
         tmpdat(2,i)=tmpdat(2,i)*g
         tmpdat(1,k)=tmpdat(1,k)*g
         tmpdat(2,k)=tmpdat(2,k)*g
         f=f+df
      end do
      k=no2+1
      g=gnorm*dexp( -(f*wline)**2/EIGHT )
      tmpdat(1,k)=tmpdat(1,k)*g
      tmpdat(2,k)=tmpdat(2,k)*g
c
      call fft(tmpdat,nft,-1)
c
      do i=1,nfld
         spectr(1,i)=tmpdat(1,i)
         spectr(2,i)=tmpdat(2,i)
      end do
c
      return
      end
