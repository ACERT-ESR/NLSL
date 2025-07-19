c Version 1.3.2
c----------------------------------------------------------------------
c                    =========================
c                       subroutine GCONVL
c                    =========================
c
c     Convolute given real-valued spectrum with a Gaussian lineshape 
c     having a specified derivative peak-to-peak linewidth.
c
c     NOTE: this differs from EPRLL version, which transforms a
c     complex-valued spectrum.      
c
c----------------------------------------------------------------------
      subroutine gconvl( spectr,wline,dfld,nfld,nft )
c
      use nlsdim
      use ftwork
      use pidef
c
      implicit none
c
      integer nfld,nft
      double precision spectr(nfld),wline,dfld
c
      integer i,no2
      double precision df,f,g,gnorm
c
      double precision EIGHT,EPS,ONE,TWO,ZERO
      parameter ( EIGHT=8.0D0,ONE=1.0D0,TWO=2.0D0,ZERO=0.0D0,
     #            EPS=1.0D-6 )
c
c######################################################################
c
      if (wline.lt.EPS) return
c
c     -----------------------------------------------
c     Store calculation in tmpdat, zero-pad, and FFT
c     -----------------------------------------------
      do i=1,nfld
         tmpdat(i)=spectr(i)
      end do
c
      do i=nfld+1,nft
         tmpdat(i)=ZERO
      end do
c
      no2=nft/2
      call realft( tmpdat,no2,+1 )
c
c     ------------------------------------------------------------
c     Convolute with Gaussian function by multiplying in
c     with a Gaussian in Fourier space.
c
c     NOTE: REALFT returns only the positive frequencies
c     since FT of a real function is symmetric. Also, the
c     first and last elements of the FT array are real-valued
c     and returned as the real and imaginary parts of tmpdat(1)
c     ------------------------------------------------------------
      df=TWO*PI/(nft*dfld)
      gnorm=ONE/dfloat(no2)
      tmpdat(1)=tmpdat(1)*gnorm
      tmpdat(2)=tmpdat(2)*gnorm*dexp(-(no2*df*wline)**2/EIGHT)
      f=df
      do i=3,nft,2
         g=gnorm*dexp( -(f*wline)**2/EIGHT )
         tmpdat(i)=tmpdat(i)*g
         tmpdat(i+1)=tmpdat(i+1)*g
         f=f+df
      end do
c
c     ------------------------------------------------------------
c     Back-transfor to obtain convoluted spectrum and restore in
c     spectr array
c     ------------------------------------------------------------
      call realft( tmpdat,no2,-1 )
c
      do i=1,nfld
         spectr(i)=tmpdat(i)
      end do
c
      return
      end
