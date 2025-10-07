! Version 1.3.2
!----------------------------------------------------------------------
!                    =========================
!                       subroutine GCONVL
!                    =========================
!
!     Convolute given real-valued spectrum with a Gaussian lineshape 
!     having a specified derivative peak-to-peak linewidth.
!
!     NOTE: this differs from EPRLL version, which transforms a
!     complex-valued spectrum.      
!
!----------------------------------------------------------------------
      subroutine gconvl( spectr,wline,dfld,nfld,nft )
!
      use nlsdim
      use ftwork
      use pidef
!
      implicit none
!
      integer nfld,nft
      double precision spectr(nfld),wline,dfld
!
      integer i,no2
      double precision df,f,g,gnorm
!
      double precision EIGHT,EPS,ONE,TWO,THIRD,ZERO
      parameter ( EIGHT=8.0D0,ONE=1.0D0,TWO=2.0D0,ZERO=0.0D0,
     #            EPS=1.0D-6,THIRD=0.33333333333333D0 )
!
!######################################################################
!
      if (wline.lt.EPS) return
!
!     -----------------------------------------------
!     Store calculation in tmpdat, zero-pad, and FFT
!     -----------------------------------------------
      do i=1,nfld
         tmpdat(i)=spectr(i)
      end do
!
      do i=nfld+1,nft
         tmpdat(i)=ZERO
      end do
!
      no2=nft/2
      call realft( tmpdat,no2,+1 )
!
!     ------------------------------------------------------------
!     Convolute with Gaussian function by multiplying in
!     with a Gaussian in Fourier space.
!
!     NOTE: REALFT returns only the positive frequencies
!     since FT of a real function is symmetric. Also, the
!     first and last elements of the FT array are real-valued
!     and returned as the real and imaginary parts of tmpdat(1)
!     ------------------------------------------------------------
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
!
!     ------------------------------------------------------------
!     Back-transfor to obtain convoluted spectrum and restore in
!     spectr array
!     ------------------------------------------------------------
      call realft( tmpdat,no2,-1 )
!
      do i=1,nfld
         spectr(i)=tmpdat(i)
      end do
!
      return
      end
