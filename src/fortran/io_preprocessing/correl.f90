! NLS Version 1.3  7/13/93
!----------------------------------------------------------------------
!                      =====================
!                        subroutine CORREL 
!                      =====================
!
!   Taken from subroutine CORREL in Numerical Recipes, Press et al.
!   Calculates the correlation function of two input arrays.
!
!
!   Inputs:
!     data     Data array
!     spctr    Array that is to be shifted to match data array
!     npt      Number of points in data and spctr 
!              NOTE: For the FFT, this must be a power of two.
!              This routine does not check for this condition.
!     ans      Contains correlation function (should be dimensioned
!              as complex*16 or twice the size of data array)
!     fft      Work array dimensioned as complex*16(npt)
!
!     Uses:
!       twofft   (from Numerical Recipes)
!       four1    (from Numerical Recipes)
!       realft   (from Numerical Recipes)
!   
!----------------------------------------------------------------------
      subroutine correl(data,spctr,npt,ans,fft)
      implicit none
      integer npt
      double precision data(npt),spctr(npt),ans(2*npt),fft(2*npt)
!
      integer i,ii,ir,no2
      double precision ansr,halfn,tmp
!
!  ----
        call twofft(data,spctr,fft,ans,npt)
        no2=npt/2
        halfn=float(no2)
        do i=1,no2+1
           ii=i+i       
           ir=ii-1
           ansr=ans(ir)
           ans(ir)=(fft(ir)*ansr+fft(ii)*ans(ii))/halfn
           ans(ii)=(fft(ii)*ansr-fft(ir)*ans(ii))/halfn
       end do
       ans(2)=ans(npt+1)
       call realft(ans,no2,-1)
! ---
!
!
!  Note that correlation function is stored from ans(1) to ans(no2)
!  for zero and positive correlations, and from ans(npt) down to 
!  ans(no2+1) for negative correlations. Rearrange the ans array 
!  according to increasing correlations.
!
        do i=1,no2
           tmp=ans(i)
           ans(i)=ans(no2+i)
           ans(no2+i)=tmp
        end do
        return
!
        end

