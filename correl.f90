c NLS Version 1.3  7/13/93
c----------------------------------------------------------------------
c                      =====================
c                        subroutine CORREL 
c                      =====================
c
c   Taken from subroutine CORREL in Numerical Recipes, Press et al.
c   Calculates the correlation function of two input arrays.
c
c
c   Inputs:
c     data     Data array
c     spctr    Array that is to be shifted to match data array
c     npt      Number of points in data and spctr 
c              NOTE: For the FFT, this must be a power of two.
c              This routine does not check for this condition.
c     ans      Contains correlation function (should be dimensioned
c              as complex*16 or twice the size of data array)
c     fft      Work array dimensioned as complex*16(npt)
c
c     Uses:
c       twofft   (from Numerical Recipes)
c       four1    (from Numerical Recipes)
c       realft   (from Numerical Recipes)
c   
c----------------------------------------------------------------------
      subroutine correl(data,spctr,npt,ans,fft)
      implicit none
      integer npt
      double precision data(npt),spctr(npt),ans(2*npt),fft(2*npt)
c
      integer i,ii,ir,no2
      double precision ansr,halfn,tmp
c
c  ----
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
c ---
c
c
c  Note that correlation function is stored from ans(1) to ans(no2)
c  for zero and positive correlations, and from ans(npt) down to 
c  ans(no2+1) for negative correlations. Rearrange the ans array 
c  according to increasing correlations.
c
        do i=1,no2
           tmp=ans(i)
           ans(i)=ans(no2+i)
           ans(no2+i)=tmp
        end do
        return
c
        end

