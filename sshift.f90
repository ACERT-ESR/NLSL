c NLSL Version 1.5 7/1/95
c----------------------------------------------------------------------
c                     =========================
c                        subroutine SSHIFT
c                     =========================
c
c     Routine to calculate the shift needed to obtain optimal overlap
c     between a spectrum (or set of spectra) and experimental data
c
c     Inputs:
c       data(npt)          Data values to be fit by shifting/scaling
c       spectr(ldspct,nsi) Spectrum/spectra to be shifted/scaled
c       ldspct             Leading dimension of spectr array
c       nsi                # spectra contained in spectr array
c       npt                # points in data/spectra arrays
c       nft                Array length for FFT (next power of 2 > npt)
c       ideriv             Flag =1 for 1st derivative data
c       srange             Allowed range for shifting (fraction of npt)
c       ctol               Tolerance for linearly independent components
c       noneg              Flag =1: zero any negative coeffients
c       tmpdat(2*nft)      Temporary storage for 0-padded data (for FFT) 
c       tmpclc(2*nft)      Temporary storage for 0-padded spectrum (for FFT) 
c       wspec(ldspct,nsi)  Work array used to store correlation functions
c       work(2*nft)        Temporary work array used in FFT and QR solution
c
c     Output:
c       sfac               Scale factors for each site
c             
c     Includes:
c      nlsdim.inc
c
c     Uses:
c       correl Calculate correlation function (from Numerical Recipes)
c       qrfac  QR factorization of normal equations (part of MINPACK)
c       rsolve Solves the equation R*x = Q(transpose)*b
c       qtbvec Calculates the vector Q(transpose)*b from the Householder
c              factors of Q
c       dchex  Updates QR factorization (from LINPACK)   
c       fmomnt Returns first moment of an array (in index units)
c
c----------------------------------------------------------------------
      function sshift( data,spct,ldspct,nsi,npt,nft,ideriv,srange,
     #     ctol,noneg,tmpdat,tmpclc,wspec,work,sfac )
c
      use nlsdim
c
      implicit none
      double precision sshift
c
      integer iarr,ideriv,ldspct,nsi,npt,nft,noneg
      double precision spct(ldspct,nsi),data(npt),tmpdat(2*nft),
     #                 tmpclc(2*nft),wspec(ldspct,nsi),work(2*nft),
     #                 sfac(nsi),srange,ctol,err
c
      integer i,irng,isi,ixmx,ixw1,ixw2,ixw3,j,jtmp,k,mneg,
     #        mnrng,mxrng
      integer jpvt(MXSITE)
      double precision a,approx,asum,b,c,c1m,d1m,dummy,ovmax,
     #                 scl,shift,smax,smin,temp,xfrac
      double precision amat(MXSITE,MXSITE),qraux(MXSITE),
     #                rdiag(MXSITE),tmpscl(MXSITE)
c
      double precision ZERO
      parameter(ZERO=0.0D0)
c
c     double precision fmomnt,lgrint
      double precision fmomnt
      double complex lgrint
      external fmomnt,lgrint
c
c######################################################################
c
c----------------------------------------------------------------------
c     Determine the allowed range for shifting. First, find an 
c     approximate shift using the first moments of the experimental 
c     and calculated spectra
c----------------------------------------------------------------------
      smax = dabs( srange*npt )
      d1m=dabs ( fmomnt( data,npt,ideriv ) )
      c1m=ZERO
      do i=1,nsi
         c1m=c1m+dabs ( fmomnt( spct(1,i),npt,ideriv ) )
      end do
c
      approx = d1m-c1m/dble(nsi)
      if (approx.gt.smax .or. approx.lt.-smax) approx=ZERO
      mnrng = max0( int(approx-smax)+nft/2, 1 )
      mxrng = min0( int(approx+smax)+nft/2, npt )
c
c----------------------------------------------------------------------
c  ***  Make a copy of the data and zero-pad it for FFT
c      (A temporary array is needed to accommodate zero-padding)
c----------------------------------------------------------------------
      do j=1,npt
         tmpdat(j)=data(j)
      end do
      do j=npt+1,nft
         tmpdat(j)=ZERO
      end do
c
c----------------------------------------------------------------------
c  *** Loop over each site for this spectrum
c----------------------------------------------------------------------
c
c  *** Copy each calculated spectrum to temporary array and zero-pad it
c      (A temporary array is needed to accommodate zero-padding)
c
      do isi=1,nsi
         do j=1,npt
            tmpclc(j)=spct(j,isi)
         end do
         do j=npt+1,nft
            tmpclc(j)=ZERO
         end do
c
         call correl( tmpdat,tmpclc,nft,wspec(1,isi),work )
c
      end do
c
c----------------------------------------------------------------------
c     *** Calculate the normal equations matrix for linear least-squares 
c----------------------------------------------------------------------
      do j=1,nsi
         do k=j,nsi
            asum=ZERO
            do i=1,npt
               asum=asum+spct(i,j)*spct(i,k)
            end do
            amat(j,k)=asum
            if (k.ne.j) amat(k,j)=asum
         end do
      end do
c
c----------------------------------------------------------------------
c     Calculate the QR decomposition of the normal equation matrix
c----------------------------------------------------------------------
      ixw1=1+nsi
      ixw2=ixw1+nsi
      ixw3=ixw2+nsi
      call qrfac(nsi,nsi,amat,MXSITE,.true.,jpvt,nsi,rdiag,
     #     work(ixw1),work(ixw2) )
c
c     ---------------------------------------------------------------
c     Store diagonal of Q (returned in diagonal of A) in work vector
c     Replace diagonal of A with diagonal of R
c     (this is the arrangement expected by qtbvec and rsolve)
c     ---------------------------------------------------------------
      do i=1,nsi
         work(i)=amat(i,i)
         amat(i,i)=rdiag(i)
      end do
c
c     -----------------------------------------------------------------
c     Determine which spectra are linearly dependent and truncate them
c     -----------------------------------------------------------------
      k=0
      do i=1,nsi
         if (dabs(amat(i,i)) .le. ctol*dabs(amat(1,1)) ) go to 21
         k=i
      end do
c
c----------------------------------------------------------------------
c      Calculate the correlation function of the data with the
c      least-squares sum of spectra by solving the (possibly truncated) 
c      linear least-squares problem for each shift value in the 
c      allowed shifting range. These will form the vector part of
c      the normal equations at the given shift value. 
c      Store the correlation functions in the tmpclc array.
c----------------------------------------------------------------------
c
c    ----------------------------------------------------------------
c     Loop through all possible shift values in discrete correlation
c     function(s)
c    ----------------------------------------------------------------
 21   do irng=mnrng,mxrng
c
c        ------------------------------------------------------------
c         Store in tmpdat the values of the correlation functions of 
c         each individual site spectrum with the data for the current
c         shift value 
c        ------------------------------------------------------------
         do j=1,nsi
            tmpdat(j)=wspec(irng,j)
         end do
c
c        -----------------------------------------------------------
c         Calculate Q(trans)*vector (returned in work(1) array) and 
c         solve for scaling coefficients 
c        -----------------------------------------------------------
         call qtbvec( nsi,k,amat,MXSITE,work,tmpdat,work(ixw1) )
         call rsolve( nsi,k,amat,MXSITE,work,work(ixw1),tmpscl,
     #               .false.,work(ixw1) ) 
c
c        ------------------------------------------------------------
c        Calculate correlation function using optimal scale factors 
c        for the current shift value
c        ------------------------------------------------------------
         tmpclc(irng)=ZERO
         do j=1,nsi
            scl=tmpscl( jpvt(j) )
            if (j.gt.k) scl=ZERO
            tmpclc(irng)=tmpclc(irng)+scl*wspec(irng,j)
         end do
c     
      end do
c
c--------------------------------------------------------------------------
c     Search correlation function stored in tmpclc for the maximum value 
c     within the allowed shifting range. 
c--------------------------------------------------------------------------
      ixmx=mnrng
      ovmax=ZERO
      do irng=mnrng,mxrng
         if (tmpclc(irng).gt.ovmax) then
            ixmx=irng
            ovmax=tmpclc(irng)
         end if
      end do
c
c     ------------------------------------------------------------------
c     Find the ordinate for the extremum (in this case a known maximum)
c     ------------------------------------------------------------------
c      if (ixmx.gt.mnrng .and. ixmx.lt.mxrng) then
c         xfrac = -(tmpclc(ixmx+1)-tmpclc(ixmx-1))/
c     #        (tmpclc(ixmx+1)+tmpclc(ixmx-1)-2.0d0*tmpclc(ixmx))
c      else
c         xfrac=ZERO
c      end if
c
      iarr=max(ixmx-2,mnrng)
      iarr=min(iarr,mxrng-4)
      xfrac=lgrint(tmpclc(iarr),err)
c
      shift=dble(ixmx-(nft/2)-1)+xfrac
c
c    ----------------------------------------------------------------------
c      Find the values of the scaling coefficients for the calculated
c      fractional shift: place interpolated values of correlation functions
c      in tmpdat
c    ----------------------------------------------------------------------
      if (ixmx.gt.mnrng.and. ixmx.lt.mxrng) then
         do j=1,nsi
            a = wspec(ixmx,j)
            b = 0.5d0*(wspec(ixmx+1,j)-wspec(ixmx-1,j))
            c = 0.5d0*(wspec(ixmx+1,j)+wspec(ixmx-1,j))-a
            tmpdat(j)= a + b*xfrac + c*xfrac*xfrac
         end do
c
      else
         do j=1,nsi
            tmpdat(j)=wspec(ixmx,j)
         end do
      end if
c
c     --------------------------------------------------
c     Now solve the linear least-squares one more time
c     to find the spectral scaling coefficients at the 
c     correlation function maximum
c     --------------------------------------------------
      call qtbvec( nsi,k,amat,MXSITE,work,tmpdat,work(ixw1) )
      call rsolve( nsi,k,amat,MXSITE,work,work(ixw1),tmpscl,.false.,
     #     dummy ) 
c
c    ------------------------------------------------------------
c    Optional check for negative coefficients: permute spectra
c     with negative coefficients into the truncated part of the QR 
c     factorization and re-solve the problem
c
c     (But don't truncate the last spectrum!)
c    ------------------------------------------------------------
c            
      if (noneg.ne.0.and.k.gt.1) then
         smin=ZERO
         mneg=0
         do i=1,k
            if (tmpscl(i).lt.smin) then
               mneg=i
               smin=tmpscl(i)
            end if
         end do
c
         if (mneg.ne.0) then
            if (mneg.lt.k) then
               call dchex(amat,MXSITE,nsi,mneg,k,work(ixw1),
     #              MXSITE,1,work(ixw2),work(ixw3),2)
c     
               jtmp=jpvt(mneg)
               do j=mneg,k-1
                  jpvt(j)=jpvt(j+1)
               end do
               jpvt(k)=jtmp
            end if
            k=k-1
            go to 21
         end if
      end if
c
c     ------------------------------------------------------------
c      Zero the unused components of sfac and initialize the 
c      permutation array jpvt for unscrambling from pivoted order
c     ------------------------------------------------------------
      do i=1,nsi
         jpvt(i) = -jpvt(i)
         sfac(i) = tmpscl(i)
         if (i.gt.k) sfac(i)=ZERO
      end do
c
c     ------------------------------------------------------------
c      Unscramble the solution from its pivoted order using jpvt
c     ------------------------------------------------------------
      do 70 i=1,nsi
         if (jpvt(i).le.0) then
            k=-jpvt(i)
            jpvt(i)=k
c     
 50         if (k.eq.i) goto 70
            temp=sfac(i)
            sfac(i)=sfac(k)
            sfac(k)=temp
            temp=rdiag(i)
            rdiag(i)=rdiag(k)
            rdiag(k)=temp
            jpvt(k)=-jpvt(k)
            k=jpvt(k)
            go to 50
         end if
 70   continue
c
      sshift=shift
      return
      end




