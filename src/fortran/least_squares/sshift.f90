! NLSL Version 1.5 7/1/95
!----------------------------------------------------------------------
!                     =========================
!                        subroutine SSHIFT
!                     =========================
!
!     Routine to calculate the shift needed to obtain optimal overlap
!     between a spectrum (or set of spectra) and experimental data
!
!     Inputs:
!       data(npt)          Data values to be fit by shifting/scaling
!       spectr(ldspct,nsi) Spectrum/spectra to be shifted/scaled
!       ldspct             Leading dimension of spectr array
!       nsi                # spectra contained in spectr array
!       npt                # points in data/spectra arrays
!       nft                Array length for FFT (next power of 2 > npt)
!       ideriv             Flag =1 for 1st derivative data
!       srange             Allowed range for shifting (fraction of npt)
!       ctol               Tolerance for linearly independent components
!       noneg              Flag =1: zero any negative coeffients
!       tmpdat(2*nft)      Temporary storage for 0-padded data (for FFT) 
!       tmpclc(2*nft)      Temporary storage for 0-padded spectrum (for FFT) 
!       wspec(ldspct,nsi)  Work array used to store correlation functions
!       work(2*nft)        Temporary work array used in FFT and QR solution
!
!     Output:
!       sfac               Scale factors for each site
!             
!     Includes:
!      nlsdim.inc
!
!     Uses:
!       correl Calculate correlation function (from Numerical Recipes)
!       qrfac  QR factorization of normal equations (part of MINPACK)
!       rsolve Solves the equation R*x = Q(transpose)*b
!       qtbvec Calculates the vector Q(transpose)*b from the Householder
!              factors of Q
!       dchex  Updates QR factorization (from LINPACK)   
!       fmomnt Returns first moment of an array (in index units)
!
!----------------------------------------------------------------------
      function sshift( data,spct,ldspct,nsi,npt,nft,ideriv,srange,
     #     ctol,noneg,tmpdat,tmpclc,wspec,work,sfac )
!
      use nlsdim
!
      implicit none
      double precision sshift
!
      integer iarr,ideriv,ldspct,nsi,npt,nft,noneg
      double precision spct(ldspct,nsi),data(npt),tmpdat(2*nft),
     #                 tmpclc(2*nft),wspec(ldspct,nsi),work(2*nft),
     #                 sfac(nsi),srange,ctol
!
      integer i,info,irng,isi,ixmx,ixw1,ixw2,ixw3,j,jtmp,k,mneg,
     #        mnrng,mxrng
      integer jpvt(MXSITE)
      double precision a,approx,asum,b,c,c1m,d1m,dummy,ovmax,
     #                 scl,shift,smax,smin,temp,xfrac
      double precision amat(MXSITE,MXSITE),qraux(MXSITE),
     #                rdiag(MXSITE),tmpscl(MXSITE)
!
      double precision ZERO
      parameter(ZERO=0.0D0)
!
!     double precision fmomnt,lgrint
      double precision fmomnt
      double complex lgrint
      external fmomnt,lgrint
!
!######################################################################
!
!----------------------------------------------------------------------
!     Determine the allowed range for shifting. First, find an 
!     approximate shift using the first moments of the experimental 
!     and calculated spectra
!----------------------------------------------------------------------
      smax = dabs( srange*npt )
      d1m=dabs ( fmomnt( data,npt,ideriv ) )
      c1m=ZERO
      do i=1,nsi
         c1m=c1m+dabs ( fmomnt( spct(1,i),npt,ideriv ) )
      end do
!
      approx = d1m-c1m/dble(nsi)
      if (approx.gt.smax .or. approx.lt.-smax) approx=ZERO
      mnrng = max0( int(approx-smax)+nft/2, 1 )
      mxrng = min0( int(approx+smax)+nft/2, npt )
!
!----------------------------------------------------------------------
!  ***  Make a copy of the data and zero-pad it for FFT
!      (A temporary array is needed to accommodate zero-padding)
!----------------------------------------------------------------------
      do j=1,npt
         tmpdat(j)=data(j)
      end do
      do j=npt+1,nft
         tmpdat(j)=ZERO
      end do
!
!----------------------------------------------------------------------
!  *** Loop over each site for this spectrum
!----------------------------------------------------------------------
!
!  *** Copy each calculated spectrum to temporary array and zero-pad it
!      (A temporary array is needed to accommodate zero-padding)
!
      do isi=1,nsi
         do j=1,npt
            tmpclc(j)=spct(j,isi)
         end do
         do j=npt+1,nft
            tmpclc(j)=ZERO
         end do
!
         call correl( tmpdat,tmpclc,nft,wspec(1,isi),work )
!
      end do
!
!----------------------------------------------------------------------
!     *** Calculate the normal equations matrix for linear least-squares 
!----------------------------------------------------------------------
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
!
!----------------------------------------------------------------------
!     Calculate the QR decomposition of the normal equation matrix
!----------------------------------------------------------------------
      ixw1=1+nsi
      ixw2=ixw1+nsi
      ixw3=ixw2+nsi
      call qrfac(nsi,nsi,amat,MXSITE,.true.,jpvt,nsi,rdiag,
     #     work(ixw1),work(ixw2) )
!
!     ---------------------------------------------------------------
!     Store diagonal of Q (returned in diagonal of A) in work vector
!     Replace diagonal of A with diagonal of R
!     (this is the arrangement expected by qtbvec and rsolve)
!     ---------------------------------------------------------------
      do i=1,nsi
         work(i)=amat(i,i)
         amat(i,i)=rdiag(i)
      end do
!
!     -----------------------------------------------------------------
!     Determine which spectra are linearly dependent and truncate them
!     -----------------------------------------------------------------
      k=0
      do i=1,nsi
         if (dabs(amat(i,i)) .le. ctol*dabs(amat(1,1)) ) go to 21
         k=i
      end do
!
!----------------------------------------------------------------------
!      Calculate the correlation function of the data with the
!      least-squares sum of spectra by solving the (possibly truncated) 
!      linear least-squares problem for each shift value in the 
!      allowed shifting range. These will form the vector part of
!      the normal equations at the given shift value. 
!      Store the correlation functions in the tmpclc array.
!----------------------------------------------------------------------
!
!    ----------------------------------------------------------------
!     Loop through all possible shift values in discrete correlation
!     function(s)
!    ----------------------------------------------------------------
 21   do irng=mnrng,mxrng
!
!        ------------------------------------------------------------
!         Store in tmpdat the values of the correlation functions of 
!         each individual site spectrum with the data for the current
!         shift value 
!        ------------------------------------------------------------
         do j=1,nsi
            tmpdat(j)=wspec(irng,j)
         end do
!
!        -----------------------------------------------------------
!         Calculate Q(trans)*vector (returned in work(1) array) and 
!         solve for scaling coefficients 
!        -----------------------------------------------------------
         call qtbvec( nsi,k,amat,MXSITE,work,tmpdat,work(ixw1) )
         call rsolve( nsi,k,amat,MXSITE,work,work(ixw1),tmpscl,
     #               .false.,work(ixw1) ) 
!
!        ------------------------------------------------------------
!        Calculate correlation function using optimal scale factors 
!        for the current shift value
!        ------------------------------------------------------------
         tmpclc(irng)=ZERO
         do j=1,nsi
            scl=tmpscl( jpvt(j) )
            if (j.gt.k) scl=ZERO
            tmpclc(irng)=tmpclc(irng)+scl*wspec(irng,j)
         end do
!     
      end do
!
!--------------------------------------------------------------------------
!     Search correlation function stored in tmpclc for the maximum value 
!     within the allowed shifting range. 
!--------------------------------------------------------------------------
      ixmx=mnrng
      ovmax=ZERO
      do irng=mnrng,mxrng
         if (tmpclc(irng).gt.ovmax) then
            ixmx=irng
            ovmax=tmpclc(irng)
         end if
      end do
!
!     ------------------------------------------------------------------
!     Find the ordinate for the extremum (in this case a known maximum)
!     ------------------------------------------------------------------
!      if (ixmx.gt.mnrng .and. ixmx.lt.mxrng) then
!         xfrac = -(tmpclc(ixmx+1)-tmpclc(ixmx-1))/
!     #        (tmpclc(ixmx+1)+tmpclc(ixmx-1)-2.0d0*tmpclc(ixmx))
!      else
!         xfrac=ZERO
!      end if
!
      iarr=max(ixmx-2,mnrng)
      iarr=min(iarr,mxrng-4)
      xfrac=dble(lgrint(tmpclc(iarr)))
!
      shift=dble(ixmx-(nft/2)-1)+xfrac
!
!    ----------------------------------------------------------------------
!      Find the values of the scaling coefficients for the calculated
!      fractional shift: place interpolated values of correlation functions
!      in tmpdat
!    ----------------------------------------------------------------------
      if (ixmx.gt.mnrng.and. ixmx.lt.mxrng) then
         do j=1,nsi
            a = wspec(ixmx,j)
            b = 0.5d0*(wspec(ixmx+1,j)-wspec(ixmx-1,j))
            c = 0.5d0*(wspec(ixmx+1,j)+wspec(ixmx-1,j))-a
            tmpdat(j)= a + b*xfrac + c*xfrac*xfrac
         end do
!
      else
         do j=1,nsi
            tmpdat(j)=wspec(ixmx,j)
         end do
      end if
!
!     --------------------------------------------------
!     Now solve the linear least-squares one more time
!     to find the spectral scaling coefficients at the 
!     correlation function maximum
!     --------------------------------------------------
      call qtbvec( nsi,k,amat,MXSITE,work,tmpdat,work(ixw1) )
      call rsolve( nsi,k,amat,MXSITE,work,work(ixw1),tmpscl,.false.,
     #     dummy ) 
!
!    ------------------------------------------------------------
!    Optional check for negative coefficients: permute spectra
!     with negative coefficients into the truncated part of the QR 
!     factorization and re-solve the problem
!
!     (But don't truncate the last spectrum!)
!    ------------------------------------------------------------
!            
      if (noneg.ne.0.and.k.gt.1) then
         smin=ZERO
         mneg=0
         do i=1,k
            if (tmpscl(i).lt.smin) then
               mneg=i
               smin=tmpscl(i)
            end if
         end do
!
         if (mneg.ne.0) then
            if (mneg.lt.k) then
               call dchex(amat,MXSITE,nsi,mneg,k,work(ixw1),
     #              MXSITE,1,work(ixw2),work(ixw3),2)
!     
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
!
!     ------------------------------------------------------------
!      Zero the unused components of sfac and initialize the 
!      permutation array jpvt for unscrambling from pivoted order
!     ------------------------------------------------------------
      do i=1,nsi
         jpvt(i) = -jpvt(i)
         sfac(i) = tmpscl(i)
         if (i.gt.k) sfac(i)=ZERO
      end do
!
!     ------------------------------------------------------------
!      Unscramble the solution from its pivoted order using jpvt
!     ------------------------------------------------------------
      do 70 i=1,nsi
         if (jpvt(i).le.0) then
            k=-jpvt(i)
            jpvt(i)=k
!     
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
!
      sshift=shift
      return
      end




