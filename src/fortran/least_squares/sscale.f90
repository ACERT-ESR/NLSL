!  NLSL Version 1.3 7/17/93
!----------------------------------------------------------------------
!                    =========================
!                      subroutine SSCALE
!                    =========================
!
!   Given the array <data> with <npt> data points and a matrix of <nsite>
!   calculated spectra, each with <npt> points, calculate a least-squares
!   set of <nsite> scaling factors. This routine uses the QR factorization
!   routine that belongs to MINPACK together with the associated routines
!   for solving the least-squares problem using QR factoriaation.
!
!   Inputs:
!       data(npt)          Data values to be fit by scaled spectra
!       spectr(ldspct,nsi) Spectrum/spectra to be scaled
!       ldspct             Leading dimension of spectr array
!       nsite              # spectra contained in spectr array
!       npt                # points in data/spectra arrays
!       ctol               Tolerance for linearly independent components
!       iscal              Flags =1 for automatic scaling of each site
!                          (otherwise used fixed value passed in sfac)
!       noneg              Flag =1: zero any negative coeffients
!       work(npt)          Temporary work array used in QR solution
!
!    Output:
!       sfac               Array of scaling factors
!                          Scale factors corresponding to spectra that
!                          were eliminated on the basis of linear dependence
!                          or were negative (if noneg=1)
!
!       resid              Residuals (data minus scaled calculated spectra)
!             
!   Includes:
!     nlsdim.inc
! 
!   Uses:
!     qrfac    Calculate QR factorization of design matrix (MINPACK)
!     qtbvec   Calculate Q(transpose)*b from Householder transformations
!     rsolve   Solve R*x = Q(transpose)*b
!     dchex    Update QR factors for column permutation of design matrix
!  
!---------------------------------------------------------------------- 
      subroutine sscale( data,spct,wspct,ldspct,work,nsite,npt,ctol,
     #                   noneg,iscal,sfac,resid )
!
      use nlsdim
!
      implicit none
!
      integer ldspct,noneg,npt,nsite
      integer iscal(nsite)
      double precision data(npt),work(npt),spct(ldspct,nsite),
     #                 wspct(ldspct,nsite),resid(npt),sfac(nsite),ctol
!
      integer ixspc(MXSITE),jpvt(MXSITE)
      double precision qraux(MXSITE),rdiag(MXSITE),
     #                 tmpfac(MXSITE),wa1(MXSITE),wa2(MXSITE)
!
      integer i,j,jtmp,k,mneg,nscl
      double precision smin,sumc2,sumdc,tmp
!
      double precision ZERO
      parameter (ZERO=0.0d0)
!
!
!######################################################################
!
      do j=1,npt
         resid(j)=data(j)
      end do
!
!     --------------------------------------------------------------
!     Copy any spectra with the autoscale flag set to the work array
!     Subtract spectra with fixed scale factors from the residuals
!     --------------------------------------------------------------
      nscl=0
      do i=1,nsite
         if (iscal(i).ne.0) then
            nscl=nscl+1
            ixspc(nscl)=i
            do j=1,npt
               wspct(j,nscl)=spct(j,i)
            end do
         else
            do j=1,npt
               resid(j)=resid(j)-sfac(i)*spct(j,i)
            end do
         end if
      end do
!
      if (nscl.eq.0) return
!
!----------------------------------------------------------------------
!     If there is only one site, calculate the least-squares
!     scaling factor directly
!----------------------------------------------------------------------
      if (nscl.eq.1) then
!
         sumdc=ZERO
         sumc2=ZERO
         do i=1,npt
            sumdc=sumdc+resid(i)*wspct(i,1)
            sumc2=sumc2+wspct(i,1)*wspct(i,1)
         end do
!
         if (sumc2.ne.ZERO) then
            tmpfac(1)=sumdc/sumc2
         else
            tmpfac(1)=ZERO
         end if
!
!---------------------------------------------------------------
!     For multiple sites, use QR decomposition to find linear 
!     least-squares scaling coefficients and rms deviation 
!---------------------------------------------------------------
      else
!
!        --------------------------------------------
!        Compute the QR factorization of the spectra
!        --------------------------------------------
         call qrfac(npt,nscl,wspct,ldspct,.true.,jpvt,nscl,rdiag,
     #        wa1,wa2)
!
         do i=1,nscl
            qraux(i)=wspct(i,i)
            wspct(i,i)=rdiag(i)
         end do
!
!        --------------------------------------------------------------
!         Determine which spectra are linearly dependent and remove them
!         from the problem
!        --------------------------------------------------------------
         k=0
         do i=1,nscl
            if (dabs(rdiag(i)) .le. ctol*dabs(rdiag(1)) ) goto 14
            k=i
         end do
!
 14      call qtbvec(npt,k,wspct,ldspct,qraux,data,work)
!
!       ---------------------------------------------------------
!        Loop to here if QR factorization has been updated 
!        to eliminatespectra entering with negative coefficients
!       ---------------------------------------------------------
!
 15      call rsolve( npt,k,wspct,ldspct,qraux,work,tmpfac,
     #               .false.,tmp ) 
!
!        ------------------------------------------------------------    
!         Optional check for negative coefficients: permute spectra
!         with negative coefficients into the truncated part of the QR 
!         factorization and re-solve the problem
!         (But don't truncate the last spectrum!)
!        ------------------------------------------------------------    
         if (noneg.ne.0.and.k.gt.1) then
            smin=ZERO
            mneg=0
            do i=1,k
               if (tmpfac(i).lt.smin) then
                  mneg=i
                  smin=tmpfac(i)
               end if
            end do
!
            if (mneg.ne.0) then
               if (mneg.lt.k) then
                  call dchex( wspct,ldspct,nscl,mneg,k,
     #                        work,ldspct,1,wa1,wa2,2 )
!     
                  jtmp=jpvt(mneg)
                  do j=mneg,k-1
                     jpvt(j)=jpvt(j+1)
                  end do
                  jpvt(k)=jtmp
               end if
               k=k-1
               go to 15
            end if
         end if
!
!      ---------------------------------------------------------------
!      Set the unused components of tmpfac to zero and initialize jpvt
!      for unscrambling scaling coefficients from pivoted order
!      ---------------------------------------------------------------
         do i=1,nscl
            jpvt(i)=-jpvt(i)
            if (i.gt.k) tmpfac(i)=ZERO
         end do
!
!      -----------------------------------------------------
!      Unscramble the solution from pivoted order using jpvt
!      -----------------------------------------------------
         do 70 i=1,nscl
            if (jpvt(i).le.0) then
               k=-jpvt(i)
               jpvt(i)=k
!     
 50            if (k.eq.i) goto 70
               tmp=tmpfac(i)
               tmpfac(i)=tmpfac(k)
               tmpfac(k)=tmp
               tmp=rdiag(i)
               rdiag(i)=rdiag(k)
               rdiag(k)=tmp
               jpvt(k)=-jpvt(k)
               k=jpvt(k)
               go to 50
            end if
 70      continue
!     
!             *** End if (nscl.eq.1) ... else ...
      end if
!
!     ------------------------------------------------------
!     Calculate residuals using least-squares scale factors
!     ------------------------------------------------------
      do i=1,nscl
         k=ixspc(i)
         sfac(k)=tmpfac(i)
         do j=1,npt
            resid(j)=resid(j)-sfac(k)*spct(j,k)
         end do
      end do
!
      return
      end
