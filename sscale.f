c  NLSL Version 1.3 7/17/93
c----------------------------------------------------------------------
c                    =========================
c                      subroutine SSCALE
c                    =========================
c
c   Given the array <data> with <npt> data points and a matrix of <nsite>
c   calculated spectra, each with <npt> points, calculate a least-squares
c   set of <nsite> scaling factors. This routine uses the QR factorization
c   routine that belongs to MINPACK together with the associated routines
c   for solving the least-squares problem using QR factoriaation.
c
c   Inputs:
c       data(npt)          Data values to be fit by scaled spectra
c       spectr(ldspct,nsi) Spectrum/spectra to be scaled
c       ldspct             Leading dimension of spectr array
c       nsite              # spectra contained in spectr array
c       npt                # points in data/spectra arrays
c       ctol               Tolerance for linearly independent components
c       iscal              Flags =1 for automatic scaling of each site
c                          (otherwise used fixed value passed in sfac)
c       noneg              Flag =1: zero any negative coeffients
c       work(npt)          Temporary work array used in QR solution
c
c    Output:
c       sfac               Array of scaling factors
c                          Scale factors corresponding to spectra that
c                          were eliminated on the basis of linear dependence
c                          or were negative (if noneg=1)
c
c       resid              Residuals (data minus scaled calculated spectra)
c             
c   Includes:
c     nlsdim.inc
c 
c   Uses:
c     qrfac    Calculate QR factorization of design matrix (MINPACK)
c     qtbvec   Calculate Q(transpose)*b from Householder transformations
c     rsolve   Solve R*x = Q(transpose)*b
c     dchex    Update QR factors for column permutation of design matrix
c  
c---------------------------------------------------------------------- 
      subroutine sscale( data,spct,wspct,ldspct,work,nsite,npt,ctol,
     #                   noneg,iscal,sfac,resid )
      implicit none
c
      integer ldspct,noneg,npt,nsite
      integer iscal(nsite)
      double precision data(npt),work(npt),spct(ldspct,nsite),
     #                 wspct(ldspct,nsite),resid(npt),sfac(nsite),ctol
c
      include 'nlsdim.inc'
c
      integer ixspc(MXSITE),jpvt(MXSITE)
      double precision qtd(MXSITE),qraux(MXSITE),rdiag(MXSITE),
     #                 tmpfac(MXSITE),wa1(MXSITE),wa2(MXSITE)
c
      integer i,info,j,jtmp,k,m,mneg,nscl
      double precision smin,sumc2,sumdc,tmp
c
      double precision ZERO
      parameter (ZERO=0.0d0)
c
c
c######################################################################
c
      do j=1,npt
         resid(j)=data(j)
      end do
c
c     --------------------------------------------------------------
c     Copy any spectra with the autoscale flag set to the work array
c     Subtract spectra with fixed scale factors from the residuals
c     --------------------------------------------------------------
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
c
      if (nscl.eq.0) return
c
c----------------------------------------------------------------------
c     If there is only one site, calculate the least-squares
c     scaling factor directly
c----------------------------------------------------------------------
      if (nscl.eq.1) then
c
         sumdc=ZERO
         sumc2=ZERO
         do i=1,npt
            sumdc=sumdc+resid(i)*wspct(i,1)
            sumc2=sumc2+wspct(i,1)*wspct(i,1)
         end do
c
         if (sumc2.ne.ZERO) then
            tmpfac(1)=sumdc/sumc2
         else
            tmpfac(1)=ZERO
         end if
c
c---------------------------------------------------------------
c     For multiple sites, use QR decomposition to find linear 
c     least-squares scaling coefficients and rms deviation 
c---------------------------------------------------------------
      else
c
c        --------------------------------------------
c        Compute the QR factorization of the spectra
c        --------------------------------------------
         call qrfac(npt,nscl,wspct,ldspct,.true.,jpvt,nscl,rdiag,
     #        wa1,wa2)
c
         do i=1,nscl
            qraux(i)=wspct(i,i)
            wspct(i,i)=rdiag(i)
         end do
c
c        --------------------------------------------------------------
c         Determine which spectra are linearly dependent and remove them
c         from the problem
c        --------------------------------------------------------------
         k=0
         do i=1,nscl
            if (dabs(rdiag(i)) .le. ctol*dabs(rdiag(1)) ) goto 14
            k=i
         end do
c
 14      call qtbvec(npt,k,wspct,ldspct,qraux,data,work)
c
c       ---------------------------------------------------------
c        Loop to here if QR factorization has been updated 
c        to eliminatespectra entering with negative coefficients
c       ---------------------------------------------------------
c
 15      call rsolve( npt,k,wspct,ldspct,qraux,work,tmpfac,
     #               .false.,tmp ) 
c
c        ------------------------------------------------------------    
c         Optional check for negative coefficients: permute spectra
c         with negative coefficients into the truncated part of the QR 
c         factorization and re-solve the problem
c         (But don't truncate the last spectrum!)
c        ------------------------------------------------------------    
         if (noneg.ne.0.and.k.gt.1) then
            smin=ZERO
            mneg=0
            do i=1,k
               if (tmpfac(i).lt.smin) then
                  mneg=i
                  smin=tmpfac(i)
               end if
            end do
c
            if (mneg.ne.0) then
               if (mneg.lt.k) then
                  call dchex( wspct,ldspct,nscl,mneg,k,
     #                        work,ldspct,1,wa1,wa2,2 )
c     
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
c
c      ---------------------------------------------------------------
c      Set the unused components of tmpfac to zero and initialize jpvt
c      for unscrambling scaling coefficients from pivoted order
c      ---------------------------------------------------------------
         do i=1,nscl
            jpvt(i)=-jpvt(i)
            if (i.gt.k) tmpfac(i)=ZERO
         end do
c
c      -----------------------------------------------------
c      Unscramble the solution from pivoted order using jpvt
c      -----------------------------------------------------
         do 70 i=1,nscl
            if (jpvt(i).le.0) then
               k=-jpvt(i)
               jpvt(i)=k
c     
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
c     
c             *** End if (nscl.eq.1) ... else ...
      end if
c
c     ------------------------------------------------------
c     Calculate residuals using least-squares scale factors
c     ------------------------------------------------------
      do i=1,nscl
         k=ixspc(i)
         sfac(k)=tmpfac(i)
         do j=1,npt
            resid(j)=resid(j)-sfac(k)*spct(j,k)
         end do
      end do
c
      return
      end
