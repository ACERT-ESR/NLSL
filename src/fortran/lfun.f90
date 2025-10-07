! NLSL Version 1.5 beta 11/23/95
!----------------------------------------------------------------------
!> @file 
!> This seems to be the central hub, which is called by the L-M
!> optimization, and which calls the functions that generate the
!> simulated spectrum
!                         =======================
!                            subroutine LFUN
!                         =======================
!
!> @brief Subroutine for interfacing EPRLL spectral calculations with the MINPACK 
!>  version of the Levenberg-Marquardt nonlinear least squares algorithm 
!>  for fitting experimental spectra. The function performed
!>  function determined by the IFLAG argument. 
!>  MOMDLS subroutine to calculate a spectrum or series of spectra, each
!>  of which consists of one or more spectral components, or "sites". 
!> 
!> @param iflag
!> @parblock
!> > IFLAG=0 (Printing of iterates)
!> >  
!> > IFLAG=1 (Evaluation of residuals function)
!> > 
!> >   For a function evaluation, this routine takes the input vector X
!> >   of search parameters, loads them into the fparm parameter arrays in 
!> >   /parcom/ (using additional indices available in /parcom/) and then 
!> >   performs spectral calculations by calling the MOMDLS subroutine.
!> > 
!> >   The experimental data to which the calculations are to be compared
!> >   is contained in the common block /expdat/. LFUN optionally shifts
!> >   the calculated spectra to minimize the squared differences between
!> >   calculation and data, and then automatically determines a least-
!> >   squares scale factor (or set of scale factors for multiple sites/spectra)
!> >   which are returned in the vector sfac in common /mspctr/. The 
!> >   The unscaled spectra are returned in the matrix spectr in common
!> >   /mspctr/. Each column of spectr corresponds to a single site;
!> >   multiple spectra are stored sequentially in the columns of spectr.
!> > 
!> >   The residuals (differences between calculation and data after shifting
!> >   and scaling) are scaled by the standard deviation of the estimated 
!> >   noise in each experimental spectrum, and returned to MINPACK through 
!> >   the fvec argument.
!> > 
!> > IFLAG=2 (Jacobian evaluation)
!> > 
!> >   For evaluation of the Jacobian, the routine uses the forward-differences
!> >   approximation to calculate the partial derivative of the spectrum with
!> >   respect to each of the parameters being varied. The routine assumes
!> >   that the individual spectra calculated at the point in parameter space
!> >   where the Jacobian is being evaluated are contained in in the spectr
!> >   matrix, and the appropriate scaling factors in the sfac array (both
!> >   in common /mspctr/ ). The Jacobian matrix is returned in the fjac
!> >   argument.
!> > 
!> > IFLAG=3
!> >     Fit is terminating. Output most recent set of parameter values and
!> >     exit.
!> @endparblock
! Uses:
!    setspc
!    sshift
!    sscale
!    momdls
!    fstplt   (X-Windows interface)
!
!---------------------------------------------------------------------- 
      subroutine single_point(iflag)

      use nlsdim
      use eprprm
      use expdat
      use tridag
      use mspctr
      use ftwork
      use parcom
      use basis
      use errmsg
      use pidef
      use stdio
      use iterat
      use rnddbl

      implicit none
      integer iflag,ise,isi,ixs,icalc,ixb,ixbp,ixt,ierr

      integer FULL,CFONLY
      double precision ZERO
      parameter (FULL=111,CFONLY=0,ZERO=0.0D0)

      logical hltchk
      external hltchk

      ierr=0

      call tdsqz()

!     -----------------------------------
!      Loop over all spectra in a series
!     -----------------------------------
      ixs=1
      do ise=1,nser
         ixsp(ise)=ixs
         tmpshft(ise)=ZERO

         do isi=1,nsite
            call tdchek(isi,ise,ierr)
            call setspc(isi,ise)
            if (hltchk(ierr,isi,ise,iflag)) return

            ixt=ixtd(isi,ise)
            if (basno(isi,ise).gt.0) then
               ixb=ixbas( basno(isi,ise) )
               ixbp=ixb
            else
               ixb=0
               ixbp=1
            end if
            icalc=CFONLY
            if (modtd(isi,ise).ne.0) icalc=FULL

            call momdls( fparm(1,isi),iparm(1,isi),icalc,
     #           alpha(ixt),beta(ixt),ibasis(1,ixbp),ixb,
     #           spectr(ixs,isi),wspec,nft(ise),ltd(isi,ise),
     #           ierr )

            modtd(isi,ise)=0
            if (hltchk(ierr,isi,ise,iflag)) return
         end do

         ixs=ixs+npts(ise)
      end do

      end subroutine single_point

      subroutine lfun( m,n,x,fvec,fjac,ldfjac,iflag )
!
      use nlsdim
      use eprprm
      use expdat
      use tridag      
      use mspctr
      use ftwork
      use parcom
      use basis
      use errmsg
      use pidef
      use stdio
!      use prmeqv
      use iterat
      use rnddbl
!
      implicit none
      integer m,n,iflag,ldfjac,ixbp
!
      double precision x(n),fvec(m),fjac(ldfjac,mxjcol)
      integer i,icalc,ierr,ise,isi,isp,ix,ixb,ixs,ixt,j,k,ld,lastsp,nj
      double precision shift,xtemp
      character dashes*132,shftnd*1
      logical shiftOK,glbscal,sppar
!
      integer FULL,CFONLY
      double precision ZERO
      parameter (FULL=111,CFONLY=0,ZERO=0.0D0)
!
      integer itrim
      double precision sshift,wtdres
      logical spcpar,hltchk
      external sshift,itrim,spcpar,hltchk,wtdres
!
!######################################################################
!
      ierr=0
      ld=0
      glbscal=nsite.gt.1 .and. nspc.gt.1
      if (iter.gt.1) warn=.false.
!
!**********************************************************************
!**********************************************************************
!     IFLAG=0: Output iteration information
!**********************************************************************
!**********************************************************************
!
!
!     ------------------------------------------
!     Draw a line of dashes for first iteration
!     ------------------------------------------
      if (iflag.eq.0) then
         if (iter.eq.1) then
            ld=12*n+16
            if (ld.gt.132) ld=132
            do i=1,ld
               dashes(i:i)='-'
            end do
!
!           -----------------------------------
!           Output names of current parameters 
!           -----------------------------------
!
            if (luout.ne.luttyo) then
               if (ld.gt.0)write (luttyo,1006) dashes(:ld)
               if (iwflag.ne.0) then
                  write (luttyo,1007) (tag(i)(:itrim(tag(i))),i=1,n)
               else
                  write (luttyo,1010) (tag(i)(:itrim(tag(i))),i=1,n)
               end if
               if (ld.gt.0)write (luttyo,1008) dashes(:ld)
            end if
!
            if (ld.gt.0)write (luout,1006)  dashes(:ld)
            if (iwflag.ne.0) then
               write (luout,1007) (tag(i)(:itrim(tag(i))),i=1,n)
            else
               write (luout,1010) (tag(i)(:itrim(tag(i))),i=1,n)
            end if
            if (ld.gt.0)write (luout,1008)  dashes(:ld)
         end if
!
         if (ishglb.ne.0 .and. (nshift.eq.0 .or.iter.le.nshift) ) then
            shftnd='S'
         else
            shftnd=' '
         end if
!
         if (iwflag.ne.0) then
            rdchsq=fnorm*fnorm/float(m-n)
            write (luout,1009) shftnd,iter,rdchsq,(x(i),i=1,n)
            if (luout.ne.luttyo) 
     #        write (luttyo,1009) shftnd,iter,rdchsq,(x(i),i=1,n)
         else
            write (luout,1009) shftnd,iter,fnorm,(x(i),i=1,n)
            if (luout.ne.luttyo) 
     #        write (luttyo,1009) shftnd,iter,fnorm,(x(i),i=1,n)
         end if
!
!        ------------------------------------------------------------
!        Output any other requested information regarding iterates
!        ------------------------------------------------------------
         if (iitrfl.ne.0) call wrfun(iter)
         if (itridg.ne.0) call wrtrdg(iter)
         if (jacobi.ne.0) call wrjac(iter)
!
!**********************************************************************
!**********************************************************************
!     IFLAG=1: Function evaluation
!**********************************************************************
!**********************************************************************
!
      else if (iflag.eq.1) then
!
         xreset=.false.
         written=0
!
         do i=1,n
            call setprm(ixpr(i),ixst(i),x(i))
         end do
!
         shiftOK = (nshift.eq.0 .or. iter.le.nshift)
         call single_point(iflag)
!
!       -----------------------------------
!        Accumulate residuals for each spectrum
!       -----------------------------------
!
         ixs=1
         do ise=1,nser
!
!           --------------------------------------------------------
!           If data are available for this calculation, determine
!           optimal shift and scaling and calculate fvec accordingly
!           (This feature permits calculations w/out data)
!           --------------------------------------------------------
            if (ise.le.nspc) then
!
!             ----------------------------------------------------------
!              If shifting is enabled, find the shift that will give 
!              optimal overlap between a calculated spectrum (or set of 
!              spectra) and the data for each element of a series
!             -----------------------------------------------------------
               if (shiftOK .and. ishft(ise).ne.0) then
!
                  shift=sshift( data(ixs),spectr(ixs,1),MXPT,nsite,
     #                  npts(ise),nft(ise),idrv(ise),srange,ctol,
     #                  noneg,tmpdat,tmpclc,wspec,work,sfac(1,ise) )
                  tmpshft(ise)=shift*sdb(ise)
!
!                 --------------------------------------------------------
!                  Re-calculate shifted spectrum with continued-fractions
!                  (this will keep the frwrd.-diff. derivative reasonable)
!                 --------------------------------------------------------
                  if (dabs(shift).gt.1.0d-3*sdb(ise)) then
                     do isi=1,nsite
                        icalc=CFONLY
                        call setspc(isi,ise)
                        if (hltchk(ierr,isi,ise,iflag)) return
!
                        ixt=ixtd(isi,ise)
                        if (basno(isi,ise).gt.0) then
                          ixb=ixbas( basno(isi,ise) )
                          ixbp=ixb
                        else
                          ixb=0
                          ixbp=1
                        endif
!
                        call momdls( fparm(1,isi),iparm(1,isi),icalc,
     #                       alpha(ixt),beta(ixt),ibasis(1,ixbp),ixb,
     #                       spectr(ixs,isi),wspec,nft(ise),
     #                       ltd(isi,ise),ierr )
!
                        if (hltchk(ierr,isi,ise,iflag)) return
                     end do
                  end if
!
!                ----------------------------------------------------------
!                 For the current experimental spectrum calculate residuals
!                 from shifted spectra and scale factors (unless sites 
!                 scales are to be determined globally)
!                ----------------------------------------------------------
                  if (.not. glbscal) then
                     do j=1,npts(ise)
                        k=ixs+j-1
                        fvec(k)=data(k)
                        do isi=1,nsite
                           fvec(k)=fvec(k)-sfac(isi,ise)*spectr(k,isi)
                        end do
                     end do
                  end if
!     
!             ---------------------------------
!              No shifting: use sscale routine
!             ---------------------------------
               else
!
!                -----------------------------------------------------------
!                 Find least-squares site scaling factors and return
!                 residuals in fvec; however, multiple site/multiple spectra
!                 problems must be scaled globally, outside loop over spectra
!                -----------------------------------------------------------
                  if (.not. glbscal)
     #                call sscale( data(ixs),spectr(ixs,1),wspec,MXPT,
     #                             work,nsite,npts(ise),ctol,noneg,
     #                             iscal,sfac(1,ise),fvec(ixs) )
!
!                      *** End of "if (shifting allowed)...else..."
               end if
!
!                   *** End of "if (ise.le.nspc)"
            else
               do j=1,npts(ise)
                  k=ixs+j-1
                  fvec(k)=ZERO
                  do isi=1,nsite
                     fvec(k)=fvec(k)-sfac(isi,ise)*spectr(k,isi)
                  end do
               end do

            end if

            ixs=ixs+npts(ise)

!                  *** End of loop over spectra in a series
         end do
!
!        -------------------------------------------------------------- 
!         For multiple-site multiple spectral fits, the scaling factors 
!         are determined by a global least-squares fit to all spectra 
!         in the series
!        -------------------------------------------------------------- 
         if (glbscal) then
!
            call sscale( data,spectr,wspec,MXPT,work,nsite,m,ctol,
     #                   noneg,iscal,sfac,fvec )
!
!          --------------------------------------------
!           Copy the global scale factors to all spectra
!          --------------------------------------------
            do ise=2,nser
               do isi=1,nsite
                  sfac(isi,ise)=sfac(isi,1)
               end do
            end do
!
         end if
!
!       ------------------------------
!        Plot the current calculation
!       ------------------------------
         do i=1,nser

            call fstplt( data(ixsp(i)),fvec(ixsp(i)), 
     #                   sbi(i)-shft(i)-tmpshft(i),sdb(i),npts(i),i )
         end do
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  New calling sequence when pltd is implemented:
!
!         do i=1,nser
!            ixs=ixsp(i)
!            call pltwin( data(ixs),fvec(ixs),spectr(ixs,1),sfac(1,i),
!     *                   MXPT,npts(i),nsite,i)
!         end do
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!     ----------------------------------------------------------
!     For weighted residual fit, normalize differences in fvec
!     to the spectral noise (then |fvec|^2 will approximate 
!     the true chi^2)
!     ----------------------------------------------------------
         if (iwflag.ne.0) then
            do isp=1,nspc
               ixs=ixsp(isp)
               do j=1,npts(isp)
                  k=ixs+j-1
                  fvec(k)=fvec(k)/rmsn(isp)
               end do
            end do
         end if
!
         if (xreset) call xpack(x,n)
!
!**********************************************************************
!**********************************************************************
! IFLAG = 2: Jacobian evaluation by forward-difference approximation:
!
!             fjac = (f(x0+dx)-f(x0))/dx
!
! Note that f(x) = data-calc(x), so that f(x0+dx)-f(x0) = calc(x0)-calc(x0+dx)
!
!**********************************************************************
!**********************************************************************
!
      else if (iflag.eq.2) then
!
!       ------------------------------------------------
!        Update shifts at the beginning of an iteration
!       ------------------------------------------------
         do isp=1,nspc
            shft(isp)=shft(isp)+tmpshft(isp)
            tmpshft(isp)=ZERO
         end do
!
!        -------------------------------------------------------------
!         Loop over all parameters, introducing the forward-difference
!         step into each parameter
!        -------------------------------------------------------------
         do ix=1,n
            xtemp=x(ix)
            x(ix)=x(ix)+xfdstp(ix)
            call setprm(ixpr(ix),ixst(ix),x(ix))
!
!         ------------------------------------
!          Loop over all spectra in a series
!         ------------------------------------
            do isp=1,nspc
               ixs=ixsp(isp)
!
!             ---------------------
!              Initialize Jacobian
!             ---------------------
               do j=1,npts(isp)
                  k=ixs+j-1
                  fjac(k,ix)=ZERO
               end do
               covarOK=.false.
!
!             ----------------------------------------
!              Loop over all sites in the spectrum
!             ----------------------------------------
               do isi=1,nsite
!
!              ---------------------------------------------------------
!               Check whether this site/spectrum depends upon the
!               given parameter before calculating the partial derivative
!              ---------------------------------------------------------
                  sppar=spcpar( ixpr(ix) )
                  if ( (abs(sfac(isi,isp)).gt.RNDOFF) .and.
     #                (ixst(ix).le.0 .or.
     #                (ixst(ix).eq.isp .and. sppar) .or.
     #                (ixst(ix).eq.isi .and. .not. sppar) ) ) then
!
                     icalc=CFONLY
                     if (modtd(isi,isp).ne.0) icalc=FULL
!
                     call setspc(isi,isp)
                     if (hltchk(ierr,isi,isp,iflag)) return
!
                     ixt=ixtd(isi,isp)
                     if(basno(isi,isp).gt.0) then
                       ixb=ixbas( basno(isi,isp) )          
                       ixbp=ixb
                     else
                       ixb=0
                       ixbp=1
                     endif
!
                     call momdls( fparm(1,isi),iparm(1,isi),icalc,
     #                    alpha(ixt),beta(ixt),ibasis(1,ixbp),ixb,
     #                    work,wspec,nft(isp),ltd(isi,isp),ierr )
!
                     if (hltchk(ierr,isi,isp,iflag)) return
!
!                    ------------------------------------------
!                     Add in difference spectrum for this site
!                    ------------------------------------------
                     if (iwflag.ne.0) then
!                                               --- weighted ---
                        do j=1,npts(isp)
                           k=ixs+j-1
                           fjac(k,ix)=fjac(k,ix)+(spectr(k,isi)-work(j))
     #                          *sfac(isi,isp)/rmsn(isp)
                        end do
!                                               --- unweighted ---
                     else
                        do j=1,npts(isp)
                           k=ixs+j-1
                           fjac(k,ix)=fjac(k,ix)+(spectr(k,isi)-work(j))
     #                          *sfac(isi,isp)
                        end do
!
                     end if
!
                  endif
!
!                     *** end of loop over sites in spectrum 
               end do
!
!                  *** end of loop over spectra in series
            end do
!
!           -------------------------------------------------------  
!            Use the step size to calculate the partial derivative 
!            w.r.t. x(ix) for the current spectrum
!           -------------------------------------------------------  
            do j=1,m
               fjac(j,ix)=fjac(j,ix)/xfdstp(ix)
            end do
!
            x(ix)=xtemp
            call setprm(ixpr(ix),ixst(ix),x(ix))
!
!                 *** end of loop over parameters
         end do
!
!        ------------------------------------------------------------
!        Place unscaled spectra in last <nsite> columns of Jacobian
!        matrix. (These are the partial derivative of the spectrum
!        with respect to each of the scaling factors)
!        These columns will not be allowed to pivot, but will be used
!        to calculate the covariance matrix of the scale factors.
!        ------------------------------------------------------------
         if (nsite.gt.1) then
            do isi=1,nsite
               nj=n+isi
               do isp=1,nspc
                  ixs=ixsp(isp)
                  do j=1,npts(isp)
                     k=ixs+j-1
                     fjac(k,nj)=spectr(k,isi)
                     if (iwflag.ne.0) fjac(k,nj)=fjac(k,nj)/rmsn(isp)
                  end do
               end do
            end do
         else
            do isp=1,nspc
               nj=n+isp
               lastsp=ixsp(isp)+npts(isp)-1
               do j=1,ixsp(isp)-1
                  fjac(j,nj)=ZERO
               end do
               do j=ixsp(isp),lastsp
                  fjac(j,nj)=spectr(j,1)
                  if (iwflag.ne.0) fjac(j,nj)=fjac(j,nj)/rmsn(isp)
               end do
               do j=lastsp+1,m
                  fjac(j,nj)=ZERO
               end do
            end do
         end if
               
!
!**********************************************************************
!**********************************************************************
!     IFLAG=3
!     Fit has terminated: print result and save x-vector 
!**********************************************************************
!**********************************************************************
!
      else if (abs(iflag).eq.3) then

!        --------------------------------------------------------
!        Print final parameters if termination is by convergence
!        --------------------------------------------------------
         if (iflag.lt.0) then

            if (iwflag.ne.0) then
               rdchsq=fnorm*fnorm/float(m-n)
               write (luout,1009) 'T',iter,rdchsq,(x(i),i=1,n)
               if (luout.ne.luttyo) 
     #           write (luttyo,1009) 'T',iter,rdchsq,(x(i),i=1,n)
            else
               write (luout,1009) 'T',iter,fnorm,(x(i),i=1,n)
               if (luout.ne.luttyo) 
     #           write (luttyo,1009) 'T',iter,fnorm,(x(i),i=1,n)
            end if

         end if
!
!        ----------------------------------------
!         Store parameters currently in x-vector
!        ----------------------------------------
         do i=1,n
            call setprm(ixpr(i),ixst(i),x(i))
         end do
!
         if (ld.gt.0) then
           write (luttyo,1008) dashes(:ld)
           if (luout.ne.luttyo) write (luout,1008) dashes(:ld)
         end if
      end if

!
      return
!
! ### format statements ##############################################
!
 1006 format(/a)
 1007 format('  Itr  RedChSq  ',10(1x,a9))
 1008 format(a)
 1009 format(a1,i3,2x,g10.5,1x,10(f10.5))
 1010 format('  Itr  RsdNorm  ',10(1x,a9))
      end



!----------------------------------------------------------------------
!                    =========================
!                         function HLTCHK
!                    =========================
!
!> @brief   Check whether a user halt (control-C) or other error 
!>    has occurred during the spectral calculation and set the
!>    iflag error flag for LMNLS accordingly. Returns .true.
!>    for fatal errors or user halts, and .false. otherwise.
!
!----------------------------------------------------------------------
      function hltchk(ierr,isite,ispec,iflag)
!
      use stdio
      use errmsg
!
      implicit none
      integer ierr,iflag,isite,ispec
      logical hltchk
!
      hltchk=.false.
      if (ierr.eq.0) return
!
!      ------------------------------------
!      Non-fatal errors: warn if requested    
!      ------------------------------------
      if (ierr.lt.FATAL) then
         if (warn) then
            write (luout,1000) 'Warning',isite,ispec,eprerr(ierr)
            if (luout.ne.luttyo) write (luttyo,1000) 'Warning',
     #                                   isite,ispec,eprerr(ierr)
            warn=.false.
         end if
!
!        ------------------------------
!        Fatal errors (non-user related)
!        ------------------------------
      else if (ierr.ge.FATAL .and. ierr.ne.MTXHLT
     #        .and. ierr.ne.CGHLT) then 
               
         write (luout,1000) 'Fatal err',isite,ispec,eprerr(ierr)
         if (luout.ne.luttyo) then
            write (luttyo,1000) 'Fatal err',isite,ispec,eprerr(ierr)
         end if
         hltchk=.true.
         iflag=-1
!
!        ------------------------------
!        Other (user halt)
!        ------------------------------
      else
         write (luout,1001) eprerr(ierr)
         if (luout.ne.luttyo) write (luttyo,1001) eprerr(ierr)
         hltchk=.true.
         iflag=-1
      end if
!
      return
!
!----------------------------------------------------------------------
!
 1000 format(/2x,'*** ',a,' site',i2,' spctrm',i2,': ',a)
 1001 format(/2x,'*** ',a)
      end

