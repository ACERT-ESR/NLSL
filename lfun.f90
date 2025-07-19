c NLSL Version 1.5 beta 11/23/95
c----------------------------------------------------------------------
c> @file 
c> This seems to be the central hub, which is called by the L-M
c> optimization, and which calls the functions that generate the
c> simulated spectrum
c                         =======================
c                            subroutine LFUN
c                         =======================
c
c> @brief Subroutine for interfacing EPRLL spectral calculations with the MINPACK 
c>  version of the Levenberg-Marquardt nonlinear least squares algorithm 
c>  for fitting experimental spectra. The function performed
c>  function determined by the IFLAG argument. 
c>  MOMDLS subroutine to calculate a spectrum or series of spectra, each
c>  of which consists of one or more spectral components, or "sites". 
c> 
c> @param iflag
c> @parblock
c> > IFLAG=0 (Printing of iterates)
c> >  
c> > IFLAG=1 (Evaluation of residuals function)
c> > 
c> >   For a function evaluation, this routine takes the input vector X
c> >   of search parameters, loads them into the fparm parameter arrays in 
c> >   /parcom/ (using additional indices available in /parcom/) and then 
c> >   performs spectral calculations by calling the MOMDLS subroutine.
c> > 
c> >   The experimental data to which the calculations are to be compared
c> >   is contained in the common block /expdat/. LFUN optionally shifts
c> >   the calculated spectra to minimize the squared differences between
c> >   calculation and data, and then automatically determines a least-
c> >   squares scale factor (or set of scale factors for multiple sites/spectra)
c> >   which are returned in the vector sfac in common /mspctr/. The 
c> >   The unscaled spectra are returned in the matrix spectr in common
c> >   /mspctr/. Each column of spectr corresponds to a single site;
c> >   multiple spectra are stored sequentially in the columns of spectr.
c> > 
c> >   The residuals (differences between calculation and data after shifting
c> >   and scaling) are scaled by the standard deviation of the estimated 
c> >   noise in each experimental spectrum, and returned to MINPACK through 
c> >   the fvec argument.
c> > 
c> > IFLAG=2 (Jacobian evaluation)
c> > 
c> >   For evaluation of the Jacobian, the routine uses the forward-differences
c> >   approximation to calculate the partial derivative of the spectrum with
c> >   respect to each of the parameters being varied. The routine assumes
c> >   that the individual spectra calculated at the point in parameter space
c> >   where the Jacobian is being evaluated are contained in in the spectr
c> >   matrix, and the appropriate scaling factors in the sfac array (both
c> >   in common /mspctr/ ). The Jacobian matrix is returned in the fjac
c> >   argument.
c> > 
c> > IFLAG=3
c> >     Fit is terminating. Output most recent set of parameter values and
c> >     exit.
c> @endparblock
c Uses:
c    setspc
c    sshift
c    sscale
c    momdls
c    fstplt   (X-Windows interface)
c
c---------------------------------------------------------------------- 
      subroutine lfun( m,n,x,fvec,fjac,ldfjac,iflag )
c
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
c      use prmeqv
      use iterat
      use rnddbl
c
      implicit none
      integer m,n,iflag,ldfjac,ixbp
c
      double precision x(n),fvec(m),fjac(ldfjac,n+nsite+nspc)
      integer i,icalc,ierr,ise,isi,isp,ix,ixb,ixs,ixt,j,k,ld,lastsp,nj
      double precision fieldi,shift,snm,wdth,xtemp,dummy
      character dashes*132,tmpnam*10,shftnd*1
      logical shiftOK,glbscal,sppar
c
      integer FULL,CFONLY
      double precision ZERO
      parameter (FULL=111,CFONLY=0,ZERO=0.0D0)
c
      integer itrim
      double precision sshift,wtdres
      logical spcpar,hltchk
      external sshift,itrim,spcpar,hltchk,wtdres
c
c######################################################################
c
      ierr=0
      ld=0
      glbscal=nsite.gt.1 .and. nspc.gt.1
      if (iter.gt.1) warn=.false.
c
c**********************************************************************
c**********************************************************************
c     IFLAG=0: Output iteration information
c**********************************************************************
c**********************************************************************
c
c
c     ------------------------------------------
c     Draw a line of dashes for first iteration
c     ------------------------------------------
      if (iflag.eq.0) then
         if (iter.eq.1) then
            ld=12*n+16
            if (ld.gt.132) ld=132
            do i=1,ld
               dashes(i:i)='-'
            end do
c
c           -----------------------------------
c           Output names of current parameters 
c           -----------------------------------
c
            if (luout.ne.luttyo) then
               if (ld.gt.0)write (luttyo,1006) dashes(:ld)
               if (iwflag.ne.0) then
                  write (luttyo,1007) (tag(i)(:itrim(tag(i))),i=1,n)
               else
                  write (luttyo,1010) (tag(i)(:itrim(tag(i))),i=1,n)
               end if
               if (ld.gt.0)write (luttyo,1008) dashes(:ld)
            end if
c
            if (ld.gt.0)write (luout,1006)  dashes(:ld)
            if (iwflag.ne.0) then
               write (luout,1007) (tag(i)(:itrim(tag(i))),i=1,n)
            else
               write (luout,1010) (tag(i)(:itrim(tag(i))),i=1,n)
            end if
            if (ld.gt.0)write (luout,1008)  dashes(:ld)
         end if
c
         if (ishglb.ne.0 .and. (nshift.eq.0 .or.iter.le.nshift) ) then
            shftnd='S'
         else
            shftnd=' '
         end if
c
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
c
c        ------------------------------------------------------------
c        Output any other requested information regarding iterates
c        ------------------------------------------------------------
         if (iitrfl.ne.0) call wrfun(iter)
         if (itridg.ne.0) call wrtrdg(iter)
         if (jacobi.ne.0) call wrjac(iter)
c
c**********************************************************************
c**********************************************************************
c     IFLAG=1: Function evaluation
c**********************************************************************
c**********************************************************************
c
      else if (iflag.eq.1) then
c
         xreset=.false.
         written=0
c
         do i=1,n
            call setprm(ixpr(i),ixst(i),x(i))
         end do
c
         shiftOK = (nshift.eq.0 .or. iter.le.nshift)
         call tdsqz()
c
c       -----------------------------------
c        Loop over all spectra in a series
c       -----------------------------------
c
         ixs=1
         do ise=1,nser
            ixsp(ise)=ixs
            tmpshft(ise)=ZERO
c     
c          --------------------------------------
c           Loop over all sites for the spectrum
c          --------------------------------------
c
            do isi=1,nsite
               call tdchek(isi,ise,ierr)
               call setspc(isi,ise)
               if (hltchk(ierr,isi,ise,iflag)) return
c
               ixt=ixtd(isi,ise)
               if(basno(isi,ise).gt.0) then
                 ixb=ixbas( basno(isi,ise) )
                 ixbp=ixb
               else
                 ixb=0
                 ixbp=1
               endif
                 icalc=CFONLY
               if (modtd(isi,ise).ne.0) icalc=FULL
c
c              ---------------------------------
c              Calculate an individual spectrum
c              ---------------------------------
c
               call momdls( fparm(1,isi),iparm(1,isi),icalc,
     #              alpha(ixt),beta(ixt),ibasis(1,ixbp),ixb,
     #              spectr(ixs,isi),wspec,nft(ise),ltd(isi,ise),
     #              ierr )
c
               modtd(isi,ise)=0
               if (hltchk(ierr,isi,ise,iflag)) return
            end do
c
c
c           --------------------------------------------------------
c           If data are available for this calculation, determine
c           optimal shift and scaling and calculate fvec accordingly
c           (This feature permits calculations w/out data)
c           --------------------------------------------------------
            if (ise.le.nspc) then
c
c             ----------------------------------------------------------
c              If shifting is enabled, find the shift that will give 
c              optimal overlap between a calculated spectrum (or set of 
c              spectra) and the data for each element of a series
c             -----------------------------------------------------------
               if (shiftOK .and. ishft(ise).ne.0) then
c
                  shift=sshift( data(ixs),spectr(ixs,1),MXPT,nsite,
     #                  npts(ise),nft(ise),idrv(ise),srange,ctol,
     #                  noneg,tmpdat,tmpclc,wspec,work,sfac(1,ise) )
                  tmpshft(ise)=shift*sdb(ise)
c
c                 --------------------------------------------------------
c                  Re-calculate shifted spectrum with continued-fractions
c                  (this will keep the frwrd.-diff. derivative reasonable)
c                 --------------------------------------------------------
                  if (dabs(shift).gt.1.0d-3*sdb(ise)) then
                     do isi=1,nsite
                        icalc=CFONLY
                        call setspc(isi,ise)
                        if (hltchk(ierr,isi,ise,iflag)) return
c
                        ixt=ixtd(isi,ise)
                        if (basno(isi,ise).gt.0) then
                          ixb=ixbas( basno(isi,ise) )
                          ixbp=ixb
                        else
                          ixb=0
                          ixbp=1
                        endif
c
                        call momdls( fparm(1,isi),iparm(1,isi),icalc,
     #                       alpha(ixt),beta(ixt),ibasis(1,ixbp),ixb,
     #                       spectr(ixs,isi),wspec,nft(ise),
     #                       ltd(isi,ise),ierr )
c
                        if (hltchk(ierr,isi,ise,iflag)) return
                     end do
                  end if
c
c                ----------------------------------------------------------
c                 For the current experimental spectrum calculate residuals
c                 from shifted spectra and scale factors (unless sites 
c                 scales are to be determined globally)
c                ----------------------------------------------------------
                  if (.not. glbscal) then
                     do j=1,npts(ise)
                        k=ixs+j-1
                        fvec(k)=data(k)
                        do isi=1,nsite
                           fvec(k)=fvec(k)-sfac(isi,ise)*spectr(k,isi)
                        end do
                     end do
                  end if
c     
c             ---------------------------------
c              No shifting: use sscale routine
c             ---------------------------------
               else
c
c                -----------------------------------------------------------
c                 Find least-squares site scaling factors and return
c                 residuals in fvec; however, multiple site/multiple spectra
c                 problems must be scaled globally, outside loop over spectra
c                -----------------------------------------------------------
                  if (.not. glbscal)
     #                call sscale( data(ixs),spectr(ixs,1),wspec,MXPT,
     #                             work,nsite,npts(ise),ctol,noneg,
     #                             iscal,sfac(1,ise),fvec(ixs) )
c
c                      *** End of "if (shifting allowed)...else..."
               end if
c
c                   *** End of "if (ise.le.nspc)"
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

c                  *** End of loop over spectra in a series
         end do
c
c        -------------------------------------------------------------- 
c         For multiple-site multiple spectral fits, the scaling factors 
c         are determined by a global least-squares fit to all spectra 
c         in the series
c        -------------------------------------------------------------- 
         if (glbscal) then
c
            call sscale( data,spectr,wspec,MXPT,work,nsite,m,ctol,
     #                   noneg,iscal,sfac,fvec )
c
c          --------------------------------------------
c           Copy the global scale factors to all spectra
c          --------------------------------------------
            do ise=2,nser
               do isi=1,nsite
                  sfac(isi,ise)=sfac(isi,1)
               end do
            end do
c
         end if
c
c       ------------------------------
c        Plot the current calculation
c       ------------------------------
         do i=1,nser

            call fstplt( data(ixsp(i)),fvec(ixsp(i)), 
     #                   sbi(i)-shft(i)-tmpshft(i),sdb(i),npts(i),i )
         end do
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c  New calling sequence when pltd is implemented:
c
c         do i=1,nser
c            ixs=ixsp(i)
c            call pltwin( data(ixs),fvec(ixs),spectr(ixs,1),sfac(1,i),
c     *                   MXPT,npts(i),nsite,i)
c         end do
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c     ----------------------------------------------------------
c     For weighted residual fit, normalize differences in fvec
c     to the spectral noise (then |fvec|^2 will approximate 
c     the true chi^2)
c     ----------------------------------------------------------
         if (iwflag.ne.0) then
            do isp=1,nspc
               ixs=ixsp(isp)
               do j=1,npts(isp)
                  k=ixs+j-1
                  fvec(k)=fvec(k)/rmsn(isp)
               end do
            end do
         end if
c
         if (xreset) call xpack(x,n)
c
c**********************************************************************
c**********************************************************************
c IFLAG = 2: Jacobian evaluation by forward-difference approximation:
c
c             fjac = (f(x0+dx)-f(x0))/dx
c
c Note that f(x) = data-calc(x), so that f(x0+dx)-f(x0) = calc(x0)-calc(x0+dx)
c
c**********************************************************************
c**********************************************************************
c
      else if (iflag.eq.2) then
c
c       ------------------------------------------------
c        Update shifts at the beginning of an iteration
c       ------------------------------------------------
         do isp=1,nspc
            shft(isp)=shft(isp)+tmpshft(isp)
            tmpshft(isp)=ZERO
         end do
c
c        -------------------------------------------------------------
c         Loop over all parameters, introducing the forward-difference
c         step into each parameter
c        -------------------------------------------------------------
         do ix=1,n
            xtemp=x(ix)
            x(ix)=x(ix)+xfdstp(ix)
            call setprm(ixpr(ix),ixst(ix),x(ix))
c
c         ------------------------------------
c          Loop over all spectra in a series
c         ------------------------------------
            do isp=1,nspc
               ixs=ixsp(isp)
c
c             ---------------------
c              Initialize Jacobian
c             ---------------------
               do j=1,npts(isp)
                  k=ixs+j-1
                  fjac(k,ix)=ZERO
               end do
               covarOK=.false.
c
c             ----------------------------------------
c              Loop over all sites in the spectrum
c             ----------------------------------------
               do isi=1,nsite
c
c              ---------------------------------------------------------
c               Check whether this site/spectrum depends upon the
c               given parameter before calculating the partial derivative
c              ---------------------------------------------------------
                  sppar=spcpar( ixpr(ix) )
                  if ( (abs(sfac(isi,isp)).gt.RNDOFF) .and.
     #                (ixst(ix).le.0 .or.
     #                (ixst(ix).eq.isp .and. sppar) .or.
     #                (ixst(ix).eq.isi .and. .not. sppar) ) ) then
c
                     icalc=CFONLY
                     if (modtd(isi,isp).ne.0) icalc=FULL
c
                     call setspc(isi,isp)
                     if (hltchk(ierr,isi,isp,iflag)) return
c
                     ixt=ixtd(isi,isp)
                     if(basno(isi,isp).gt.0) then
                       ixb=ixbas( basno(isi,isp) )          
                       ixbp=ixb
                     else
                       ixb=0
                       ixbp=1
                     endif
c
                     call momdls( fparm(1,isi),iparm(1,isi),icalc,
     #                    alpha(ixt),beta(ixt),ibasis(1,ixbp),ixb,
     #                    work,wspec,nft(isp),ltd(isi,isp),ierr )
c
                     if (hltchk(ierr,isi,isp,iflag)) return
c
c                    ------------------------------------------
c                     Add in difference spectrum for this site
c                    ------------------------------------------
                     if (iwflag.ne.0) then
c                                               --- weighted ---
                        do j=1,npts(isp)
                           k=ixs+j-1
                           fjac(k,ix)=fjac(k,ix)+(spectr(k,isi)-work(j))
     #                          *sfac(isi,isp)/rmsn(isp)
                        end do
c                                               --- unweighted ---
                     else
                        do j=1,npts(isp)
                           k=ixs+j-1
                           fjac(k,ix)=fjac(k,ix)+(spectr(k,isi)-work(j))
     #                          *sfac(isi,isp)
                        end do
c
                     end if
c
                  endif
c
c                     *** end of loop over sites in spectrum 
               end do
c
c                  *** end of loop over spectra in series
            end do
c
c           -------------------------------------------------------  
c            Use the step size to calculate the partial derivative 
c            w.r.t. x(ix) for the current spectrum
c           -------------------------------------------------------  
            do j=1,m
               fjac(j,ix)=fjac(j,ix)/xfdstp(ix)
            end do
c
            x(ix)=xtemp
            call setprm(ixpr(ix),ixst(ix),x(ix))
c
c                 *** end of loop over parameters
         end do
c
c        ------------------------------------------------------------
c        Place unscaled spectra in last <nsite> columns of Jacobian
c        matrix. (These are the partial derivative of the spectrum
c        with respect to each of the scaling factors)
c        These columns will not be allowed to pivot, but will be used
c        to calculate the covariance matrix of the scale factors.
c        ------------------------------------------------------------
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
               
c
c**********************************************************************
c**********************************************************************
c     IFLAG=3
c     Fit has terminated: print result and save x-vector 
c**********************************************************************
c**********************************************************************
c
      else if (abs(iflag).eq.3) then

c        --------------------------------------------------------
c        Print final parameters if termination is by convergence
c        --------------------------------------------------------
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
c
c        ----------------------------------------
c         Store parameters currently in x-vector
c        ----------------------------------------
         do i=1,n
            call setprm(ixpr(i),ixst(i),x(i))
         end do
c
         if (ld.gt.0) then
           write (luttyo,1008) dashes(:ld)
           if (luout.ne.luttyo) write (luout,1008) dashes(:ld)
         end if
      end if

c
      return
c
c ### format statements ##############################################
c
 1006 format(/a)
 1007 format('  Itr  RedChSq  ',10(1x,a9))
 1008 format(a)
 1009 format(a1,i3,2x,g10.5,1x,10(f10.5))
 1010 format('  Itr  RsdNorm  ',10(1x,a9))
      end



c----------------------------------------------------------------------
c                    =========================
c                         function HLTCHK
c                    =========================
c
c> @brief   Check whether a user halt (control-C) or other error 
c>    has occurred during the spectral calculation and set the
c>    iflag error flag for LMNLS accordingly. Returns .true.
c>    for fatal errors or user halts, and .false. otherwise.
c
c----------------------------------------------------------------------
      function hltchk(ierr,isite,ispec,iflag)
c
      use stdio
      use errmsg
c
      implicit none
      integer ierr,iflag,isite,ispec
      logical hltchk
c
      hltchk=.false.
      if (ierr.eq.0) return
c
c      ------------------------------------
c      Non-fatal errors: warn if requested    
c      ------------------------------------
      if (ierr.lt.FATAL) then
         if (warn) then
            write (luout,1000) 'Warning',isite,ispec,eprerr(ierr)
            if (luout.ne.luttyo) write (luttyo,1000) 'Warning',
     #                                   isite,ispec,eprerr(ierr)
            warn=.false.
         end if
c
c        ------------------------------
c        Fatal errors (non-user related)
c        ------------------------------
      else if (ierr.ge.FATAL .and. ierr.ne.MTXHLT
     #        .and. ierr.ne.CGHLT) then 
               
         write (luout,1000) 'Fatal err',isite,ispec,eprerr(ierr)
         if (luout.ne.luttyo) then
            write (luttyo,1000) 'Fatal err',isite,ispec,eprerr(ierr)
         end if
         hltchk=.true.
         iflag=-1
c
c        ------------------------------
c        Other (user halt)
c        ------------------------------
      else
         write (luout,1001) eprerr(ierr)
         if (luout.ne.luttyo) write (luttyo,1001) eprerr(ierr)
         hltchk=.true.
         iflag=-1
      end if
c
      return
c
c----------------------------------------------------------------------
c
 1000 format(/2x,'*** ',a,' site',i2,' spctrm',i2,': ',a)
 1001 format(/2x,'*** ',a)
      end

