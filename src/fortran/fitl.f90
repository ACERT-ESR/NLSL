! NLSL Version 1.5 beta 11/23/95
!----------------------------------------------------------------------
!                    =========================
!                         program FITL
!                    =========================
!
! Carry out nonlinear least-squares analysis of slow-motional spectra
! the EPRLL calculation. Uses a modification of the Levenberg-Marquardt
! trust-region algorithm available in MINPACK. 
!
!  Uses
!    getdat
!    setfil
!    setdat
!    xpack (coded below)
!    wtdres (coded below)
!----------------------------------------------------------------------
!
      subroutine fitl
!
      use nlsdim
      use nlsnam
      use eprprm
      use expdat
      use parcom
      use mspctr
      use lpnam
      use lmcom
      use errmsg
      use iterat
      use stdio
      use timer
      use rnddbl
!
      implicit none
      integer i,ifile,iflag,iret,ix,nfev,njev
      real cpu
      logical snglclc,convrg
!
      integer getdat,itrim
!      real dtime
      double precision enorm,corrl,residx,wtdres
      external getdat,lfun,enorm,itrim,corrl,residx
!      external dtime
!
!######################################################################
!
        snglclc = maxitr.lt.1 .or. maxev.le.1 
     #           .or. nprm.le.0 .or. nspc.le.0
!
!      ------------------------------
!       Open trace file (if desired)
!      ------------------------------
        if (itrace.ne.0) then
          open(lutrc,file=trname(:lthfnm),status='unknown',
     #         access='sequential',form='formatted')
          itrace=lutrc
        end if
!
        call catchc( hltfit )
        warn=.true.
        call xpack(x,nprm)
!
!      -----------------------------------------------------------
!       Call lfun once if only a single calculation is specified
!      -----------------------------------------------------------
        if (snglclc) then
           iflag=1
!          cpu=dtime(tarray)
           call lfun(ndatot,nprm,x,fvec,fjac,mxpt,iflag)
!          cpu=dtime(tarray)
           lmflag=1
           info=11
           if (itridg.ne.0) call wrtrdg(0)

!
           if (nspc.gt.0) then
              do i=1,nspc
                 shft(i)=shft(i)+tmpshft(i)
                 tmpshft(i)=0.0D0
              end do
!
              fnorm=enorm(ndatot,fvec)
!                                         --- weighted residual fit ---
              if (iwflag.ne.0) then
                 chisqr=fnorm*fnorm
                 rdchsq=chisqr/float(ndatot)
!                                         --- unweighted residual fit ---
              else
                 chisqr=wtdres( fvec,ndatot,nspc,ixsp,npts,rmsn )**2
                 rdchsq=chisqr/float(ndatot)
              end if
!
              write(luout,1046) cpu
              write(luout,2036) fnorm,chisqr,rdchsq,corrl(),residx()
!
              if (luout.ne.luttyo) then
                 write(luttyo,1046) cpu
                 write(luttyo,2036) fnorm,chisqr,rdchsq,corrl(),residx()
              end if
           else
              write (luout,2037)
              if (luout.ne.luttyo) write(luttyo,2037)
           end if
!
        else
!
!======================================================================
!     Call Marquardt-Levenberg nonlinear least squares routine
!======================================================================
!
           nprint=1
           if (nsite.gt.1) then
              njcol=nprm+nsite
           else
              njcol=nprm+nspc
           end if
!          cpu=dtime(tarray)
!     
           call lmnls( lfun,ndatot,nprm,x,fvec,fjac,MXPT,ftol,xtol,gtol,
     1          maxev,maxitr,diag,prscl,factor,nprint,info,nfev,njev,
     2          ipvt,qtf,gnvec,gradf,work1,work2,work3,work4 )
!
           lmflag=1
!          cpu=dtime(tarray)
!
!                                         *** MINPACK error return
           convrg = info.ne.0 .and. info .le.3
           if (.not. convrg) then
              write (luout,2034) minerr(info)
              write (luout,2037) nfev,cpu 
              if (luout.ne.luttyo) then
                 write (luttyo,2034) minerr(info)
                 write (luttyo,2037) nfev,cpu 
              end if
              if (itrace.ne.0)  write (lutrc,2034) minerr(info)
!
!                                         *** Normal return
           else
!
              fnorm=enorm(ndatot,fvec)
              if (iwflag.ne.0) then
                 chisqr=fnorm*fnorm
                 rdchsq=chisqr/float(ndatot-nprm)
              else
                 chisqr=wtdres( fvec,ndatot,nspc,ixsp,npts,rmsn )
                 rdchsq=chisqr/float(ndatot)
              end if
!
              write(luout,1000)
              write(luout,2035) minerr(info),nprm,ndatot
              write(luout,2037) nfev,cpu 
              write(luout,2036) fnorm,chisqr,rdchsq,corrl(),residx()
              write(luout,1000)
!
              if (luout.ne.luttyo) then
                 write(luttyo,2035) minerr(info),nprm,ndatot
                 write(luttyo,2037) nfev,cpu 
                 write(luttyo,2036) fnorm,chisqr,rdchsq,corrl(),
     #                              residx()
              end if
!
              if (itrace.ne.0) then
                 write(lutrc,2035) minerr(info),nprm,ndatot
                 write(lutrc,2037) nfev,cpu 
                 write(lutrc,2036) fnorm,chisqr,rdchsq,corrl(),
     #                             residx()
              end if
!
              call covrpt( luout )
              if (luout.ne.luttyo) call covrpt( luttyo )

!              *** end 'if info.eq.0.or.info.gt.3 ... else'
         end if
!
!              *** end 'if single calculation...else'
      end if
!
!     ------------------------------------------------------------------
!     Optionally output splined data and calculated spectrum/spectra
!     for each datafile
!     ------------------------------------------------------------------
!     ----------------------------------------
!     Report shifting and scaling information 
!     ----------------------------------------
      call sclstt(luout)
      if (luout.ne.luttyo) call sclstt(luttyo)
!
!     ----------------------
!     Close files and return
!     ----------------------
      if (itrace.ne.0) close (lutrc)
      call uncatchc( hltfit )
      return
!   
!# format statements ##################################################
!
 1000 format(/,2x,70('-'),/)
 1046 format(//10x,'Single calculation: CPU time= ',f9.3,' sec')
!
 2034 format(/2x,'MINPACK could not find a solution: ', a )
 2035 format(/9x,'MINPACK completed: ', a/
     #       12x,i2,' parameters;',i5,' data points') 
 2036 format(10x,'Residual norm= ',g13.6/
     #       10x,'Chi-squared=',g13.6,2x,'Reduced Chi-sq=',g13.6/
     #       10x,'Correlation= ',f8.5,7x,'Residual index=',f8.5/)
 2037 format(12x,i4,' function evals (CPU time: ',f9.3,' sec)')
 2038 format(10x,'No data files were specified')
      end
!
!----------------------------------------------------------------------
!                    =========================
!                       subroutine XPACK
!                    =========================
!
!   Pack parameters from fparm array (and spectral parameter arrays)
!   into x vector.
!   
!----------------------------------------------------------------------
      subroutine xpack(x,n)
!
      use nlsdim
      use eprprm
      use expdat
      use parcom
      use rnddbl
!
      implicit none
      integer i,ix2,n
      double precision x(n)
!
! --- loop over all parameters in x-vector
!
      do 10 i=1,n
!
!----------------------------------------------------------------------
!     For certain parameters which cannot vary from site to site
!     in the same spectrum, the index refers to a spectrum in the
!     series rather than a site within a spectrum.
!     These parameters include LB, PHASE, B0, and PSI.
!----------------------------------------------------------------------
         ix2=ixst(i)
         if (ix2.le.0) ix2=1
         if (ixpr(i).eq.LB) then
            x(i)=slb(ix2)
!
!   --- B0: Spectrometer field
!
         else if (ixpr(i).eq.IB0) then
            x(i)=sb0(ix2)
!
!  --- PHASE: Microwave phase
!
         else if (ixpr(i).eq.IPHASE) then
            x(i)=sphs(ix2)
!
!  --- PSI: Director tilt
!
         else if (ixpr(i).eq.IPSI) then
            x(i)=spsi(ix2)
!     
!----------------------------------------------------------------------
!     All other parameters: index refers to a site in the spectrum.
!----------------------------------------------------------------------
         else
            x(i)=fparm(ixpr(i),ix2)
         end if
         
 10   continue
!
      return
      end

!----------------------------------------------------------------------
!                    =========================
!                         function WTDRES
!                    =========================
!
!     Returns weighted residuals in fvec (weighted according to
!     rms noise for each spectrum given in rmsn array)
!----------------------------------------------------------------------

      function wtdres( fvec, m, nspc, ixsp, npts, rmsn )
      implicit none
      integer m,nspc,ixsp(nspc),npts(nspc)
      double precision fvec(m),rmsn(nspc),wtdres
!
      integer isp,ixs,j,k
!
      wtdres=0.0d0
      do isp=1,nspc
         ixs=ixsp(isp)
         do j=1,npts(isp)
            k=ixs+j-1
            wtdres=wtdres+(fvec(k)/rmsn(isp))**2
         end do
      end do
      return
      end


