
c  VERSION 1.0  (NLSPMC version)   2/5/99
c----------------------------------------------------------------------
c                         =======================
c                            subroutine PFUN
c                         =======================
c
c Subroutine for interfacing EPRCGF spectral calculations with the
c MINPACK version of the Levenberg-Marquardt nonlinear least squares
c algorithm for fitting experimental spectra. This routine loads the
c necessary parameters into common /eprprm/ and then calls the SIM2D
c subroutine to calculate a spectrum or evaluate a Jacobian as required
c by the calling program.
c
c *** NLSPMC only ***
c For evaluation of a series of spectra (which share the eigenvalues and
c eigenvectors), the routine re-evaluates the spectrum using the saved
c eigenvalues and vectors.  Therefore any attempt to fit a series of 
c parameters that affects the eigenvalues is forbidden (e.g. series in
c tilt angle or temperature).
c
c For evaluation of the Jacobian, the routine uses the forward-differences
c approximation to calculate the partial derivative of the spectrum with
c respect to each of the parameters being varied. The routine assumes
c that the function values at the point where the Jacobian is being
c evaluated are contained in fvec and the scaling factor for each spectrum
c in series option is contained in sfac ( common /expdat/ )
c
c In the Jacobian evaluation, the routine makes a special check for 
c variation of the parameters relevent in series option (now gib,lib,
c hwid,mwid).  If line-broadening parameter (gib) is varied, it simply
c re-evaluates the spectrum using the saved time domain data (xspec).
c
c NOTE that in order to follow this strategy, the Jacobian should
c be calculated with respect to any line-broadening parameters before 
c any other parameter is varied (which would change the time domain
c signal in the absence of Gaussian inhomogeneous broadening).
c
c In the future, one might hope to recognize when only the real or 
c imaginary part of the matrix changes, and ultimately to utilize a
c perturbation approach to calculating the Jacobian.
c
c Modified calculation of fvec and fjac to include weighting option for
c total data sets.
c
c iflag=0 just print variables
c iflag=1 function evaluation
c iflag=2 or 3 Jacobian
c iflag=3 output Jacobian matrix
c
c Includes:
c    nlsdim.inc
c    stdio.inc
c    eprprm.inc
c    expdat.inc
c    parcom.inc
c    tdspec.inc
c    wkspcm.inc
c    iterat.inc
c
c Uses
c    sim2d
c    xshft
c    sscale
c---------------------------------------------------------------------- 
c
      subroutine pfun( totpts,n,xxx,fvec,fjac,ldfjac,iflag )
      implicit none
      logical vchange
      double precision enorm
      external vchange,enorm
      integer totpts,n,iflag,ldfjac,ispec,isite
      double precision xxx(n),fvec(totpts),fjac(ldfjac,n),xxxo(n)
c
      include 'limits.inc'
      include 'stdio.inc'
      include 'lpnam.inc'
      include 'simparm.inc'
      include 'datas.inc'
      include 'parms.inc'
      include 'tdspec.inc'
      include 'wkspcm.inc'
      include 'miscel.inc'
c      include 'iterat.inc'
c for storing the individual simulations:
      character*20 cchar
      data cchar/'1234567890abcdefghij'/
      character*(WORDLG) fnmtst/'simss??.tst'/
      real*8 amin,amax
      logical sitewgc
c
      integer i,icalc,ierr,npt1,npt2,isp,ixi,ixs,ixscw,j,k,ld
      double precision rmsdv,shift,resol,ovlap,xtemp,delx,ccc
      character dashes*132,jacnam*10
      integer jspec,jsite,ixprmx(MXVAR),iparmptr,ivar
c
      double precision zero,totwgt,siteweight,totwgt2(MXSPEC)
      parameter (zero=0.0D0)
c
c######################################################################
c
c----------------------------------------------------------------------
c    Print parameter values for the current iteration
c----------------------------------------------------------------------
c
      if (iflag.eq.0) then
         ld=min(132,12*n+16)
         do 1 i=1,ld
            dashes(i:i)='-'
    1       continue
         if (iter.eq.1) then
            if (luout.ne.luttyo) then
               write (luttyo,1006) dashes(:ld)
               write (luttyo,1007) (tag(i),i=1,n)
               write (luttyo,1008) dashes(:ld)
            end if
            write (luout,1006)  dashes(:ld)
            write (luout,1007)  (tag(i),i=1,n)
            write (luout,1008)  dashes(:ld)
         end if
         rmsdv=fnorm/dsqrt(dfloat(nptot))
         sitewgc=.false.
         do 2 i=1,n
           iparmptr=ixpr(i)
           if (iparmptr.eq.ISWGT) then
             sitewgc=.true.
             totwgt=0.0d0
             do 3 j=1,ncomps
               totwgt=totwgt+fparm(ISWGT,1,j)
 3           continue
             xxxo(i)=xxx(i)/totwgt
            else
             xxxo(i)=xxx(i) 
            end if
 2        continue
cc
         if (sitewgc) then
          if (luout.ne.luttyo)
     #        write (luttyo,1009) iter,rmsdv,(xxxo(i),i=1,n)
           write (luout,1009) iter,rmsdv,(xxxo(i),i=1,n)
         else
          if (luout.ne.luttyo)
     #        write (luttyo,1009) iter,rmsdv,(xxx(i),i=1,n)
           write (luout,1009) iter,rmsdv,(xxx(i),i=1,n)
         end if
cc
cc         if (luout.ne.luttyo)
cc     #                write (luttyo,1009) iter,rmsdv,(xxx(i),i=1,n)
cc         write (luout,1009) iter,rmsdv,(xxx(i),i=1,n)
         if (ihltcmd.ne.0) return
c
c**********************************************************************
c**********************************************************************
c IFLAG=1: Function evaluation
c**********************************************************************
c**********************************************************************
c
      else if (iflag.eq.1) then
c
c Here, the variables contained in xxx determine the new spectra to 
c be calculated.  If there is nothing changed, we don't repeat the 
c spectral calculation, if something is changed, we calculate only 
c the necessary parts.  We require the same parameters to apply 
c to all spectra. 
c 
         do 41 j=1,n	! for eack xxx variable,
c           ixprmx(j)=0	! keep track of variable severity
           iparmptr=ixpr(j)	! get variable id set in addprm, srchc
c for all sites,spectra affected
           do 42 isite=1,ncomps
             do 43 ispec=1,nspectra
c if this variable applies to this site/spectrum:
               if (ixx(iparmptr,ispec,isite).eq.j) then
c                 if (ispec.eq.1)write(*,*)'affecting site ',isite
c check if this parameter is changed
                 if (vchange(xxx(j),fparm(iparmptr,ispec,isite))) then
c if so, update fparm and basinfo to define simulation
                   matrxok(isite)=min(matrxok(isite),specid(iparmptr))
c                   if (ispec.eq.1)write(*,*)'changed to ',xxx(j)
                   fparm(iparmptr,ispec,isite)=xxx(j)
                 end if
               end if
               if (ixx(iparmptr,ispec,isite).ne.ixx(iparmptr,1,isite))
     #			then
                 stop
               end if
 43          continue
 42        continue
 41      continue
c
c Now matrxok tell if need new matrix etc.  After first matrix 
c calculation, update this information if that matrix applies to more
c than one site.  Must apply to all sites in this version
c     ...Loop over all spectra in a series and all sites....
c
c       call ftest2(fparm(1,1,1),NFPRM,1,'pfun,fparm 1c   ',14)
c       call ftest2(fparm(1,1,2),NFPRM,1,'pfun,fparm 2c   ',14)

         do 5 jsite=1,ncomps
          cursite=jsite	   ! used to keep track of basis, matrix etc.
         do 5 jspec=1,nspectra
c
c this is checked in datac:
c          if(ixsp(nser).gt.mxpt)then
c            write(*,*)'error in pfun, inc. mxpt',ixsp(nser),mxpt
c            stop
c          end if
          ixs=ixsp(jspec)	! pointer to this data set
          ixscw=ixspcw(jspec)	! this cw data set
c
c replace by new arrays recording what to do:
          icalc=2
          if(matrxok(jsite).ge.0.and.matrxok(jsite).le.2) then
	    icalc = 2-matrxok(jsite)  ! 2=need matrix,1=use eval
          elseif (matrxok(jsite).eq.4) then
c here if matrxok = 4
            icalc=0	! here if vary site weighting
          end if
c
c----------------------------------------------------------------------
c If desired, find the shift that will give optimal overlap between 
c a spectrum and the data.  This is achieved by comparing cw magnitude
c spectrum with the cw-equivalent spectrum extracted from 2D data,
c which was obtain in data command and saved in cwdata array in common
c block /expdat/.
c if got dat and want shift, must do cw calculation.  Here for multiple
c components, do for each time and sum to give composite cw spec.
c In the case of varying site weight, redo shift, but don't do it for 
c case of vary gib or lib.
c----------------------------------------------------------------------
c  Don't redo shift if varying gib.
c  For site weighting, don't call sim2d, to preserve xspec, but check 
c  shift.
c
          if (dataOK.and.(icalc.ne.0).and.(matrxok(jsite).ne.4).and.
     #	     (sishft(jspec).ne.0)) then
c set this spectrum parameters                  ** cw spectrum
             call setspc( jspec,jsite,npt1,npt2 )
c specify cw calculation
             fparm(ISHFT,jspec,jsite)=0.0d0
             iparm(IIEXP,jspec,jsite)=0
             iparm(IICOMB,jspec,jsite)=1
             resol=1.0d3/(sstept2(jspec)*float(npt2-1))
             npt1=npt2
             npt2=1	! specify CW calculation, result to 
c
             if(matrxok(jsite).ne.4) then
c need to mod sim2d to use site basis
                call sim2d( icalc,xspectr(ixs,jsite),cspectr(ixs,jsite),
     #                     npt1,npt2,n,iflag,ixi,delx,jspec,jsite,ierr )
c
               if(icalc.lt.1)stop 'icalc=0 error in pfunnew'
c
c return here means have eval, evec, and matrix for this spec,
c reset matrxok to indicate  this.
c note 2D spec requires calculation of xspec if not present.
c
               matrxok(jsite)=max(1,matrxok(jsite))  ! have matrix for site.
               if ((ihltcmd.eq.0).and. .false.) then
                 if (luout.ne.luttyo) write(luttyo,1004)
                 write(luout,1004)
                 iflag=-1
                 return
               else if (ierr.ne.0) then
                 if (luout.ne.luttyo) write(luttyo,1003)
                 write(luout,1003)
                 iflag=-1
               end if
             end if
c
c save the cw spectra for all sites,spectra. Later, sum with weights
c  do this so that if changing only weights, we don't recalculate cwsim
c
             do 6 j=1,npt1
               ccwsim(ixscw-1+j,jsite)=cspectr(ixs-1+j,jsite)
 6           continue
          end if	! end CW calculation.
 5       continue
c
c----------------------------------------------------------------------
c        Calculate shift term (shft).
c        First, sum over sites with site weighting
c----------------------------------------------------------------------
c
c we use absolute values of the weights normalized:
c
        do 71 jspec=1,nspectra
          if (dataOK.and.(icalc.ne.0).and.(sishft(jspec).ne.0)) then
             ixscw=ixspcw(jspec)	! this cw data set
             totwgt=0.0D0
            do 72 jsite=1,ncomps
              siteweight=abs(fparm(ISWGT,jspec,jsite))
              totwgt=totwgt+siteweight
              do 8 j=1,npt1
                if(jsite.eq.1) ccwsimtot(ixscw-1+j)=(0.0D0,0.0D0)
                ccwsimtot(ixscw-1+j)=ccwsimtot(ixscw-1+j)+
     #		 ccwsim(ixscw-1+j,jsite)*siteweight
                if(jsite.eq.ncomps) cwsimtot(ixscw-1+j)=
     #	         cdabs(ccwsimtot(ixscw-1+j))/totwgt
 8            continue
c
c  Then compare this spectrum with CW data to determine best shift
c
 72         continue
            call xshft(cwdata(ixscw),cwsimtot(ixscw),npt1,
     #               cwsp1,cwsp2,shift,ovlap)
            sshft(jspec)=resol*shift	! save shift for setspc
          end if
 71     continue
c
c----------------------------------------------------------------------
c        Calculate 2D spectrum with desired shift terms (sshft).
c----------------------------------------------------------------------
c
c       call ftest2(fparm(1,1,1),NFPRM,1,'pfun,fparm 1b   ',14)
c       call ftest2(fparm(1,1,2),NFPRM,1,'pfun,fparm 2b   ',14)
c
         do 11 jsite=1,ncomps
c          cursite=jsite	; used to keep track of basis, matrix etc.
         do 10 jspec=1,nspectra
          ixs=ixsp(jspec)	! pointer to this data set
c          ixscw=ixspcw(jspec)	! this cw data set
          icalc=2
          if(matrxok(jsite).ge.0.and.matrxok(jsite).le.2) then
	    icalc = 2-matrxok(jsite)  ! 2=need matrix,1=use eval, 0 for gib
          elseif (matrxok(jsite).eq.4) then
c here if matrxok = 4
            icalc=0	! here if vary site weighting
          end if
          call setspc( jspec,jsite,npt1,npt2 )
          if(matrxok(jsite).ne.4) then
c
c            time1=mclock()
            call sim2d( icalc,xspectr(ixs,jsite),cspectr(ixs,jsite),
     #                  npt1,npt2,n,iflag,ixi,delx,jspec,jsite,ierr )
c       call ftest3(xspectr(1,1),16384,'pfun,xsptr 1b   ',14)
c       call ftest3(xspectr(1,2),16384,'pfun,xsptr 2b   ',14)
            matrxok(jsite)=max(1,matrxok(jsite))  ! have matrix for site.
c
c note don't set matrxok = 0 since it won't be correct for next spectrum
c
            if ((ihltcmd.eq.0 .and. .false.)) then
               if (luout.ne.luttyo) write(luttyo,1004)
               write(luout,1004)
               iflag=-1
               return
            else if (ierr.ne.0) then
               if (luout.ne.luttyo) write(luttyo,1003)
               write(luout,1003)
               iflag=-1
               return
            end if
          end if
c
c               ... End of loop over spectra in a series and sites
c
 10      continue
          matrxok(jsite)=4	! got everything for this site.  If 
c                            nothing is changed, next calc will be fast
 11      continue
c
c Accumulate 2D spectra over sites, weighting appropriately
c We assume the site weights are identical for all spectra.
c This may not be true for temperature variation in the future.
c
         do 115 jspec=1,nspectra
           totwgt2(jspec)=0.0D0
           do 115 jsite=1,ncomps
             totwgt2(jspec)=totwgt2(jspec)+
     #		abs(fparm(ISWGT,jspec,jsite))
 115     continue
         do 12 jsite=1,ncomps
           i=0
           do 12 jspec=1,nspectra
c  for testing, write the site,spectrum simulations to files
	     fnmtst(6:6)=cchar(jspec:jspec)
	     fnmtst(7:7)=cchar(jsite:jsite)
             amin=0.0d0
             amax=0.0d0
             do 13 j=1,ndata(jspec)
               spectr(j,jsite)=cdabs(cspectr(j,jsite))
 13          continue
             call wrfit(spectr(1+(jspec-1)*ndata(jspec),jsite),
     #		jspec,amin,amax,fnmtst)
c
             siteweight=abs(fparm(ISWGT,jspec,jsite))
             do 12 j=1,ndata(jspec)
               i=i+1
               if(jsite.eq.1)ctotspec(i)=(0.0D0,0.0D0)
               ctotspec(i)=ctotspec(i)+cspectr(i,jsite)*siteweight
               if(jsite.eq.ncomps)totspec(i)=cdabs(ctotspec(i))/totwgt2(
     #			jspec)
c be careful, from this point on we work with totspec instead of spectr
 12      continue
c
c----------------------------------------------------------------------
c Calculate the scale factor(s) for each spectrum
c according to the dependencies (idepnd) that user specified in data
c command.  The subroutine returns the scaled spectra and the scale
c factors (sfac in common block /expdat/).
c----------------------------------------------------------------------
c
         if (dataOK) call sscale(totspec)
c
c----------------------------------------------------------------------
c Calculate the difference between the experimental and scaled
c calculated spectra and store the result in fvec.
c----------------------------------------------------------------------
c
           i=0
           do 20 jspec=1,nspectra
             do 20 j=1,ndata(jspec) 
               i=i+1
               if(dataOK) then
                 fvec(i)=(totspec(i)-data(i))*datwgt(jspec)
c               else
c                 fvec(i)=totspec(i)	! not used
               end if
 20        continue
c
         if (idebug.ne.0) then
            write (ludeb,2010)
            do 21 j=1,nser
               write (ludeb,2020) j,sfac(j),sshft(j)
 21         continue
         end if
c
c        ** Update X vector for the case where the real part of
c           the eigenvalues are negative.
c this is no longer done.  RC 11/11/98.
c         do 6 j=1,n
c          if(xxx(j)-fparm(ixpr(j)) .gt. 1.0D-10) then
c            if (luout.ne.luttyo) then
c              write(luttyo,*)'remapping x: ',j,xxx(j),fparm(ixpr(j))
c              write(luttyo,*)'This should not happen!'
c            end if
c              write(luout,*)'remapping x: ',j,xxx(j),fparm(ixpr(j))
c              write(luout,*)'This should not happen!'
c          end if
c            xxx(j)=fparm(ixpr(j))
c 6       continue
c
c**********************************************************************
c IFLAG >= 2: Jacobian evaluation by forward-difference approximation:
c             fjac = (f(x0+dx)-f(x0))/dx
c          where f(x0) = fvec+data, and dx is assumed to be positive
c We assume we have a previously calculated spectrum, for the current
c set of parameters.  Then each is varied by dx and a new spectrum
c calculated.  
c
c**********************************************************************
c
      else if (iflag.ge.2) then
c
c----------------------------------------------------------------------
c Loop over all parameters, introducing the forward-difference step
c into each parameter.  On entry, have spectra calculated
c from the current values of the parms and matrxok=4, save this info.
c As vary each parm, do only calculation necessary for that parameter.
c----------------------------------------------------------------------
c
         do 30 ivar=1,n	        ! for eack xxx variable, increment it
           xtemp=xxx(ivar)		! save current value of xxx(i)
           xxx(ivar)=xxx(ivar)+xfdstp(ivar)	! introduce step into this parm
           ixi=ixpr(ivar) 	! get variable id set in addprm, srchc
           delx=xfdstp(ivar)	! save step size
c
c for all sites,spectra affected, map only this xxx into fparm.
c  The others should not be changed.
c
           do 82 isite=1,ncomps
             do 83 ispec=1,nspectra
c if this variable applies to this site/spectrum:
               if (ixx(ixi,ispec,isite).eq.ivar) then
c update fparm to define simulation.  We know variable is changed.
                 matrxok(isite)=min(matrxok(isite),specid(ixi))
                 fparm(ixi,ispec,isite)=xxx(ivar)
               end if
               if (ixx(ixi,ispec,isite).ne.ixx(ixi,1,isite))
     #			 then
                 write(*,*)'different parameters for different '
                 write(*,*)'spectra not allowed yet (in pfun).'
                 stop
               end if
 83          continue
 82        continue
c 
c  for multiple variables, restore the previous one since we only
c  vary one parameter at a time in computing the Jacobian
c
           if (ivar.gt.1) then
             do 821 isite=1,ncomps
             do 831 ispec=1,nspectra
c if this variable applies to this site/spectrum:
               if (ixx(ixi,ispec,isite).eq.ivar-1) then
c update fparm to define simulation.  We know variable is changed.
                 matrxok(isite)=min(matrxok(isite),specid(ixi))
                 fparm(ixi,ispec,isite)=xxx(ivar-1)
               end if
 831         continue
 821         continue
           end if
c
c       call ftest2(fparm(1,1,1),NFPRM,1,'pfun,fparm 1a   ',14)
c       call ftest2(fparm(1,1,2),NFPRM,1,'pfun,fparm 2a   ',14)
c
c        ...Loop over all spectra in a series, same as do 10 loop
c
           do 28 jsite=1,ncomps
c            cursite=jsite        ; used to keep track of basis, matrix etc.
           do 27 jspec=1,nspectra
             ixs=ixsp(jspec)       ! pointer to this data set
c             ixscw=ixspcw(jspec)   ! this cw data set
             icalc=2
             if(matrxok(jsite).ge.0.and.matrxok(jsite).le.2) then
               icalc = 2-matrxok(jsite)  ! 2=need matrix,1=use eval
             elseif (matrxok(jsite).eq.4) then
c here if matrxok = 4
               icalc=0	! here if vary site weighting
             end if
             call setspc( jspec,jsite,npt1,npt2 )
             if(matrxok(jsite).ne.4) then
c               time1=mclock()
               call sim2d( icalc,xspectr(ixs,jsite),cspectr(ixs,jsite),
     #           npt1,npt2,n,iflag,ixi,delx,jspec,jsite,ierr )
c       call ftest3(xspectr(1,1),16384,'pfun,xsptr 1a   ',14)
c       call ftest3(xspectr(1,2),16384,'pfun,xsptr 2a   ',14)
               matrxok(jsite)=max(1,matrxok(jsite))  ! have matrix for site.
c
               if ((ihltcmd.ne.0 .and. .false. )) then
                 if (luout.ne.luttyo) write(luttyo,1004)
                 write(luout,1004)
                 iflag=-1
                 return
               else if (ierr.ne.0) then
                 if (luout.ne.luttyo) write(luttyo,1003)
                 write(luout,1003)
                 iflag=-1
                 return
               end if
             end if
c
c               ... End of loop over spectra in a series and sites
c
 27        continue
             matrxok(jsite)=4      ! got everything for this site.  If
c                            nothing is changed, next calc will be fast
 28        continue
c
c sum over the sites with weighting.  Shift already in spectra.  Use
c previous scale factor.
c
           do 215 jspec=1,nspectra
             totwgt2(jspec)=0.0D0
             do 215 jsite=1,ncomps
               totwgt2(jspec)=totwgt2(jspec)+
     #	         abs(fparm(ISWGT,jspec,jsite))
 215       continue

           do 92 jsite=1,ncomps
             i=0
             do 92 jspec=1,nspectra
               do 92 j=1,ndata(jspec)
                 i=i+1
                 if(jsite.eq.1)ctotspec(i)=(0.0D0,0.0D0)
                 ctotspec(i)=ctotspec(i)+cspectr(i,jsite)*
     %                      abs(fparm(ISWGT,jspec,jsite))
                 if(jsite.eq.ncomps)totspec(i)=
     #               cdabs(ctotspec(i))/totwgt2(jspec)
c be careful, from this point on we work with totspec instead of spectr
 92        continue
c           call ftest1(totspec,i,'92 pfunnew   ',11 )
c
c          ... Store function value before forward-differences step
c
           do 39 jspec=1,nspectra
             ixs=ixsp(jspec)
             do 25 j=1,ndata(jspec)
               k=ixs+j-1
c since data and fvec are not changed here, use to calc. prev simulation
c at the original values of xxx.
               fjac(k,ivar)=-(fvec(k)+data(k)*datwgt(jspec))
               fjac(k,ivar)=(fjac(k,ivar)+sfac(jspec)*totspec(k)*
     #	         datwgt(jspec))/delx
 25          continue
 39	   continue
c           call ftest2(fjac(1,ivar),mxpt,1,'39 pfunnew fjac  ',16 )
c
           xxx(ivar)=xtemp		! restore x
c
c Restore fparms array with this value, readying for next variable:
c
           do 820 isite=1,ncomps
             do 830 ispec=1,nspectra
c if this variable applies to this site/spectrum:
               if (ixx(ixi,ispec,isite).eq.ivar) then
c update fparm to define simulation.  We know variable is changed.
                 matrxok(isite)=min(matrxok(isite),specid(ixi))
                 fparm(ixi,ispec,isite)=xxx(ivar)
               end if
               if (ixx(ixi,ispec,isite).ne.ixx(ixi,1,isite))
     #			 then
                 write(*,*)'different parameters for different '
                 write(*,*)'spectra not allowed yet (in pfun).'
                 stop
               end if
 830          continue
 820        continue
c                     *** end of loop over parameters
 30      continue
c
c note here, the x values are as they were on entry to pfun, while
c  the spectra are calculated from the last parameter incremented.
c  This should not be a problem, since the fparm values do
c  correspond to the matrix.
c
c----------------------------------------------------------------------
c     Optional output of Jacobian matrix for this iteration
c----------------------------------------------------------------------  
c
         if (iflag.eq.3) then
            write(jacnam,2000) iter
            open(unit=ludisk,file=jacnam,status='unknown',
     #           access='sequential',form='formatted')
            do 40 i=1,totpts
               write (ludisk,2003) data(i),fvec(i),(fjac(i,j),j=1,n)
 40         continue
            close(ludisk)
         end if
c
      end if
c
      return
c
c ### format statements ##############################################
c
 1003 format(/20x,'*** Error: fit or search procedure terminating',
     #       ' ***'/)
 1004 format(/20x,'*** Fit halted by user ***')
 1006 format(/a)
 1007 format('Iter',5x,'RmsDv',10(6x,a))
 1008 format(a)
 1009 format(i4,4x,g10.4,2x,10(g12.6))
 2000 format('nlsjac.',i3.3)
 2003 format(10(g12.5,1x))
 2010 format(/' ** ',6x,'Spectrum',5x,'Scale Factor',5x,
     #       'Shift Factor **')
 2020 format(13x,i1,7x,g14.7,3x,g14.7)
c
      end

c----------------------------------------------------------------------
c                    =========================
c                       subroutine SETSPC
c                    =========================
c
c  Sets the series-dependent parameters in the fparm and iparm arrays.
c  These parameters are not allowed to change the eigenvalues.
c
c  Specifically, the parameters allowed in series option are
c
c    Integer :  iexp, icomb, npt1, npt2
c    Floating:  init1, stept1, init2, stept2, tfix, shft
c
c  on entry, ispc,isit specify the current spectrum,site.
c
c----------------------------------------------------------------------
c
      subroutine setspc( ispc,isit,npt1,npt2 )
      implicit none
      integer ispc,isit,npt1,npt2
      logical vchange
      external vchange
      logical chnge
c
      include 'limits.inc'
      include 'simparm.inc'
      include 'datas.inc'
      include 'parms.inc'
c
      if (ispc.le.0 .or. ispc.gt.nser) return
c
      chnge=.false.	! are we changing anything?
      chnge=chnge.or.vchange(fparm(IINIT1,ispc,isit),sinit1(ispc))
      fparm(IINIT1,ispc,isit) = sinit1(ispc)
      chnge=chnge.or.vchange(fparm(ISTEPT1,ispc,isit), sstept1(ispc))
      fparm(ISTEPT1,ispc,isit) = sstept1(ispc)
      chnge=chnge.or.vchange(fparm(IINIT2,ispc,isit), sinit2(ispc))
      fparm(IINIT2,ispc,isit) = sinit2(ispc)
      chnge=chnge.or.vchange(fparm(ISTEPT2,ispc,isit),sstept2(ispc))
      fparm(ISTEPT2,ispc,isit) = sstept2(ispc)
      chnge=chnge.or.vchange(fparm(ITFIX,ispc,isit),stfix(ispc))
      fparm(ITFIX,ispc,isit) = stfix(ispc)
      chnge=chnge.or.vchange(fparm(ISHFT,ispc,isit),sshft(ispc))
      fparm(ISHFT,ispc,isit) = sshft(ispc)
c
c
      chnge=chnge.or.(iparm(IIEXP,ispc,isit).ne.siexp(ispc))
      iparm(IIEXP,ispc,isit) = siexp(ispc)
      chnge=chnge.or.(iparm(IICOMB,ispc,isit).ne.sicomb(ispc))
      iparm(IICOMB,ispc,isit) = sicomb(ispc)
c any changes here can use Eval,Evec already calculated:
      if(chnge) matrxok(isit)=min(matrxok(isit),1)
c
      npt1 = snpt1(ispc)
      npt2 = snpt2(ispc)
c
      return
      end
