c Version 1.5.1 beta 2/3/96
c----------------------------------------------------------------------
c                    =========================
c                       subroutine STATC
c                    =========================
c
c  Prints out the current status of the fitting procedure to the given
c  logical unit number. Includes information on fitting variables and
c  user-specified convergence tolerances, etc.
c
c  Includes:
c     nlsdim.inc
c     eprprm.inc
c     expdat.inc
c     parcom.inc
c     lmcom.inc
c
c----------------------------------------------------------------------
      subroutine statc( line )
      implicit none
      character*80 line
c
      include 'nlsdim.inc'
      include 'eprprm.inc'
      include 'expdat.inc'
      include 'parcom.inc'
      include 'lmcom.inc'
      include 'stdio.inc'
c
      integer NKEYWD
      parameter(NKEYWD=7)
c
      integer i,lth
      logical fitflg,datflg,basflg,tdgflg,sclflg
      character*8 keywrd(NKEYWD)
      character*30 token
c
      integer itrim
      external itrim
c
      fitflg=.false.
      datflg=.false.
      basflg=.false.
      tdgflg=.false.
      sclflg=.false.
c
      data keywrd /'FIT','DATA','BASIS','TRIDIAG','SHIFT','SCALE','ALL'/
c
c----------------------------------------------------------------------
c         Look for a keyword
c----------------------------------------------------------------------
 5    call gettkn(line,token,lth)
      if (lth.ne.0) then
         lth=min(lth,8)
         call touppr(token,lth)
         do i=1,NKEYWD
            if (token(:lth).eq.keywrd(i)(:lth)) go to 7
         end do
         write (luout,1000) token(:itrim(token))
         if (luout.ne.luttyo) write (luttyo,1000) token(:itrim(token))
         go to 5
c
 7       if (i.eq.1) then
            fitflg=.true.
         else if (i.eq.2) then
            datflg=.true.
         else if (i.eq.3) then
            basflg=.true.
         else if (i.eq.4) then
            tdgflg=.true.
         else if (i.eq.5 .or. i.eq.6) then
            sclflg=.true.
         else if (i.eq.7) then
            fitflg=.true.
            datflg=.true.
            basflg=.true.
            tdgflg=.true.
            sclflg=.true.
         end if
         go to 5
      end if
c
      fitflg = .not.(datflg.or.basflg.or.tdgflg.or.sclflg)
      if (fitflg) then
         call fitstt( luout )
         call sclstt( luout )
         if (luout.ne.luttyo) then
            call fitstt( luttyo )
            call sclstt( luttyo )
         end if
      end if
c
      if (sclflg.and..not.fitflg) then
         call sclstt( luout )
         if (luout.ne.luttyo) call sclstt( luttyo )
      end if
c
      if (datflg) then
         call datstt( luout )
         if (luout.ne.luttyo) call datstt( luttyo )
      end if
c
      if (basflg) then
         call basstt( luout )
         if (luout.ne.luttyo) call basstt( luttyo )
      end if
c
      if (tdgflg) then
         call tdgstt( luout )
         if (luout.ne.luttyo) call tdgstt( luttyo )
      end if
c
      return
c
c######################################################################
c
 1000 format(' *** unrecognized STATUS keyword: ''',a,''' ***')
      end

c----------------------------------------------------------------------
c                    =========================
c                       subroutine FITSTT  
c                    =========================
c
c   Output status report for NLS fitting parameters on logical unit lu
c----------------------------------------------------------------------

      subroutine fitstt( lu )
      implicit none
      integer lu
c
      include 'nlsdim.inc'
      include 'eprprm.inc'
      include 'expdat.inc'
      include 'parcom.inc'
      include 'lmcom.inc'
      include 'errmsg.inc'
      include 'iterat.inc'
      include 'stdio.inc'
c
      integer i
c
      double precision wtdres
      external wtdres
c
c----------------------------------------------------------------------
c     parameters for fit command
c----------------------------------------------------------------------
      write(lu,1000)
      write(lu,1001) xtol,ftol,gtol
      write(lu,1002) maxitr,maxev,factor
c
c----------------------------------------------------------------------
c     Options for fit command
c----------------------------------------------------------------------
      if (nshift.lt.0) write (lu,1005)
      if (nshift.eq.0) write (lu,1006) 100.0d0*srange
      if (nshift.gt.0) write (lu,1007) nshift,nshift,100.0d0*srange
      if (noneg.eq.0) write (lu,1008)
      if (noneg.ne.0) write (lu,1009)
      if (output.eq.1) write (lu,1014)
      if (itrace.eq.1) write (lu,1015)
c
      if (iwflag.ne.0) then
         write (lu,1013) 'WEIGHTED'
      else
         write (lu,1013) 'UNWEIGHTED'
      end if
c
c----------------------------------------------------------------------
c     variable information
c----------------------------------------------------------------------
      if (nprm.ge.1) then
         write(lu,1010) nprm
         write(lu,1011)
c
         call xpack(x,nprm)
         do i=1,nprm
            write(lu,1012) tag(i),x(i),prscl(i),xfdstp(i)
         end do
      else
         write(lu,1020)
      end if
c
      if (info.lt.0) info=0
      if (info.gt.10) info=10
      if (lmflag.ne.0) then
         if (iwflag.ne.0) then
            chisqr=fnorm*fnorm
            rdchsq=chisqr/dfloat(ndatot-nprm)
         else
            chisqr=wtdres( fvec,ndatot,nspc,ixsp,npts,rmsn )
            rdchsq=chisqr/dfloat(ndatot-nprm)
         end if

         write(lu,1030) fnorm,chisqr,rdchsq,minerr(info)
      else
         write(lu,1031)
      end if
      return
c
c# formats ##############################################################
c
 1000 format(/,2x,20('#'),' FIT status ',20('#'))
 1001 format(5x,'xtol = ',1p,e8.1,3x,'ftol = ',1p,e8.1,3x,'gtol = ',
     #       1p,e8.1)
 1002 format(5x,'maxitr=',i3,'  maxfun=',i4,' bound=',f8.1//
     #2x,'Options in effect:')
 1005 format(5x,'NOSHIFT: No spectral shifting' )
 1006 format(5x,'SHIFT: Spectral shifting (srange= +/-',f5.1,'%)')
 1007 format(5x,'SHIFT',i4,': Spectral shifting first',i4,
     #' iterations (srange= +/-',f5.1,'%)')
 1008 format(5x,'NEG: negative scaling coefficients allowed')
 1009 format(5x,'NONEG: negative scaling coefficients not allowed')
 1010 format(/,2x,'VARY parameters : ',i1,' variables')
 1011 format(2x,56('-'),/,6x,'Parameter',8x,'Value',10x,'Scale',5x,
     #       'FD step',/,2x,56('-'))
 1012 format(7x,a9,5x,g14.7,3x,f6.1,5x,1p,e8.1)
 1013 format(5x,a,' residuals used for minimization')
 1014 format(5x,'WRITE: always write .spc file')
 1015 format(5x,'TRACE: write .trc file')
 1020 format(/,2x,'No parameters are being varied.')
 1030 format(/,2x,'Rsd. norm=',g14.7,2x,'Chi-sqr= ',g14.7,
     #         2x,'Red. chi-sqr=',g14.7,
     #       /,2x,'Status: ',a/)
 1031 format(/,2x,'No fit has been completed')
      end

c----------------------------------------------------------------------
c                    =========================
c                       subroutine DATSTT  
c                    =========================
c
c   Output status report for NLS datafile space on logical unit lu
c----------------------------------------------------------------------

      subroutine datstt( lu )
      implicit none
      integer lu
c
      include 'nlsdim.inc'
      include 'expdat.inc'
      include 'parcom.inc'
c
      character*2 shftopt,normopt
      integer i
c
      character*6 formopt(2)
      data formopt/'ASCII','BINARY'/
c
      write (lu,1000) nspc,ndatot,MXPT-ndatot
      if (nspc.lt.nser) write (lu,1004) nser
      if (nspc.gt.0) then
         write (lu,1001)
         do i=1,nspc
            shftopt='N'
            normopt='N'
            if (ishft(i).ne.0) shftopt='Y'
            if (nrmlz(i).ne.0) normopt='Y'
            write (lu,1002) i,dataID(i),npts(i),ibase(i),idrv(i),
     #                      shftopt,normopt,sb0(i),shft(i),rmsn(i) 
         end do
      end if
c
      shftopt=' '
      normopt=' '
      if (shftflg.eq.0) shftopt='NO'
      if (normflg.eq.0) normopt='NO'
      write (lu,1003) bcmode,nspline,shftopt,normopt,formopt(inform+1)
      return
c
c######################################################################
c
 1000 format(/2x,25('#'),' Datafile status ',25('#')//
     # i2,' data files and',i4,' data points in buffer'/
     # 'Space available for',i5,' points ')
 1001 format(/2x,'File   Filename   # Pts  Bsln  Drv',
     #'  Shft Nrm   B0     Offset  RMS noise'/70('-'))
 1002 format(2x,i3,4x,a,t20,2(2x,i4),3x,i1,2(4x,a1),f9.1,f9.2,g9.3)
 1003 format(/'Options in effect:'/5x,'bcmode=',i3,',  nspline=',i4,
     #', ',a2,'SHIFT, ',a2,'NORM, ',a/)
 1004 format(i2,' data files are required for the series fit')
      end






c----------------------------------------------------------------------
c                    =========================
c                       subroutine BASSTT  
c                    =========================
c
c   Output status report for NLS datafile space on logical unit lu
c----------------------------------------------------------------------

      subroutine basstt( lu )
      implicit none
      integer lu
c
      include 'nlsdim.inc'
      include 'expdat.inc'
      include 'parcom.inc'
      include 'basis.inc'
      include 'mtsdef.inc'
c
      integer i,j
c
      write (lu,1000) nbas,MXDIM-nextbs+1
      if (nbas.gt.0) then
         write (lu,1001)
         do i=1,nbas
            write (lu,1002) i,basisID(i),ltbas(i),(mts(j,i),j=1,NTRC)
     #                  
         end do
c
         write (lu,1003) (i,i=1,nsite)
         do j=1,nser
            write(lu,1004) dataID(j),(basno(i,j),i=1,nsite)
         end do
c
      end if
c
      write(lu,'(/)')
      return
c
c######################################################################
c
 1000 format(/2x,i2,' BASIS SETS in buffer; Space for ',i6,' elements'/)
 1001 format(2x,'Set   Identification       ndim    Trunc. indices'/
     #2x,60('-'))
 1002 format(2x,i3,3x,a,t29,i5,4x,6(i3,','),i2)
 1003 format(/2x,'Basis set usage:'/t24,5('Site',i2,2x:))
 1004 format(2x,a,t24,5(i3,5x)) 
      end



c----------------------------------------------------------------------
c                    =========================
c                        subroutine TDGSTT
c                    =========================
c
c
c   Output status report for NLS datafile space on logical unit lu
c----------------------------------------------------------------------

      subroutine tdgstt( lu )
      implicit none
      integer lu
c
      include 'nlsdim.inc'
      include 'expdat.inc'
      include 'parcom.inc'
      include 'tridag.inc'
c
      integer i,isi,isp
c
      write (lu,1000) ntd,MXTDG-nexttd+1
      if (ntd.gt.0) then
         write (lu,1001)
         do i=1,ntd
            isi=tdsite(i)
            isp=tdspec(i)
            if (modtd(isi,isp).ne.0) then
               write (lu,1002) i,isi,dataID(isp)
            else
               write (lu,1003) i,isi,dataID(isp),ltd(isi,isp)
            end if
         end do
      end if
c
      write(lu,'(/)')
      return
c
c######################################################################
c
 1000 format(/2x,i2,' TRIDIAGONAL MATRICES in buffer; Space for',i6,
     #       ' elements'/)
 1001 format(2x,'Matrix  Site  Spectrum           Lth'/2x,40('-'))
 1002 format(2x,i3,4x,i3,4x,a,t34,'(modified)')
 1003 format(2x,i3,4x,i3,4x,a,t34,i5)
      end



c----------------------------------------------------------------------
c                    =========================
c                       subroutine SCLSTT
c                    =========================
c
c   Output shifting and scaling information to specified logical unit
c----------------------------------------------------------------------
      subroutine sclstt(lu)
      implicit none
      integer lu
c
      include 'nlsdim.inc'
      include 'parcom.inc'
      include 'expdat.inc'
      include 'mspctr.inc'
      include 'nlsnam.inc'
c
      integer i,j,k,nsp
      double precision stot(MXSPC)
      character*1 shflg(MXSPC),scflg(MXSITE)
c
c     ----------------------------------------------------------------
c     Calculate total of scaling factor in order to report percentages 
c     Set flags showing which spectra are enabled for shifting
c     ----------------------------------------------------------------
      do j=1,nspc
         shflg(j)='N'
         if (ishft(j).ne.0) shflg(j)='Y'
c
         stot(j)=0.0d0
         do i=1,nsite
            stot(j)=stot(j)+sfac(i,j)
         end do
         stot(j)=stot(j)/1.0d2
      end do
c
c     ---------------------------------------------
c     Set flags showing which sites are autoscaled
c     ---------------------------------------------
      do i=1,nsite
         scflg(i)=' '
         if (iscal(i).ne.0) scflg(i)='*'
      end do
c
      if (nsite.gt.1.and.nspc.gt.1) then
         nsp=1
      else
         nsp=nspc
      end if
c
c     -----------------------------------------------------
c     Print out scaling information for each spectrum/site
c     -----------------------------------------------------
      if (nspc.gt.0) then
         if (nsite.gt.1) then
              write(lu,1047) (i,i=1,nsite)
            do j=1,nsp
              write(lu,1048) (sfac(i,j),sfac(i,j)/stot(j), i=1,nsite)
            end do
         end if
c
c     --------------------------------------------------
c      Print out shifting information for each spectrum
c     --------------------------------------------------
         k=lthdnm-4
         write (lu,1049) (dataid(j)(:k),shft(j),shflg(j),sb0(j),
     #                    sb0(j)+shft(j),j=1,nspc)
      endif
      return
c
c ##### formats ######################################################
c
 1047 format(/2x,'Scales:',t10,5(11x,'site',i2,2x))
 1048 format(12x,t16,5(g12.4,'(',f5.1,'%)  '))
 1049 format(/2x,'Shifts:  file',t33,'shift',2x,'auto',2x,4x,'B0',
     #       6x,'B0(eff)'/(12x,a,t30,f8.2,3x,a1,3x,f9.2,1x,f9.2))
      end

