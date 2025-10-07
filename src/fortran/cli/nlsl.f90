! Version 1.5.1b 11/6/95
      subroutine procline(line)
!> @brief
!>  processes the command line for the main nlsl program
!> @details
!>  it's separated out to facilitate wrapping with python
!>
!>
!>
      use stdio
      use nlsnam
      use symdef
      implicit none
      character(80), intent(in) :: line
      character fileID*30,token*30,scratch*30
      logical fexist
      integer ioerr,lth
!######################################################################
!
         call gettkn(line,token,lth)
         call touppr(token,lth)
         if (lth.eq.0 .or. token.eq.'C' .or. token.eq.'/*') return
!
!----------------------------------------------------------------------
         if (token.eq.'ASSIGN') then
            call assgnc(line)
!
!----------------------------------------------------------------------
!         AXIAL command
!----------------------------------------------------------------------
         else if (token.eq.'AXIAL') then
            call convtc(line,AXIAL)
!
!----------------------------------------------------------------------
!        BASIS command
!----------------------------------------------------------------------
         else if (token.eq.'BASIS') then
            call basisc(line)
!
!----------------------------------------------------------------------
!         CARTESIAN command
!----------------------------------------------------------------------
         else if (token.eq.'CARTESIAN' .or. token.eq.'CART') then
            call convtc(line,CARTESIAN)
!
!----------------------------------------------------------------------
!         CONFIDENCE command
!----------------------------------------------------------------------
         else if (token.eq.'CONFIDENCE' .or. token .eq. 'CONF') then
            call confc(line)
!
!----------------------------------------------------------------------
!         CONFIDENCE command
!----------------------------------------------------------------------
         else if (token.eq.'CORRELATION' .or. token .eq. 'CORR') then
            call covrpt(luttyo)
!
!----------------------------------------------------------------------
!         DATA command
!----------------------------------------------------------------------
         else if (token.eq.'DATA') then
            call datac(line)
!
!----------------------------------------------------------------------
!         DELETE command
!----------------------------------------------------------------------
         else if (token.eq.'DELETE' .or. token.eq.'DEL' ) then
            call deletc(line)

!
!----------------------------------------------------------------------
!         ECHO command
!----------------------------------------------------------------------
         else if (token.eq.'ECHO') then
            call gettkn(line,token,lth)
            if (lth.gt.0) then
               scratch=token
            else
               scratch=' '
            end if
            call touppr(scratch,3)
            if (scratch.eq.'ON') then
               luecho=luttyo
            else if (scratch.eq.'OFF') then
               luecho=0
            else
               call ungett(token,lth,line)
               write(luttyo,1060) line
               if (luout.ne.luttyo) write (luout,1060) line
            end if
!
!----------------------------------------------------------------------
!        FIT command
!----------------------------------------------------------------------
         else if (token.eq.'FIT') then
            call fitc(line)
!
!----------------------------------------------------------------------
!         FIX command (alias REMOVE)
!----------------------------------------------------------------------
         else if (token.eq.'FIX' .or. token.eq.'REMOVE') then
            call fixc(line)
!
!----------------------------------------------------------------------
!         HELP command
!----------------------------------------------------------------------
         else if (token.eq.'HELP') then
            call helpc(line)
!
!----------------------------------------------------------------------
!         LET command
!----------------------------------------------------------------------
         else if (token.eq.'LET') then
            call letc(line)
!
!----------------------------------------------------------------------
!           LOG command: Set log file identification
!----------------------------------------------------------------------
         else if (token.eq.'LOG') then
            call gettkn(line,fileID,lth)
            if (lth.eq.0) then
!
!                              *** File name not specified
               write (luttyo,1020)
!     
            else if (fileID.eq.'END' .or. fileID.eq.'end') then
               if (luout.eq.luttyo) then
                  write (luttyo,1021)
               else
                  close(lulog)
                  luout=luttyo
               end if
!     
            else
               call setfil( fileID )
               open(lulog,file=lgname(:lthfnm),status='unknown',
     #              access='sequential',form='formatted',
     #              iostat=ioerr)
               if (ioerr.ne.0) then
                  write (luttyo,1022) ioerr,lgname(:lthfnm)
               else
                  luout=lulog
               end if
            end if
!
!----------------------------------------------------------------------
!           PARMS command
!----------------------------------------------------------------------
         else if (token.eq.'PARMS') then
            call parc(line)
!
!----------------------------------------------------------------------
!           QUIT command (alias EXIT)
!----------------------------------------------------------------------
         else if (token.eq.'QUIT'.or.token.eq.'EXIT') then
             exitprog = .true.
             return
!
!----------------------------------------------------------------------
!           READ command (alias CALL)
!           open a new input file and set I/O units appropriately
!----------------------------------------------------------------------
!
         else if (token.eq.'READ' .or. token.eq.'CALL') then
!
!       --- get filename
!
            call gettkn(line,fileID,lth)
!
            if (lth.ne.0) then
!     
!        --- open input file if possible
!
               if (nfiles.ge.MXFILE) then
                  write (luttyo,1050) fileID(:lth),MXFILE
               else
                  nfiles=nfiles+1
                  lucmd=ludisk+nfiles
                  inquire(file=fileID(:lth),exist=fexist)
                  if (fexist) open(lucmd,file=fileID(:lth),
     #                 status='old',access='sequential',
     #                 form='formatted',iostat=ioerr)
!
                  if ((.not.fexist) .or. ioerr.ne.0) then
!
!                                               *** open error
                     write (luttyo,1030) fileID(:lth)
                     nfiles=nfiles-1
                     if (nfiles.eq.0) then
                        lucmd=luttyi
                     else
                        lucmd=lucmd-1
                     end if
!
                  else
                     files(nfiles)=fileID
                     call catchc( hltcmd )
                  end if
               end if
!
!                              *** File name not specified
            else
               write (luttyo,1020) 
            end if
!
!----------------------------------------------------------------------
!        RESET command
!----------------------------------------------------------------------
         else if (token.eq.'RESET') then
            call nlsinit
!
!----------------------------------------------------------------------
!        SCALE command
!----------------------------------------------------------------------
         else if (token.eq.'SCALE') then
            call scalec(line)
!
!----------------------------------------------------------------------
!        SEARCH command
!----------------------------------------------------------------------
         else if (token.eq.'SEARCH') then
            call srchc(line)

!----------------------------------------------------------------------
!        SERIES command
!----------------------------------------------------------------------
         else if (token.eq.'SERIES') then
            call series(line)
!     
!----------------------------------------------------------------------
!        SHIFT command
!----------------------------------------------------------------------
         else if (token.eq.'SHIFT') then
            call shiftc(line)
!
!----------------------------------------------------------------------
!        SITES command
!----------------------------------------------------------------------
         else if (token.eq.'SITES') then
            call sitec(line)
!
!----------------------------------------------------------------------
!         SPHERICAL command
!----------------------------------------------------------------------
         else if (token.eq.'SPHERICAL' .or. token.eq.'SPHER') then
            call convtc(line,SPHERICAL)
!
!----------------------------------------------------------------------
!        STATUS command 
!----------------------------------------------------------------------
         else if (token.eq.'STATUS') then
            call statc(line)
!
!----------------------------------------------------------------------
!        VARY command
!----------------------------------------------------------------------
         else if (token.eq.'VARY') then
            call varyc(line)
!
!----------------------------------------------------------------------
!        WRITE command
!----------------------------------------------------------------------
         else if (token.eq.'WRITE') then
            call writec( line )
!
!----------------------------------------------------------------------
!        Unknown command  
!----------------------------------------------------------------------
         else
            write(luttyo,1040) token(:lth)
         end if
!
!## format statements ###############################################
!
 1020 format('*** File name must be specified ***'/)
 1021 format('*** Log file is not open ***')
 1022 format('*** Error',i3,' opening file ',a,' ***')
 1030 format('*** Error opening or reading file ''',a,''' ***'/)
 1040 format('*** Unknown command : ''',a,''' ***')
 1050 format('*** Cannot open ''',a,''': more than',i2,
     #       ' read files ***')
 1060 format(a)
      end
!----------------------------------------------------------------------
!                  =========================
!                        program NLSL
!                  =========================
!
!     Main program for a nonlinear least-squares fit using an 
!     EPRLL-family slow-motional calculation. The options in running
!     this program are too numerous to detail here. Read the manual...
!     or better yet, wait until the movie comes  out.  DB 
!
!     Updated by DC and DEB to include graphical interface by
!     calls in this program, subroutine GETDAT, and subroutine LFUN.
!
!
!     Modules needed:      Additional modules needed by NLSINIT:
!        nlsdim.mod           basis.mod
!        nlsnam.mod           iterat.mod
!        eprprm.mod           mspct.mod
!        expdat.mod           tridag.mod
!        lmcom.mod
!        parcom.mod
!        symdef.mod
!        stdio.mod
!
!######################################################################
!
      program nlsl
!
      use nlsdim
      use nlsnam
      use eprprm
      use expdat
      use lmcom
      use parcom
      use symdef
      use stdio
!
      implicit none
!
      integer i
      character line*80
!
      logical getlin
      external getlin
!
      write (luttyo,1000)
      write (luttyo,1010)
!
!    ----------------------
!    Initialize NLS system
!    ----------------------
!
      call nlsinit
!
!----------------------------------------------------------------------
!     Get next command from input stream (skip comments)
!----------------------------------------------------------------------
!
! Check for ^C during command input from indirect files.
! If detected, close all files and return to tty mode
!
! The following line is the start of the main loop
 25   if (hltcmd.ne.0 .and. nfiles.gt.0) then
         do 26 i=1,nfiles
            close(ludisk+i)
 26      continue
         nfiles=0
         lucmd=luttyi
         call uncatchc( hltcmd )
      end if
!
      if (getlin(line)) then
          call procline(line)
!----------------------------------------------------------------------
!    No more lines (getlin returned .false.)
!    Close current input unit; if there are no open files, stop program
!----------------------------------------------------------------------
!
      else
         if (nfiles.eq.0) then
            write(luttyo,1000)
            stop 'end of program NLSL'
         else
            close(lucmd)
            nfiles=nfiles-1
            if (nfiles.eq.0) then
               lucmd=luttyi
               call uncatchc( hltcmd )
            else
               lucmd=lucmd-1
            end if
         end if
      end if
      go to 25
!
!----------------------------------------------------------------------
!       Exit program
!----------------------------------------------------------------------
!
!
!
!     -----------------------------------
!      Close all windows before exiting
!     -----------------------------------
!      call shutwindows()
      call shtwndws()
      stop
!
!
!## format statements ###############################################
!
 1000 format(//,2x,70('#'),//)
 1010 format(25x,'PROGRAM : NLSL'/20x,'*** Version 1.5.1 beta ***'/
     #26x,'Mod 05/18/96'/
     #15x,'Recompiled by Zhichun Liang, 12/13/07'/
     #25x,'---------------',//)
      end



!----------------------------------------------------------------------
!                    =========================
!                       subroutine NLSINIT
!                    =========================
!
!     Initializes the following:
!       Data arrays
!       NLS parameter arrays
!       NLS convergence criteria
!----------------------------------------------------------------------
      subroutine nlsinit
!
      use nlsdim
      use nlsnam
      use eprprm
      use expdat
      use parcom
      use mspctr
      use tridag
      use basis
      use lmcom
      use iterat
      use parsav
      use stdio
!
!      implicit none
!
      integer i,j
!
!----------------------------------------------------------------------
!    Initializations
!----------------------------------------------------------------------
      exitprog=.false.
      warn=.false.
      hltcmd=0
      hltfit=0
      nfiles=0
      lthfnm=0
      lthdnm=0
      prname=' '
      lgname=' '
      trname=' '
      dtname=' '
      spname=' '
      files=' '
      lucmd=luttyi
      luout=luttyo
      luecho=luttyo
      nspc=0
      ndatot=0
      nprm=0
      iser=0
      nser=1
      nsite=1
      nwin=0
      nspline=0
!
!----------------------------------------
!     Initialize parameter arrays
!----------------------------------------
      do j=1,MXSITE
         do i=1,NFPRM
            fparm(i,j)=0.0d0
            ixx(i,j)=0
         end do
!
         do i=1,NIPRM
            iparm(i,j)=0
         end do
!
         do i=1,MXSPC
            sfac(j,i)=1.0d0
         end do
!
!----------------------------------------------------------
!     Initialize pointers to serve as meaningful aliases
!     for the various parameters stored in lengthy arrays
!----------------------------------------------------------
      call prm_ptr_init
!
!--------------------------------------------------
!  Put in defaults for often-forgotten parameters
!--------------------------------------------------
         fparm(ICGTOL,j)=1.0D-3
         fparm(ISHIFT,j)=1.0D-3
!
!------------------------------------------------------------
!  Define an all-purpose default MTS 
!  (caveat: this is a conservative set, corresponding
!   to pretty slow motions at X-band!, and long calculation
!   times!)
!------------------------------------------------------------
         iparm(ILEMX,j)=10
         iparm(ILEMX+1,j)=9
         iparm(ILEMX+3,j)=6
         iparm(ILEMX+5,j)=2
         iparm(ILEMX+6,j)=2
      end do

      do i=1,MXSPC
         serval(i)=0.0d0
      end do

!--------------------------------------------------
!  Reset saved-parameter storage
!--------------------------------------------------
      xsaved=.false.
      prsaved=.false.
      nstsav=0
      nspsav=0
      nprsav=0
      fxsave = 0.0d0
      fprsav = 0.0d0
      spsave = 0.0d0
      ixsave = 0
      iprsav = 0
      ixxsav = 0
      tagsav = ' '

!--------------------------------------------------
!  Reset variable-parameter bookkeeping arrays
!--------------------------------------------------
      do i=1,MXVAR
         prmax(i)=0.0d0
         prmin(i)=0.0d0
         prscl(i)=0.0d0
         xfdstp(i)=0.0d0
         ibnd(i)=0
         ixpr(i)=0
         ixst(i)=0
      end do
      do i=1,MXJCOL
         xerr(i)=0.0d0
         tag(i)=" "
      end do
      mtxclc=.false.
!
!--------------------------------------------------
!  Reset working parameter copies
!--------------------------------------------------
      do i=1,NFPRM
         fepr(i)=0.0d0
      end do
      do i=1,NIPRM
         iepr(i)=0
      end do
      a0=0.0d0
      g0=0.0d0
      w0=0.0d0
      expl=0.0d0
      expkxy=0.0d0
      expkzz=0.0d0
      faa=0.0d0
      fgm=0.0d0
      fwm=0.0d0
      fam=0.0d0
      fgd=0.0d0
      fad=0.0d0
      fwd=0.0d0
      cpot=0.0d0
      xlk=0.0d0
      itype=0
      ipt=0
      itm=0
      itd=0
      ipsi0=0
      lband=0
      kband=0
      ldelta=0
      kdelta=0
      lptmx=0
      kptmx=0
      neltot=0
      nelv=0
      nelre=0
      nelim=0
      ncgstp=0

!--------------------------------------------------
!  Reset iteration and fit bookkeeping
!--------------------------------------------------
      iter   = 0
      iwflag = 1
      fnorm  = 0.0d0
      chisqr = 0.0d0
      rdchsq = 0.0d0
      ch2bnd = 0.0d0
      chnbnd = 0.0d0
      delchi1 = 0.0d0
      qfit   = 0.0d0
      f2bnd  = 0.0d0
      fnbnd  = 0.0d0
      tbound = 0.0d0
      fnmin  = 0.0d0
      newitr = .false.
      covarOK = .false.
      xreset = .false.
!
!-------------------------------------------------------
!  Initialize tridiagonal matrix and basis index space
!-------------------------------------------------------
      do j=1,MXSPC
         do i=1,MXSITE
            modtd(i,j)=1
            ltd(i,j)=0
            ixtd(i,j)=0
            basno(i,j)=0
         end do
      end do
      do i=1,MXTDM
         tdspec(i)=0
         tdsite(i)=0
      end do
      do i=1,MXTDG
         alpha(i)=(0.0d0,0.0d0)
         beta(i)=(0.0d0,0.0d0)
      end do
      do i=1,MXDIM
         stv(i)=(0.0d0,0.0d0)
         y(i)=(0.0d0,0.0d0)
      end do
      nexttd=1
      nextbs=1
      ntd=0
      nbas=0
      ibasis = 0
      mts    = 0
      ixbas  = 0
      ltbas  = 0
      bsused = 0
      basisID = ' '
!
!------------------------------------------------------------
!    Initialize data array parameters
!------------------------------------------------------------
      do i=1,MXSPC
         sbi(i)=0.0d0
         sdb(i)=0.0d0
         shft(i)=0.0d0
         tmpshft(i)=0.0d0
         sb0(i)=0.0d0
         spsi(i)=0.0d0
         sphs(i)=0.0d0
         slb(i)=0.0d0
         srng(i)=0.0d0
         iform(i)=0
         ibase(i)=0
         nft(i)=0
         nrmlz(i)=0
         npts(i)=0
         ishft(i)=0
         ixsp(i)=1
         idrv(i)=1
         dataid(i)=' '
         wndoid(i)=' '
      end do
      ishglb=0
      shftflg=0
      normflg=0
      written=0
      inform=0
      bcmode=0
      drmode=0
      data   = 0.0d0
      spltmp = 0.0d0
      rmsn   = 0.0d0

!------------------------------------------------------------
!  Clear stored spectra
!------------------------------------------------------------
      spectr = 0.0d0
      wspec  = 0.0d0
!
!----------------------------------------
!  -- Enable autoscaling for all sites
!----------------------------------------
      do i=1,MXSITE
         iscal(i)=1
      end do
      iscglb=1
!
!------------------------------------------------------
!  -- Set initial values for NLS convergence criteria
!  -- First, call the routine to initalize F90 pointers
!------------------------------------------------------
      flmprm = 0.0d0
      ilmprm = 0
      call lmcom_init
      fjac = 0.0d0
      fvec = 0.0d0
      x    = 0.0d0
      diag = 0.0d0
      qtf  = 0.0d0
      corr = 0.0d0
      work1 = 0.0d0
      work2 = 0.0d0
      work3 = 0.0d0
      work4 = 0.0d0
      gnvec = 0.0d0
      gradf = 0.0d0
      tcov  = 0.0d0
      ipvt  = 0
      nprint=0
      lmflag=0
      xtol=1.0d-4
      ftol=1.0d-4
      gtol=1.0d-6
      factor=1.0d2
      maxitr=10
      nshift=8
      noneg=1
      srange=0.5d0
      maxev=100
      itrace=0
      lmflag=0
      iwflag=1
      confid=0.683
      ctol=1.0d-3
      njcol=0
      itridg=0
      iitrfl=0
      jacobi=0
      ixp1p=0
      ixs1p=0
      output=0
!
!--------------------------------------------------
!  -- Set initial values for line search parameters
!--------------------------------------------------
      pstep=0.05d0
      ptol=0.01d0
      pftol=0.1d0
      pbound=5.0d0
      mxpitr=100
!
!--------------------------------------------------
!  -- Set initial values for data input
!--------------------------------------------------
      inform=0
      bcmode=0
      shftflg=1
      drmode=1
      normflg=0
!
!      call shutwindows
        call shtwndws()
!      call initwindows
!
      return
      end
