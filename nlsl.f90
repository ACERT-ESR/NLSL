c Version 1.5.1b 11/6/95
c----------------------------------------------------------------------
c                  =========================
c                        program NLSL
c                  =========================
c
c     Main program for a nonlinear least-squares fit using an 
c     EPRLL-family slow-motional calculation. The options in running
c     this program are too numerous to detail here. Read the manual...
c     or better yet, wait until the movie comes  out.  DB 
c
c     Updated by DC and DEB to include graphical interface by
c     calls in this program, subroutine GETDAT, and subroutine LFUN.
c
c
c     Modules needed:      Additional modules needed by NLSINIT:
c        nlsdim.mod           basis.mod
c        nlsnam.mod           iterat.mod
c        eprprm.mod           mspct.mod
c        expdat.mod           tridag.mod
c        lmcom.mod
c        parcom.mod
c        symdef.mod
c        stdio.mod
c
c######################################################################
c
      program nlsl
c
      use nlsdim
      use nlsnam
      use eprprm
      use expdat
      use lmcom
      use parcom
      use symdef
      use stdio
c
      implicit none
c
      integer i,iflg,ioerr,j,lth
      logical fexist
      character line*80,token*30,scratch*30,fileID*30,chr1*1
c
      logical getlin
      external getlin
c
      write (luttyo,1000)
      write (luttyo,1010)
c
c    ----------------------
c    Initialize NLS system
c    ----------------------
c
      call nlsinit
c
c----------------------------------------------------------------------
c     Get next command from input stream (skip comments)
c----------------------------------------------------------------------
c
c Check for ^C during command input from indirect files.
c If detected, close all files and return to tty mode
c
 25   if (hltcmd.ne.0 .and. nfiles.gt.0) then
         do 26 i=1,nfiles
            close(ludisk+i)
 26      continue
         nfiles=0
         lucmd=luttyi
         call uncatchc( hltcmd )
      end if
c
      if (getlin(line)) then
c
c######################################################################
c
         call gettkn(line,token,lth)
         call touppr(token,lth)
         if (lth.eq.0 .or. token.eq.'C' .or. token.eq.'/*') go to 25
c
c----------------------------------------------------------------------
         if (token.eq.'ASSIGN') then
            call assgnc(line)
c
c----------------------------------------------------------------------
c         AXIAL command
c----------------------------------------------------------------------
         else if (token.eq.'AXIAL') then
            call convtc(line,AXIAL)
c
c----------------------------------------------------------------------
c        BASIS command
c----------------------------------------------------------------------
         else if (token.eq.'BASIS') then
            call basisc(line)
c
c----------------------------------------------------------------------
c         CARTESIAN command
c----------------------------------------------------------------------
         else if (token.eq.'CARTESIAN' .or. token.eq.'CART') then
            call convtc(line,CARTESIAN)
c
c----------------------------------------------------------------------
c         CONFIDENCE command
c----------------------------------------------------------------------
         else if (token.eq.'CONFIDENCE' .or. token .eq. 'CONF') then
            call confc(line)
c
c----------------------------------------------------------------------
c         CONFIDENCE command
c----------------------------------------------------------------------
         else if (token.eq.'CORRELATION' .or. token .eq. 'CORR') then
            call covrpt(luttyo)
c
c----------------------------------------------------------------------
c         DATA command
c----------------------------------------------------------------------
         else if (token.eq.'DATA') then
            call datac(line)
c
c----------------------------------------------------------------------
c         DELETE command
c----------------------------------------------------------------------
         else if (token.eq.'DELETE' .or. token.eq.'DEL' ) then
            call deletc(line)

c
c----------------------------------------------------------------------
c         ECHO command
c----------------------------------------------------------------------
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
c
c----------------------------------------------------------------------
c        FIT command
c----------------------------------------------------------------------
         else if (token.eq.'FIT') then
            call fitc(line)
c
c----------------------------------------------------------------------
c         FIX command (alias REMOVE)
c----------------------------------------------------------------------
         else if (token.eq.'FIX' .or. token.eq.'REMOVE') then
            call fixc(line)
c
c----------------------------------------------------------------------
c         HELP command
c----------------------------------------------------------------------
         else if (token.eq.'HELP') then
            call helpc(line)
c
c----------------------------------------------------------------------
c         LET command
c----------------------------------------------------------------------
         else if (token.eq.'LET') then
            call letc(line)
c
c----------------------------------------------------------------------
c           LOG command: Set log file identification
c----------------------------------------------------------------------
         else if (token.eq.'LOG') then
            call gettkn(line,fileID,lth)
            if (lth.eq.0) then
c
c                              *** File name not specified
               write (luttyo,1020)
c     
            else if (fileID.eq.'END' .or. fileID.eq.'end') then
               if (luout.eq.luttyo) then
                  write (luttyo,1021)
               else
                  close(lulog)
                  luout=luttyo
               end if
c     
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
c
c----------------------------------------------------------------------
c           PARMS command
c----------------------------------------------------------------------
         else if (token.eq.'PARMS') then
            call parc(line)
c
c----------------------------------------------------------------------
c           QUIT command (alias EXIT)
c----------------------------------------------------------------------
         else if (token.eq.'QUIT'.or.token.eq.'EXIT') then
            goto 9999
c
c----------------------------------------------------------------------
c           READ command (alias CALL)
c           open a new input file and set I/O units appropriately
c----------------------------------------------------------------------
c
         else if (token.eq.'READ' .or. token.eq.'CALL') then
c
c       --- get filename
c
            call gettkn(line,fileID,lth)
c
            if (lth.ne.0) then
c     
c        --- open input file if possible
c
               if (nfiles.ge.MXFILE) then
                  write (luttyo,1050) fileID(:lth),MXFILE
               else
                  nfiles=nfiles+1
                  lucmd=ludisk+nfiles
                  inquire(file=fileID(:lth),exist=fexist)
                  if (fexist) open(lucmd,file=fileID(:lth),
     #                 status='old',access='sequential',
     #                 form='formatted',iostat=ioerr)
c
                  if ((.not.fexist) .or. ioerr.ne.0) then
c
c                                               *** open error
                     write (luttyo,1030) fileID(:lth)
                     nfiles=nfiles-1
                     if (nfiles.eq.0) then
                        lucmd=luttyi
                     else
                        lucmd=lucmd-1
                     end if
c
                  else
                     files(nfiles)=fileID
                     call catchc( hltcmd )
                  end if
               end if
c
c                              *** File name not specified
            else
               write (luttyo,1020) 
            end if
c
c----------------------------------------------------------------------
c        RESET command
c----------------------------------------------------------------------
         else if (token.eq.'RESET') then
            call nlsinit
c
c----------------------------------------------------------------------
c        SCALE command
c----------------------------------------------------------------------
         else if (token.eq.'SCALE') then
            call scalec(line)
c
c----------------------------------------------------------------------
c        SEARCH command
c----------------------------------------------------------------------
         else if (token.eq.'SEARCH') then
            call srchc(line)

c----------------------------------------------------------------------
c        SERIES command
c----------------------------------------------------------------------
         else if (token.eq.'SERIES') then
            call series(line)
c     
c----------------------------------------------------------------------
c        SHIFT command
c----------------------------------------------------------------------
         else if (token.eq.'SHIFT') then
            call shiftc(line)
c
c----------------------------------------------------------------------
c        SITES command
c----------------------------------------------------------------------
         else if (token.eq.'SITES') then
            call sitec(line)
c
c----------------------------------------------------------------------
c         SPHERICAL command
c----------------------------------------------------------------------
         else if (token.eq.'SPHERICAL' .or. token.eq.'SPHER') then
            call convtc(line,SPHERICAL)
c
c----------------------------------------------------------------------
c        STATUS command 
c----------------------------------------------------------------------
         else if (token.eq.'STATUS') then
            call statc(line)
c
c----------------------------------------------------------------------
c        VARY command
c----------------------------------------------------------------------
         else if (token.eq.'VARY') then
            call varyc(line)
c
c----------------------------------------------------------------------
c        WRITE command
c----------------------------------------------------------------------
         else if (token.eq.'WRITE') then
            call writec( line )
c
c----------------------------------------------------------------------
c        Unknown command  
c----------------------------------------------------------------------
         else
            write(luttyo,1040) token(:lth)
         end if
c
c----------------------------------------------------------------------
c    No more lines (getlin returned .false.)
c    Close current input unit; if there are no open files, stop program
c----------------------------------------------------------------------
c
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
c
c----------------------------------------------------------------------
c       Exit program
c----------------------------------------------------------------------
c
 9999 continue
c
c
c     -----------------------------------
c      Close all windows before exiting
c     -----------------------------------
c      call shutwindows()
      call shtwndws()
      stop
c
c
c## format statements ###############################################
c
 1000 format(//,2x,70('#'),//)
 1010 format(25x,'PROGRAM : NLSL'/20x,'*** Version 1.5.1 beta ***'/
     #26x,'Mod 05/18/96'/
     #15x,'Recompiled by Zhichun Liang, 12/13/07'/
     #25x,'---------------',//)
 1020 format('*** File name must be specified ***'/)
 1021 format('*** Log file is not open ***')
 1022 format('*** Error',i3,' opening file ',a,' ***')
 1030 format('*** Error opening or reading file ''',a,''' ***'/)
 1040 format('*** Unknown command : ''',a,''' ***')
 1050 format('*** Cannot open ''',a,''': more than',i2,
     #       ' read files ***')
 1060 format(a)
      end



c----------------------------------------------------------------------
c                    =========================
c                       subroutine NLSINIT
c                    =========================
c
c     Initializes the following:
c       Data arrays
c       NLS parameter arrays
c       NLS convergence criteria
c----------------------------------------------------------------------
      subroutine nlsinit
c
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
      use stdio
c
      implicit none
c
      integer i,j
c
c----------------------------------------------------------------------
c    Initializations
c----------------------------------------------------------------------
      nfiles=0
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
c
c----------------------------------------
c     Initialize parameter arrays
c----------------------------------------
      do j=1,MXSITE
         do i=1,NVPRM
            fparm(i,j)=0.0d0
            ixx(i,j)=0
         end do
c
         do i=1,NIPRM
            iparm(i,j)=0
         end do
c
         do i=1,MXSPC
            sfac(j,i)=1.0d0
         end do
c
c--------------------------------------------------
c  Put in defaults for often-forgotten parameters
c--------------------------------------------------
         fparm(ICGTOL,j)=1.0D-3
         fparm(ISHIFT,j)=1.0D-3
c
c------------------------------------------------------------
c  Define an all-purpose default MTS 
c  (caveat: this is a conservative set, corresponding
c   to pretty slow motions at X-band!, and long calculation
c   times!)
c------------------------------------------------------------
         iparm(ILEMX,j)=10
         iparm(ILEMX+1,j)=9
         iparm(ILEMX+3,j)=6
         iparm(ILEMX+5,j)=2
         iparm(ILEMX+6,j)=2
      end do
c
c-------------------------------------------------------
c  Initialize tridiagonal matrix and basis index space
c-------------------------------------------------------
      do j=1,MXSPC
         do i=1,MXSITE
            modtd(i,j)=1
            ltd(i,j)=0
            basno(i,j)=0
         end do
      end do
      nexttd=1
      nextbs=1
      ntd=0
      nbas=0
c
c------------------------------------------------------------
c    Initialize data array parameters
c------------------------------------------------------------
      do i=1,MXSPC
         sbi(i)=0.0d0
         sdb(i)=0.0d0
         shft(i)=0.0d0
         tmpshft(i)=0.0d0
         sb0(i)=0.0d0
         spsi(i)=0.0d0
         sbi(i)=0.0d0
         iform(i)=0
         ibase(i)=0
         npts(i)=0
         ishft(i)=0
         ixsp(i)=1
         idrv(i)=1
      end do
      ishglb=0
c
c----------------------------------------
c  -- Enable autoscaling for all sites
c----------------------------------------
      do i=1,MXSITE
         iscal(i)=1
      end do
      iscglb=1
c
c------------------------------------------------------
c  -- Set initial values for NLS convergence criteria
c  -- First, call the routine to initalize F90 pointers
c------------------------------------------------------
      call lmcom_init
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
c
c--------------------------------------------------
c  -- Set initial values for line search parameters
c--------------------------------------------------
      pstep=0.05d0
      ptol=0.01d0
      pftol=0.1d0
      pbound=5.0d0
      mxpitr=100
c
c--------------------------------------------------
c  -- Set initial values for data input
c--------------------------------------------------
      inform=0
      bcmode=0
      shftflg=1
      drmode=1
      normflg=0
c
c      call shutwindows
        call shtwndws()
c      call initwindows
c
      return
      end
