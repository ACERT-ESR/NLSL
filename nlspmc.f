c  VERSION 1.0  (NLSPMC version)   2/5/99
c----------------------------------------------------------------------
c                  ===================
c                     program NLSPMC
c                  ===================
c
c     MAIN PROGRAM FOR NON-LINEAR LEAST SQUARES FIT OF 2D-SPECTRA
c
c     This program reads in parameters and data for a nonlinear least-
c     squares fit using an EPRCGF-family slow-motional 2D spectral
c     calculation.
c
c     Modified to allow multiple components.  RC 7/98
c
c     ** Needs graphical interface **
c
c     Includes:
c        nlsdim.inc
c        nlsnam.inc
c        eprprm.inc
c        expdat.inc
c        parcom.inc
c        parmequ.inc
c        lmcomm.inc 
c        stdio.inc
c        iterat.inc
c 	 miscel.inc	new variables controlling sites/spectra choice
c
c######################################################################
c
      program nlspmc
cc	include 'CXML_include.F90'
c
      implicit none
c
c      include 'nlsdim.inc'
c      include 'nlsnam.inc'
c      include 'eprprm.inc'
c      include 'expdat.inc'
c      include 'parcom.inc'
c      include 'parmequ.inc'
c      include 'lmcomm.inc'
c      include 'stdio.inc'
c      include 'iterat.inc'
      include 'stdio.inc'
      include 'limits.inc'
      include 'names.inc'
      include 'parms.inc'
      include 'parmequ.inc'
      include 'simparm.inc'
      include 'miscel.inc'
c
      integer i,iflg,ioerr,j,lth
      logical fexist
      character line*(linelg), token*(wordlg), scratch*(wordlg),
     #	fileID*(wordlg), chr*2
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
      write (luttyo,1011)
c
c----------------------------------------------------------------------
c     Get next command from input stream (skip comments)
c----------------------------------------------------------------------
c
 25   if ((ihltcmd.ne.0).and.(nfiles.gt.0)) then
c
c if ^c is struck, close all macro files and start over:
c
         if (luout.ne.luttyo) write(luttyo,1004)
         write(luout,1004)
         do 27 i=1,nfiles
            close (ludisk+i)
 27      continue
         nfiles=0
         lucmd=luttyi
c         call uncatchc(ihltcmd)
c  reinitalize everything (later, maybe only do some things):
	 call nlsinit
      end if
c input a line:
      if (getlin(line)) then
c
c######################################################################
c get the first (next) token: 
         call gettkn(line,token,lth)
         call touppr(token,lth)
c if done (empty line or comment), get next line:
         if (lth.eq.0 .or. token.eq.'C' .or. token.eq.'/*') go to 25
c
c Now analyze the line received, first word is command, then extra 
c information may follow or not.
c
c----------------------------------------------------------------------
c        ASSIGN command
c  Assigns a given basis set to a given set of spectrum, site indices 
c command removed since basis(i,j) selects site and spectrum info.
c----------------------------------------------------------------------
c         if (token.eq.'ASSIGN') then
c            call assgnc(line)
c
c----------------------------------------------------------------------
c        AXIAL command
c----------------------------------------------------------------------
         if (token.eq.'AXIAL') then
            call convtc(line,3)
c
c----------------------------------------------------------------------
c        BASIS command
c----------------------------------------------------------------------
         else if (token.eq.'BASIS') then
c
            call basis(line,luout)
c
c----------------------------------------------------------------------
c        CARTESIAN command
c----------------------------------------------------------------------
         else if (token.eq.'CARTESIAN' .or. token.eq.'CART') then
            call convtc(line,1)
c
c----------------------------------------------------------------------
c        COMPONENTS command (alias COMP)
c----------------------------------------------------------------------
	 else if (token.eq.'COMPONENTS' .or. token.eq.'COMP') then
	    call comps(line)
c
c----------------------------------------------------------------------
c        DATA command
c----------------------------------------------------------------------
         else if (token.eq.'DATA') then
            call datac(line)
c
c----------------------------------------------------------------------
c        DEBUG command: open debug file identification
c----------------------------------------------------------------------
         else if (token.eq.'DEBUG') then
            idebug=1
            open(ludeb,file=dbname(:lthfnm),status='unknown',
     #           access='sequential',form='formatted',iostat=ioerr)
c
            ievec=0
            call gettkn(line,token,lth)
            call touppr(token,lth)
	    if (token.eq.'VECTOR') ievec=1
c
c----------------------------------------------------------------------
c        ECHO command
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
            end if
c
c----------------------------------------------------------------------
c        FIT command
c----------------------------------------------------------------------
         else if (token.eq.'FIT') then
            call fitc(line)
c
c----------------------------------------------------------------------
c        FIX command (alias REMOVE)
c----------------------------------------------------------------------
         else if (token.eq.'FIX' .or. token.eq.'REMOVE') then
            call fixc(line)
c
c----------------------------------------------------------------------
c        HELP command
c----------------------------------------------------------------------
         else if (token.eq.'HELP') then
            call helpc(line)
c
c----------------------------------------------------------------------
c        LET command
c----------------------------------------------------------------------
         else if (token.eq.'LET') then
            call letcmc(line)
c
c----------------------------------------------------------------------
c        LOG command: Set log file identification
c----------------------------------------------------------------------
         else if (token.eq.'LOG') then
            call gettkn(line,fileID,lth)
            if (lth.eq.0) then
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
c
               open(lulog,file=lgname(:lthfnm),status='unknown',
     #         access='sequential',form='formatted',iostat=ioerr)
c
                write(lulog,*)'log file of program nlspc.rc'
c
c  remove appending to existing log file (RC 7/25/96)
c
c               inquire(file=lgname(:lthfnm),exist=fexist)
c               if (fexist) then
c                                   * append the log file
c 30               read (lulog,'(a2)',end=32) chr
c                  go to 30
c               end if
 32            luout=lulog
            end if              
c
c----------------------------------------------------------------------
c        PARMS command
c----------------------------------------------------------------------
         else if (token.eq.'PARMS') then
            call parc(line)
c
c----------------------------------------------------------------------
c        QUIT command (alias EXIT)
c----------------------------------------------------------------------
         else if (token.eq.'QUIT'.or.token.eq.'EXIT') then
            goto 9999
c
c----------------------------------------------------------------------
c        READ command (alias CALL)
c        open a new input file and set I/O units appropriately
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
c       --- open input file if possible
c
               if (nfiles.ge.MXSPEC) then
                  write (luttyo,1050) fileID(:lth),MXSPEC
               else
                  nfiles=nfiles+1
                  lucmd=ludisk+nfiles
                  inquire(file=fileID(:lth),exist=fexist)
                  if (fexist) open(lucmd,file=fileID(:lth),
     #                    status='old',access='sequential',
     #                    form='formatted',iostat=ioerr)
                  if ((.not.fexist) .or. ioerr.ne.0) then
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
c
c                     call catchc( ihltcmd )
c                     
                  end if
               end if
c                              *** File identification not specified
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
c        SEARCH command
c----------------------------------------------------------------------
         else if (token.eq.'SEARCH') then
            call srchc(line)
c
c----------------------------------------------------------------------
c        SERIES command
c----------------------------------------------------------------------
         else if (token.eq.'SERIES') then
            call series(line)
c
c----------------------------------------------------------------------
c        SPECTRA command (alias SPEC)
c        specify multiple experimental data sets
c----------------------------------------------------------------------
c
	 else if (token.eq.'SPECTRA' .or. token.eq.'SPEC') then
	   call spectra(line)
c
c----------------------------------------------------------------------
c        SPHERICAL command
c----------------------------------------------------------------------
         else if (token.eq.'SPHERICAL' .or. token.eq.'SPHER') then
            call convtc(line,2)
c
c----------------------------------------------------------------------
c        STATUS command 
c----------------------------------------------------------------------
         else if (token.eq.'STATUS') then
            call statc(line,luttyo)
c
c
c----------------------------------------------------------------------
c        TEMPERATURE command (alias TEMP)
c----------------------------------------------------------------------
c	 else if (token.eq.'TEMPERATURE' .or. token.eq.'TEMP')
c     #		then
c	    tempvar = .true.
c
c----------------------------------------------------------------------
c        TILT command	Not needed yet. 
c----------------------------------------------------------------------
c	 else if (token.eq.'TILT') then
c	   psivar = .true.
c
c----------------------------------------------------------------------
c        UNIFORM command
c----------------------------------------------------------------------
         else if (token.eq.'UNIFORM') then
            uniflg=.true.
c
c----------------------------------------------------------------------
c        VARY command
c----------------------------------------------------------------------
         else if (token.eq.'VARY') then
            call varyc(line)
c
c----------------------------------------------------------------------
c        Unknown command  
c----------------------------------------------------------------------
         else
            write(luttyo,1040) token(:lth)
         end if
c
c----------------------------------------------------------------------
c    No more lines (getlin returned .false. or it received ^d).
c    Close current input unit; if there are no open macro files,
c    stop program.  Interactive usage will not stop as getlin reads 
c    from keyboard by default, in bach mode, the end of the last macro
c    file will terminate here.
c----------------------------------------------------------------------
c
      else
         if (nfiles.eq.0) then
            write(luttyo,1000)
            stop 'end of program NLSPMC'
         else
            close(lucmd)
            nfiles=nfiles-1
            if (nfiles.eq.0) then
               lucmd=luttyi
c               call uncatchc(ihltcmd)
            else
               lucmd=lucmd-1
            end if
         end if
      end if
      go to 25
c
c----------------------------------------------------------------------
c     Exit program
c----------------------------------------------------------------------
c
 9999 continue
      stop
c
c## format statements ###############################################
c
 1000 format(//,2x,70('#'),//)
 1004 format(/20x,'*** MACRO execution halted by user ***')
 1010 format(25x,'PROGRAM : NLSPMC',/,25x,'---------------',//)
c
 1011 format(/,'Program defaults to single spectra, single component.',
     #/,'You must first specify other choice with commands:',/,
     #'spectra # (number of data sets), components # (# in ',
     #'each spectra).',//,'Further, temperature variation', 
     #' is not allowed ',/,'unless first specified with temperature',
     #' command.',/ )
 1020 format('*** File identification must be specified ***'/)
 1021 format('*** Log file is not open ***')
 1030 format('*** Error opening or reading file ''',a,''' ***'/)
 1040 format('*** Unknown command : ''',a,''' ***')
 1050 format('*** Cannot open ''',a,''': more than',i2,
     #       ' read files ***')
 1060 format(a)
      end
