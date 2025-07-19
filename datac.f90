c NLSL Version 1.5 beta 11/25/95
c----------------------------------------------------------------------
c                    =========================
c                       subroutine DATAC
c                    =========================
c
c Interprets a "data" command from the given line. The format of
c the command is as follows:
c
c data <dataid> { format <iform> bcmode <ibc> nspline <n> deriv <idrv> }
c   dataid : base name for datafile and associated output files:
c            <dataid>.DAT -- datafile
c            <dataid>.SPC -- fit spectrum
c   iform  : input format 
c               0 = ASCII 2-column format (default)
c               1 = PC-Genplot binary
c   ibc    : baseline correction mode 
c               0 = no baseline correction (default)
c               1 = subtract constant to give zero integral
c               n = fit a line to the n points at each end of
c                   of the spectrum
c   n      : number of splined points, or 0 for no spline
c            (default=zero)
c   idrv   : derivative mode
c               0 = 0th derivative (absorption) mode
c               1 = 1st derivative mode (default)
c
c----------------------------------------------------------------------
      subroutine datac( line )
c
      use nlsdim
      use expdat
      use lmcom
      use parcom
      use nlsnam
      use mspctr
      use stdio
      use rnddbl
c
      implicit none
      character line*80,token*30
c
      integer NKEYWD
      parameter (NKEYWD=9)
c
      integer i,iret,ispc,ival,ix,ixs,j,lth,ncmts
      character*80 comment(MXCMT)
c
      logical itoken
      integer getdat,itrim,newwindow
      external itoken,itrim,getdat,newwindow
c
      character*8 keywrd(NKEYWD)
      data keywrd /'ASCII','BINARY','BCMODE','DERIV','NSPLINE',
     #             'NOSHIFT','SHIFT','NORM','NONORM'/
c
c----------------------------------------------------------------------
c  Get the name of the datafile
c----------------------------------------------------------------------
      call gettkn(line,token,lth)
c
c----------------------------------------------------------------------
c  If (1) no data name was specified, or (2) the number of datafiles
c  already read in equals the number of spectra to be calculated (nser)
c  then reset the data buffer
c----------------------------------------------------------------------
      if (lth.eq.0.or.nspc.ge.nser) then
         nspc=0
         ndatot=0
         write(luttyo,1020)
         if (lth.eq.0) return
      end if
c
      nspc=nspc+1
      dataid(nspc)=token
c
c     -----------------------------------------------------
c      Force normalization for multisite-multispectrum fits
c     -----------------------------------------------------
      if (nsite.gt.1.and.nser.gt.1) normflg=1
c
c     --------------------
c      Look for a keyword 
c     --------------------
 5    call gettkn(line,token,lth)
      if (lth.ne.0) then
         lth=min(lth,8)
         call touppr(token,lth)
         do i=1,NKEYWD
            if (token(:lth).eq.keywrd(i)(:lth)) go to 7
         end do
c
         write (luttyo,1000) token(:lth)
         go to 5
c
c----------------------------------------------------------------------
c  Keyword found: assign appropriate value using next token
c----------------------------------------------------------------------
c
c        ---------------------------------------------------
c         The BCMODE, DERIV, and NSPLINE keywords require an 
c         integer argument
c        -------------------------------------------------------
 7       if (i.ge.3 .and. i.le.5) then
c
            call gettkn(line,token,lth)
c                                            *** No value given
            if (lth.eq.0) then
               write(luttyo,1003) keywrd(i)(:itrim(keywrd(i)))
               return
            end if
c
            if (itoken(token,lth,ival)) then 
c                                         *** BCMODE keyword
               if (i.eq.3) then
                  bcmode=ival
c                                         *** DERIV keyword
               else if (i.eq.4) then
                  drmode=ival
c                                         *** NSPLINE keyword
               else if (i.eq.5) then
                  nspline=ival
               end if
c                                         *** Illegal integer token!
            else
               write(luttyo,1010) token(:lth)
            end if
c
c----------------------------------------------------------------------
c --- Non-argument keywords: 
c     ASCII, BINARY, SHIFT, NOSHIFT, NORM, NONORM
c----------------------------------------------------------------------
c
         else
c                                         *** ASCII keyword
            if (i.eq.1) then
               inform=0
c                                         *** BINARY keyword
            else if (i.eq.2) then
               inform=1
c                                         *** NOSHIFT keyword
            else if (i.eq.6) then
               shftflg=0
c                                         *** SHIFT keyword
            else if (i.eq.7) then
               shftflg=1
c                                         *** NORM keyword
            else if (i.eq.8) then
               normflg=1
c                                         *** NONORM keyword
            else if (i.eq.9) then
               normflg=0
            end if
         end if
         go to 5
c
c----------------------------------------------------------------------
c     No more tokens on the line: read in datafile
c----------------------------------------------------------------------
      else
         iform(nspc)=inform
         ibase(nspc)=bcmode
         ishft(nspc)=shftflg
         npts(nspc)=nspline
         idrv(nspc)=drmode
         nrmlz(nspc)=normflg
c     
         call setdat( dataid(nspc) )
c
c        -----------------------------------------------------
c         Check whether an acceptable number of spline points
c         has been specified and modify if necessary 
c        -----------------------------------------------------

         if ( npts(nspc).ne.0 .and.
     #       (npts(nspc).lt.4 .or. npts(nspc).gt.MXSPT)
     #      ) then
            npts(nspc)=max(4,npts(nspc))
            npts(nspc)=min(MXSPT,npts(nspc))
            write(luttyo,1040) npts(nspc)
         end if
c
c        -------------------------------------------------------
c         Check whether there is enough storage for the new data
c        --------------------------------------------------------
         ix=ndatot+1
         if ( (npts(nspc).eq.0 .and. ix+MXSPT.gt.MXPT)
     #      .or.(npts(nspc).gt.0 .and. ix+npts(nspc).gt.MXPT) ) then
            write(luttyo,1050) MXPT
            nspc=nspc-1
            return
         end if
c
         iret=getdat(dtname,iform(nspc),ibase(nspc),npts(nspc),
     #          idrv(nspc),nrmlz(nspc),luout,comment,ncmts,data(ix),
     #          sbi(nspc),sdb(nspc),rmsn(nspc),spltmp,spltmp(1,2),
     #          spltmp(1,3) )
c
c                                        *** Error opening/reading datafile
         if (iret.ne.0) then
            write(luttyo,1060) dtname(:lthdnm)
            nspc=nspc-1
            return
         end if
c
c        -------------------------------------------------------------
c         Find smallest power of 2 greater than or equal to the number
c         of data points (for Fourier-Transform applications)
c        -------------------------------------------------------------
         nft(nspc)=1
 9       nft(nspc)=nft(nspc)*2
         if (nft(nspc).lt.npts(nspc)) go to 9
c
         ixsp(nspc)=ix
         shft(nspc)=0.0D0
         slb(nspc)=0.0D0
         srng(nspc)=sdb(nspc)*(npts(nspc)-1)
         if (ishft(nspc).ne.0) ishglb=1
c
         write (wndoid(nspc),1070) nspc,
     #        dataid(nspc)(:itrim(dataid(nspc))),char(0)
         if (rmsn(nspc).le.RNDOFF) rmsn(nspc)=1.0d0
         ndatot=ndatot+npts(nspc)
c
c        -----------------------------------
c         If nspc > nwin create a new window
c        -----------------------------------
c         if (nspc.gt.nwin) then
c            nwin = newwindow( wndoid )
c
c
c    ------------------------------------------------------------------
c    Call getwindows and plot files when the last datafile is read in
c      nser = number of spectra in the series
c      wndoid = character array of window I.D.'s
c               defined as character*20 wndoid(mxspc) in expdat.inc
c    ------------------------------------------------------------------
          if (nspc.eq.nser) then
              call getwndws( nspc, wndoid )

c           ------------------------------------------
c            Set fvec array to equal the data array 
c           (calculated spectrum=0) and plot the data  
c           -------------------------------------------
            do ispc=1,nspc
               do i=1,npts(ispc)
                  j=ixsp(ispc)+i-1
                  fvec(j)=data(j)
               end do
c
               ixs=ixsp(ispc)
               sfac(1,ispc)=1.0d0
c               call pltwin( data(ixs),fvec(ixs),spectr(ixs,1),
c     *                    sfac(1,ispc),MXPT,npts(ispc),nsite,ispc )
            call fstplt( data(ixsp(ispc)), fvec(ixsp(ispc)), 
     #           sbi(ispc), sdb(ispc), npts(ispc), ispc )
c
            end do
c
         end if
c
      end if
      return
c
c #### format statements ########################################
c
 1000 format('*** Unrecognized DATA keyword: ''',a,''' ***')
 1003 format('*** No value given for ''',a,''' ***')
 1010 format('*** Integer value expected: ''',a,''' ***')
 1020 format('*** Data buffer has been reset *** ')
 1040 format('*** Number of splined points reset to ',i4,' ***')
 1050 format('*** Maximum number of data points (',i4,') exceeded ***')
 1060 format(/13x,'*** Error opening or reading datafile ''',a,
     #       ''' ***'/)
 1070 format(i2,': ',a,a1)
      end

