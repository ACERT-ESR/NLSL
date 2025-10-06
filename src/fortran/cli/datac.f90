! NLSL Version 1.5 beta 11/25/95
!----------------------------------------------------------------------
!                    =========================
!                       subroutine DATAC
!                    =========================
!
! Interprets a "data" command from the given line. The format of
! the command is as follows:
!
! data <dataid> { format <iform> bcmode <ibc> nspline <n> deriv <idrv> }
!   dataid : base name for datafile and associated output files:
!            <dataid>.DAT -- datafile
!            <dataid>.SPC -- fit spectrum
!   iform  : input format 
!               0 = ASCII 2-column format (default)
!               1 = PC-Genplot binary
!   ibc    : baseline correction mode 
!               0 = no baseline correction (default)
!               1 = subtract constant to give zero integral
!               n = fit a line to the n points at each end of
!                   of the spectrum
!   n      : number of splined points, or 0 for no spline
!            (default=zero)
!   idrv   : derivative mode
!               0 = 0th derivative (absorption) mode
!               1 = 1st derivative mode (default)
!
!----------------------------------------------------------------------
      subroutine datac( line )
!
      use nlsdim
      use expdat
      use lmcom
      use parcom
      use nlsnam
      use mspctr
      use stdio
      use rnddbl
!
      implicit none
      character line*80,token*30
!
      integer NKEYWD
      parameter (NKEYWD=9)
!
      integer i,iret,ispc,ival,ix,ixs,j,lth,ncmts,normlz
      character*80 comment(MXCMT)
!
      logical itoken
      integer getdat,itrim,newwindow
      external itoken,itrim,getdat,newwindow
!
      character*8 keywrd(NKEYWD)
      data keywrd /'ASCII','BINARY','BCMODE','DERIV','NSPLINE',
     #             'NOSHIFT','SHIFT','NORM','NONORM'/
!
!----------------------------------------------------------------------
!  Get the name of the datafile
!----------------------------------------------------------------------
      call gettkn(line,token,lth)
!
!----------------------------------------------------------------------
!  If (1) no data name was specified, or (2) the number of datafiles
!  already read in equals the number of spectra to be calculated (nser)
!  then reset the data buffer
!----------------------------------------------------------------------
      if (lth.eq.0.or.nspc.ge.nser) then
         nspc=0
         ndatot=0
         write(luttyo,1020)
         if (lth.eq.0) return
      end if
!
      nspc=nspc+1
      dataid(nspc)=token
!
!     -----------------------------------------------------
!      Force normalization for multisite-multispectrum fits
!     -----------------------------------------------------
      if (nsite.gt.1.and.nser.gt.1) normflg=1
!
!     --------------------
!      Look for a keyword 
!     --------------------
 5    call gettkn(line,token,lth)
      if (lth.ne.0) then
         lth=min(lth,8)
         call touppr(token,lth)
         do i=1,NKEYWD
            if (token(:lth).eq.keywrd(i)(:lth)) go to 7
         end do
!
         write (luttyo,1000) token(:lth)
         go to 5
!
!----------------------------------------------------------------------
!  Keyword found: assign appropriate value using next token
!----------------------------------------------------------------------
!
!        ---------------------------------------------------
!         The BCMODE, DERIV, and NSPLINE keywords require an 
!         integer argument
!        -------------------------------------------------------
 7       if (i.ge.3 .and. i.le.5) then
!
            call gettkn(line,token,lth)
!                                            *** No value given
            if (lth.eq.0) then
               write(luttyo,1003) keywrd(i)(:itrim(keywrd(i)))
               return
            end if
!
            if (itoken(token,lth,ival)) then 
!                                         *** BCMODE keyword
               if (i.eq.3) then
                  bcmode=ival
!                                         *** DERIV keyword
               else if (i.eq.4) then
                  drmode=ival
!                                         *** NSPLINE keyword
               else if (i.eq.5) then
                  nspline=ival
               end if
!                                         *** Illegal integer token!
            else
               write(luttyo,1010) token(:lth)
            end if
!
!----------------------------------------------------------------------
! --- Non-argument keywords: 
!     ASCII, BINARY, SHIFT, NOSHIFT, NORM, NONORM
!----------------------------------------------------------------------
!
         else
!                                         *** ASCII keyword
            if (i.eq.1) then
               inform=0
!                                         *** BINARY keyword
            else if (i.eq.2) then
               inform=1
!                                         *** NOSHIFT keyword
            else if (i.eq.6) then
               shftflg=0
!                                         *** SHIFT keyword
            else if (i.eq.7) then
               shftflg=1
!                                         *** NORM keyword
            else if (i.eq.8) then
               normflg=1
!                                         *** NONORM keyword
            else if (i.eq.9) then
               normflg=0
            end if
         end if
         go to 5
!
!----------------------------------------------------------------------
!     No more tokens on the line: read in datafile
!----------------------------------------------------------------------
      else
         iform(nspc)=inform
         ibase(nspc)=bcmode
         ishft(nspc)=shftflg
         npts(nspc)=nspline
         idrv(nspc)=drmode
         nrmlz(nspc)=normflg
!     
         call setdat( dataid(nspc) )
!
!        -----------------------------------------------------
!         Check whether an acceptable number of spline points
!         has been specified and modify if necessary 
!        -----------------------------------------------------

         if ( npts(nspc).ne.0 .and.
     #       (npts(nspc).lt.4 .or. npts(nspc).gt.MXSPT)
     #      ) then
            npts(nspc)=max(4,npts(nspc))
            npts(nspc)=min(MXSPT,npts(nspc))
            write(luttyo,1040) npts(nspc)
         end if
!
!        -------------------------------------------------------
!         Check whether there is enough storage for the new data
!        --------------------------------------------------------
         ix=ndatot+1
         if ( (npts(nspc).eq.0 .and. ix+MXSPT.gt.MXPT)
     #      .or.(npts(nspc).gt.0 .and. ix+npts(nspc).gt.MXPT) ) then
            write(luttyo,1050) MXPT
            nspc=nspc-1
            return
         end if
!
         iret=getdat(dtname,iform(nspc),ibase(nspc),npts(nspc),
     #          idrv(nspc),nrmlz(nspc),luout,comment,ncmts,data(ix),
     #          sbi(nspc),sdb(nspc),rmsn(nspc),spltmp,spltmp(1,2),
     #          spltmp(1,3) )
!
!                                        *** Error opening/reading datafile
         if (iret.ne.0) then
            write(luttyo,1060) dtname(:lthdnm)
            nspc=nspc-1
            return
         end if
!
!        -------------------------------------------------------------
!         Find smallest power of 2 greater than or equal to the number
!         of data points (for Fourier-Transform applications)
!        -------------------------------------------------------------
         nft(nspc)=1
 9       nft(nspc)=nft(nspc)*2
         if (nft(nspc).lt.npts(nspc)) go to 9
!
         ixsp(nspc)=ix
         shft(nspc)=0.0D0
         slb(nspc)=0.0D0
         srng(nspc)=sdb(nspc)*(npts(nspc)-1)
         if (ishft(nspc).ne.0) ishglb=1
!
         write (wndoid(nspc),1070) nspc,
     #        dataid(nspc)(:itrim(dataid(nspc))),char(0)
         if (rmsn(nspc).le.RNDOFF) rmsn(nspc)=1.0d0
         ndatot=ndatot+npts(nspc)
!
!        -----------------------------------
!         If nspc > nwin create a new window
!        -----------------------------------
!         if (nspc.gt.nwin) then
!            nwin = newwindow( wndoid )
!
!
!    ------------------------------------------------------------------
!    Call getwindows and plot files when the last datafile is read in
!      nser = number of spectra in the series
!      wndoid = character array of window I.D.'s
!               defined as character*20 wndoid(mxspc) in expdat.inc
!    ------------------------------------------------------------------
          if (nspc.eq.nser) then
              call getwndws( nspc, wndoid )

!           ------------------------------------------
!            Set fvec array to equal the data array 
!           (calculated spectrum=0) and plot the data  
!           -------------------------------------------
            do ispc=1,nspc
               do i=1,npts(ispc)
                  j=ixsp(ispc)+i-1
                  fvec(j)=data(j)
               end do
!
               ixs=ixsp(ispc)
               sfac(1,ispc)=1.0d0
!               call pltwin( data(ixs),fvec(ixs),spectr(ixs,1),
!     *                    sfac(1,ispc),MXPT,npts(ispc),nsite,ispc )
            call fstplt( data(ixsp(ispc)), fvec(ixsp(ispc)), 
     #           sbi(ispc), sdb(ispc), npts(ispc), ispc )
!
            end do
!
         end if
!
      end if
      return
!
! #### format statements ########################################
!
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

