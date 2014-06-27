c NLSL Version 1.5 beta 11/24/95
c----------------------------------------------------------------------
c                    =========================
c                       subroutine SERIES
c                    =========================
c
c Interpret a "series" command from the given line.
c The command is of the form:
c
c    series <name> { {=} <value>, <value> { , <value>...} }
c 
c         name  : Name of a variable parameter
c         value : List of values assumed by the named parameter
c                 for each spectrum in the series
c
c   If a "series" command is given without an argument, it resets 
c   the calculation to a single spectum  
c----------------------------------------------------------------------
      subroutine series( line )
      implicit none
      character line*80
c
      include 'nlsdim.inc'
      include 'eprprm.inc'
      include 'expdat.inc'
      include 'mspctr.inc'
      include 'parcom.inc'
      include 'stdio.inc'
c
      integer i,ix,lth,nser1,iser1
      double precision fval
      character token*30
c
      integer ipfind
      logical ftoken
      external ftoken,ipfind
c
      iser1=iser
      nser1=nser
c
c----------------------------------------------------------------------
c     Get name of series parameter
c----------------------------------------------------------------------
      call gettkn(line,token,lth)
      lth=min(lth,6)
c                               *** No series variable specified: reset
      if (lth.eq.0) then
         iser=0
         nser=1
         write(luout,1002)
        return
      end if
c
c----------------------------------------
c Check whether parameter may be varied
c----------------------------------------
      call touppr(token,lth)
      ix=ipfind(token,lth)
c
      if (ix.eq.0.or.ix.gt.IB0) then
       write (luttyo,1001) token(:lth)
       return
      end if
c
      iser=ix
      nser=0
c
c----------------------------------------------------------------------
c         Get a list of values for the series parameter
c----------------------------------------------------------------------
   10 call gettkn(line,token,lth)
c
c------------------------------
c No more values--exit
c------------------------------
      if (lth.eq.0) then
c                               *** Must have been at least 2 in series
        if (nser.lt.2) then
           write(luttyo,1004)
           iser=iser1
           nser=nser1
        end if  
c
c    --- Copy spectral parameters from spectrum 1 for all new
c        series members (Note that B0, FIELDI, DFLD, RANGE, NFIELD, 
c        IDERIV supplied here as defaults will be obtained from the
c        datafiles if they are available 

        if (nser.gt.nser1) then
           do ix=nser1+1,nser
              slb(ix)=slb(1)
              sphs(ix)=sphs(1)
              spsi(ix)=spsi(1)
              sb0(ix)=sb0(1)
              sbi(ix)=sbi(1)
              sdb(ix)=sdb(1)
              srng(ix)=srng(1)
              npts(ix)=npts(1)
              nft(ix)=nft(1)
              idrv(ix)=idrv(1)
c
              do i=1,MXSITE
                 sfac(i,ix)=1.0d0
              end do
c
           end do
c
        end if
        return
c
      end if
c
c     Check for optional '=' after parameter name
c
      if (token(:lth).eq.'='.and.nser.eq.0) goto 10
c
      if (ftoken(token,lth,fval)) then
        if (nser.ge.MXSPC) then
           write (luttyo,1005) MXSPC
           return
        end if
c
c                                *** Increment series list
        nser=nser+1
        serval(nser)=fval
c                                *** Illegal real number
      else
        write(luttyo,1003) token(:lth)
      end if
c
      go to 10   
c
c ###### format statements ######################################
c
 1001 format('*** ''',a,''' is not a variable parameter ***')
 1002 format('*** SERIES reset to single spectrum ***')
 1003 format('*** Real value expected: ''',a,''' ***')
 1004 format('*** SERIES must have at least 2 values ***')
 1005 format('*** SERIES may not have more than',i2,
     #       ' values: remainder ignored ***')
      end
