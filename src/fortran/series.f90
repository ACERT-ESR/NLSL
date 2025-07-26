! NLSL Version 1.5 beta 11/24/95
!----------------------------------------------------------------------
!                    =========================
!                       subroutine SERIES
!                    =========================
!
! Interpret a "series" command from the given line.
! The command is of the form:
!
!    series <name> { {=} <value>, <value> { , <value>...} }
! 
!         name  : Name of a variable parameter
!         value : List of values assumed by the named parameter
!                 for each spectrum in the series
!
!   If a "series" command is given without an argument, it resets 
!   the calculation to a single spectum  
!----------------------------------------------------------------------
      subroutine series( line )
!
      use nlsdim
      use eprprm
      use expdat
      use mspctr
      use parcom
      use stdio
!
      implicit none
      character line*80
!
      integer i,ix,lth,nser1,iser1
      double precision fval
      character token*30
!
      integer ipfind
      logical ftoken
      external ftoken,ipfind
!
      iser1=iser
      nser1=nser
!
!----------------------------------------------------------------------
!     Get name of series parameter
!----------------------------------------------------------------------
      call gettkn(line,token,lth)
      lth=min(lth,6)
!                               *** No series variable specified: reset
      if (lth.eq.0) then
         iser=0
         nser=1
         write(luout,1002)
        return
      end if
!
!----------------------------------------
! Check whether parameter may be varied
!----------------------------------------
      call touppr(token,lth)
      ix=ipfind(token,lth)
!
      if (ix.eq.0.or.ix.gt.IB0) then
       write (luttyo,1001) token(:lth)
       return
      end if
!
      iser=ix
      nser=0
!
!----------------------------------------------------------------------
!         Get a list of values for the series parameter
!----------------------------------------------------------------------
   10 call gettkn(line,token,lth)
!
!------------------------------
! No more values--exit
!------------------------------
      if (lth.eq.0) then
!                               *** Must have been at least 2 in series
        if (nser.lt.2) then
           write(luttyo,1004)
           iser=iser1
           nser=nser1
        end if  
!
!    --- Copy spectral parameters from spectrum 1 for all new
!        series members (Note that B0, FIELDI, DFLD, RANGE, NFIELD, 
!        IDERIV supplied here as defaults will be obtained from the
!        datafiles if they are available 

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
!
              do i=1,MXSITE
                 sfac(i,ix)=1.0d0
              end do
!
           end do
!
        end if
        return
!
      end if
!
!     Check for optional '=' after parameter name
!
      if (token(:lth).eq.'='.and.nser.eq.0) goto 10
!
      if (ftoken(token,lth,fval)) then
        if (nser.ge.MXSPC) then
           write (luttyo,1005) MXSPC
           return
        end if
!
!                                *** Increment series list
        nser=nser+1
        serval(nser)=fval
!                                *** Illegal real number
      else
        write(luttyo,1003) token(:lth)
      end if
!
      go to 10   
!
! ###### format statements ######################################
!
 1001 format('*** ''',a,''' is not a variable parameter ***')
 1002 format('*** SERIES reset to single spectrum ***')
 1003 format('*** Real value expected: ''',a,''' ***')
 1004 format('*** SERIES must have at least 2 values ***')
 1005 format('*** SERIES may not have more than',i2,
     #       ' values: remainder ignored ***')
      end
