! NLSL Version 1.5 beta 11/23/95
!----------------------------------------------------------------------
!                    =========================
!                       subroutine SCALEC
!                    =========================
!
!  Interprets a line containing the "scale" command. Used for adjusting
!  scale factor and automatic scaling for a specific site or range of
!  sites.
!
!  Syntax:
!
!       scale <index>|ALL|*  <value>|AUTO|FIX
!  
!          <index>     Site for which scale factor is to be adjusted
!
!          <value>     New value of scale factor
!
!          ALL|*       Apply command to all currently defined sites
!
!          AUTO        Enables automatic scaling of the site
!
!          FIX         Fixes scaling at present value 
!
!  Includes:
!     nlsdim.inc
!     expdat.inc
!     parcom.inc
!     stdio.inc 
!----------------------------------------------------------------------
      subroutine scalec( line )
!
      use nlsdim
      use expdat
      use lmcom
      use parcom
      use mspctr
      use iterat
      use stdio
!
      implicit none
      character*80 line
!
      integer i,iact,iflag,ival,ixsm,jx,jx1,jx2,lth
      double precision fval
      character*30 token
!
      integer NKEYWD
      parameter(NKEYWD=2)
      character*8 keywrd(NKEYWD)
!
      double precision ZERO
      parameter (ZERO=0.0d0)
!
      integer IAUTO,IFIX,IADJST
      parameter (IAUTO=1,IFIX=2,IADJST=3)
!
      integer itrim
      logical ftoken,itoken
      double precision enorm
      external enorm,ftoken,itoken,itrim
!
      data keywrd /'AUTO', 'FIX' /
!
!
!    ----------------------------------------
!     Look for an index identifying spectrum
!    ----------------------------------------
      call gettkn(line,token,lth)
!
!                          *** Site index expected
      if (lth.eq.0) then
         write (luout,1000)
         return
      end if
!
!     ----------------------------------------
!     Look for ALL keyword or wildcard index
!     ----------------------------------------
      if (.not.itoken(token,lth,ival)) then
         if (token(:lth).eq.'ALL'.or.token(:lth).eq.'*') then
            ival=-1
!                                             *** Illegal index
         else
            write(luout,1001) token(:lth)
            return
         end if
      end if
!
      if (ival.le.0) then
         jx1=1
         jx2=nsite
      else
         jx1=ival
         jx2=ival
      endif
!
!     --------------------
!      Look for a keyword
!     --------------------
 5    call gettkn(line,token,lth)
      lth=min(lth,8)
      if (lth.eq.0) go to 10
!
      call touppr(token,lth)
      do i=1,NKEYWD
         if (token(:lth).eq.keywrd(i)(:lth)) then
            iact=i
            go to 5
         end if
      end do
!
!     ---------------------------------------------
!      Not a keyword: is token a floating pt number?
!     ---------------------------------------------
      if (.not.ftoken(token,lth,fval)) then
         write (luout,1003) token(:itrim(token))
         if (luout.ne.luttyo) write(luttyo,1003) token(:itrim(token))
         go to 5
      else
         iact=IADJST
      end if
!
!     -----------------------------------------------
!      Do not fix scale factor if shifting is enabled
!     -----------------------------------------------
      if ((iact.eq.IFIX .or. iact.eq.IADJST) .and. ishglb.ne.0) then
         write (luout,1005)
         if (luout.ne.luttyo) write (luttyo,1005)
         return
      end if
!
!     --------------------------------------------
!      Perform action for specified range of sites
!     --------------------------------------------
 10   do jx=jx1,jx2
!
!                                        Float scale factor
         if (iact.eq.IAUTO) then
            iscal(jx)=1
!
!                                        Fix scale factor at present value
         else if (iact.eq.IFIX) then
            iscal(jx)=0
!
!                                        Set scale factor to given value
         else if (iact.eq.IADJST) then
            do i=1,nspc
               sfac(jx,i)=fval
            end do
            iscal(jx)=0
         end if
      end do
!
!    ------------------------------------------------------------------
!     Set flag indicating whether autoscaling is enabled for all sites
!    ------------------------------------------------------------------
      iscglb=1
      do i=1,nsite
         if (iscal(i).eq.0) iscglb=0
      end do
!
!     ------------------------------------------------------------
!      Repeat function calculation if a scale factor was "floated"
!      or changed by the user, and report the results
!     ------------------------------------------------------------
      if (iact.eq.IAUTO .or. iact.eq.IADJST) then
         call catchc( hltfit )
         call xpack( x, nprm )
!
         iflag=1
         call lfun(ndatot,nprm,x,fvec,fjac,MXPT,iflag)
!
         fnorm=enorm(ndatot,fvec)
         write(luout,1004) fnorm
         if (luout.ne.luttyo) write(luttyo,1004) fnorm
         call sclstt( luout )
!
         lmflag=1
         info=11
         call uncatchc( hltfit )
      end if
!
      return
!
!######################################################################
!
 1000 format('*** Site index expected ***')
 1001 format('*** Illegal index: ''',a,''' ***')
 1003 format('*** Unrecognized SCALE keyword: ''',a,''' ***')
 1004 format(/10x,'Recalculated RMS deviation =',g12.5/)
 1005 format('*** Scale cannot be fixed when shifting is enabled ***')
      end
