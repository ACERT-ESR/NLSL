c NLSL Version 1.5 beta 11/23/95
c----------------------------------------------------------------------
c                    =========================
c                       subroutine SCALEC
c                    =========================
c
c  Interprets a line containing the "scale" command. Used for adjusting
c  scale factor and automatic scaling for a specific site or range of
c  sites.
c
c  Syntax:
c
c       scale <index>|ALL|*  <value>|AUTO|FIX
c  
c          <index>     Site for which scale factor is to be adjusted
c
c          <value>     New value of scale factor
c
c          ALL|*       Apply command to all currently defined sites
c
c          AUTO        Enables automatic scaling of the site
c
c          FIX         Fixes scaling at present value 
c
c  Includes:
c     nlsdim.inc
c     expdat.inc
c     parcom.inc
c     stdio.inc 
c----------------------------------------------------------------------
      subroutine scalec( line )
      implicit none
      character*80 line
c
      include 'nlsdim.inc'
      include 'expdat.inc'
      include 'lmcom.inc'
      include 'parcom.inc'
      include 'mspctr.inc'
      include 'iterat.inc'
      include 'stdio.inc'
c
      integer i,iact,iflag,ival,ixsm,jx,jx1,jx2,lth
      double precision fval
      character*30 token
c
      integer NKEYWD
      parameter(NKEYWD=2)
      character*8 keywrd(NKEYWD)
c
      double precision ZERO
      parameter (ZERO=0.0d0)
c
      integer IAUTO,IFIX,IADJST
      parameter (IAUTO=1,IFIX=2,IADJST=3)
c
      integer itrim
      logical ftoken,itoken
      double precision enorm
      external enorm,ftoken,itoken,itrim
c
      data keywrd /'AUTO', 'FIX' /
c
c
c    ----------------------------------------
c     Look for an index identifying spectrum
c    ----------------------------------------
      call gettkn(line,token,lth)
c
c                          *** Site index expected
      if (lth.eq.0) then
         write (luout,1000)
         return
      end if
c
c     ----------------------------------------
c     Look for ALL keyword or wildcard index
c     ----------------------------------------
      if (.not.itoken(token,lth,ival)) then
         if (token(:lth).eq.'ALL'.or.token(:lth).eq.'*') then
            ival=-1
c                                             *** Illegal index
         else
            write(luout,1001) token(:lth)
            return
         end if
      end if
c
      if (ival.le.0) then
         jx1=1
         jx2=nsite
      else
         jx1=ival
         jx2=ival
      endif
c
c     --------------------
c      Look for a keyword
c     --------------------
 5    call gettkn(line,token,lth)
      lth=min(lth,8)
      if (lth.eq.0) go to 10
c
      call touppr(token,lth)
      do i=1,NKEYWD
         if (token(:lth).eq.keywrd(i)(:lth)) then
            iact=i
            go to 5
         end if
      end do
c
c     ---------------------------------------------
c      Not a keyword: is token a floating pt number?
c     ---------------------------------------------
      if (.not.ftoken(token,lth,fval)) then
         write (luout,1003) token(:itrim(token))
         if (luout.ne.luttyo) write(luttyo,1003) token(:itrim(token))
         go to 5
      else
         iact=IADJST
      end if
c
c     -----------------------------------------------
c      Do not fix scale factor if shifting is enabled
c     -----------------------------------------------
      if ((iact.eq.IFIX .or. iact.eq.IADJST) .and. ishglb.ne.0) then
         write (luout,1005)
         if (luout.ne.luttyo) write (luttyo,1005)
         return
      end if
c
c     --------------------------------------------
c      Perform action for specified range of sites
c     --------------------------------------------
 10   do jx=jx1,jx2
c
c                                        Float scale factor
         if (iact.eq.IAUTO) then
            iscal(jx)=1
c
c                                        Fix scale factor at present value
         else if (iact.eq.IFIX) then
            iscal(jx)=0
c
c                                        Set scale factor to given value
         else if (iact.eq.IADJST) then
            do i=1,nspc
               sfac(jx,i)=fval
            end do
            iscal(jx)=0
         end if
      end do
c
c    ------------------------------------------------------------------
c     Set flag indicating whether autoscaling is enabled for all sites
c    ------------------------------------------------------------------
      iscglb=1
      do i=1,nsite
         if (iscal(i).eq.0) iscglb=0
      end do
c
c     ------------------------------------------------------------
c      Repeat function calculation if a scale factor was "floated"
c      or changed by the user, and report the results
c     ------------------------------------------------------------
      if (iact.eq.IAUTO .or. iact.eq.IADJST) then
         call catchc( hltfit )
         call xpack( x, nprm )
c
         iflag=1
         call lfun(ndatot,nprm,x,fvec,fjac,MXPT,iflag)
c
         fnorm=enorm(ndatot,fvec)
         write(luout,1004) fnorm
         if (luout.ne.luttyo) write(luttyo,1004) fnorm
         call sclstt( luout )
c
         lmflag=1
         info=11
         call uncatchc( hltfit )
      end if
c
      return
c
c######################################################################
c
 1000 format('*** Site index expected ***')
 1001 format('*** Illegal index: ''',a,''' ***')
 1003 format('*** Unrecognized SCALE keyword: ''',a,''' ***')
 1004 format(/10x,'Recalculated RMS deviation =',g12.5/)
 1005 format('*** Scale cannot be fixed when shifting is enabled ***')
      end
