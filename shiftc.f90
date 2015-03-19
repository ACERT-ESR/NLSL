c NLSL Version 1.5.1 beta 1/20/96
c----------------------------------------------------------------------
c                    =========================
c                       subroutine SHIFTC
c                    =========================
c
c  Interprets a line containing the "shift" command, which adjusts
c  the shifting parameters for a selected spectrum, or for all spectra.
c
c  Syntax:
c
c       shift spectrum|index|ALL|*  <value>|ON|OFF|RESET
c  
c          <spectrum>|<index>   Spectrum name or index for which shifting 
c                               is to be changed
c
c          <value>              New value of shift parameter
c
c          ON                   Enable automatic shifting
c
c          OFF                  Disable automatic shifting (fix at current
c                               value)
c
c          RESET                Add the shift parameter to B0 and reset
c                               the shift parameter to zero
c
c  Includes:
c     nlsdim.inc
c     expdat.inc
c     parcom.inc
c     stdio.inc 
c----------------------------------------------------------------------
      subroutine shiftc( line )
c
      use nlsdim
      use expdat
      use lmcom
      use parcom
      use iterat
      use stdio
c
      implicit none
      character*80 line
c
      integer NKEYWD
      parameter(NKEYWD=3)
c
      double precision ZERO
      parameter (ZERO=0.0D0)
c
      integer ION,IOFF,IRESET,ISET,IADJST
      parameter (ION=1,IOFF=2,IRESET=3,ISET=4,IADJST=5)
c
      integer i,iact,iflag,ival,ixsm,jx,jx1,jx2,lth
      double precision fval
      character*8 keywrd(NKEYWD)
      character*30 token
c
      integer isfind,itrim,itemp(MXSPC)
      logical ftoken,itoken
      double precision enorm
      external enorm,ftoken,isfind,itoken,itrim
c
      data keywrd /'ON','OFF','RESET'/
c
c    ----------------------------------------
c     Look for an index identifying spectrum
c    ----------------------------------------
      call gettkn(line,token,lth)
c
c                             Spectrum name/index expected
      if (lth.eq.0) then
         write (luout,1000)
         return
      end if
c
      ixsm=isfind(token,lth)
      if (ixsm.eq.0) then
         if (.not.itoken(token,lth,ival)) then
            if (token(:lth).eq.'ALL'.or.token(:lth).eq.'*') then
               ival=-1
c                                             *** Illegal index
            else
               write(luout,1001) token(:lth)
               return
            end if
         end if
      else
         ival=abs( ixsm )
      end if
c
      if (ival.le.0) then
         jx1=1
         jx2=nspc
      else
         jx1=ival
         jx2=ival
      endif
c
c     ---------------------
c      Look for keywords
c     ---------------------
 5    call gettkn(line,token,lth)
      lth=min(lth,8)
      if (lth.eq.0) then
         iact=IADJST
         go to 10
      end if
c
      call touppr(token,lth)
      do i=1,NKEYWD
         if (token(:lth).eq.keywrd(i)(:lth)) then
            iact=i
            go to 5
         end if
      end do
c
c     -----------------------------------------------
c      Not a keyword: is token a floating pt number?
c     -----------------------------------------------
      if (.not.ftoken(token,lth,fval)) then
         write (luout,1003) token(:itrim(token))
         if (luout.ne.luttyo) write(luttyo,1003) token(:itrim(token))
         go to 5
      else
         iact=ISET
      end if
c
c    ----------------------------------------------
c     Make adjustment for specified range of sites
c    ----------------------------------------------
 10   do jx=jx1,jx2
         if (iact.eq.IRESET) then
            sb0(jx)=sb0(jx)+shft(jx)
            shft(jx)=ZERO
         else if (iact.eq.IADJST) then
            itemp(jx)=ishft(jx)
            ishft(jx)=1
         else if (iact.eq.ION) then
            ishft(jx)=1
         else if (iact.eq.IOFF) then
            ishft(jx)=0
         else if (iact.eq.ISET) then
            shft(jx)=fval
            ishft(jx)=0
         end if
      end do
c
c     ----------------------------------------------------------------
c      Set flag indicating whether shifting is enabled for any spectrum
c     ----------------------------------------------------------------
      ishglb=0
      do i=1,nspc
         if (ishft(i).ne.0) ishglb=1
      end do
c
c     --------------------------------------------------
c     If shifting was turned on or a shift was adjusted, 
c     recalculate the spectra
c     --------------------------------------------------
      if (iact.eq.ION .or.iact.eq.ISET .or.iact.eq.IADJST) then
         call catchc( hltfit )
         call xpack( x, nprm )
c
         iflag=1
         call lfun(ndatot,nprm,x,fvec,fjac,MXPT,iflag)
c
         fnorm=enorm(ndatot,fvec)
         write(luout,1046) fnorm
         if (luout.ne.luttyo) write(luttyo,1046) fnorm
         call sclstt( luout )
c
c    --------------------------------------------------
c     Reset shift parameters after single calculation
c    --------------------------------------------------
         do i=1,nspc
            shft(i)=shft(i)+tmpshft(i)
            tmpshft(i)=ZERO
         end do
         lmflag=1
         info=11
         call uncatchc( hltfit )
c
         if (iact.eq.IADJST) then
            do jx=jx1,jx2
               ishft(jx)=itemp(jx)
            end do
         end if
      end if
c
      return
c
c######################################################################
c
 1000 format('*** Spectrum ID expected ***')
 1001 format('*** Illegal index: ''',a,''' ***')
 1003 format('*** Illegal SHIFT value: ''',a,''' ***')
 1046 format(/10x,'Recalculated RMS deviation =',g12.5/)
      end
