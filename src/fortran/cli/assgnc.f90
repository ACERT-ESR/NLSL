! NLSL Version 1.3.2 2/27/94
!----------------------------------------------------------------------
!                    =========================
!                         subroutine ASSGNC
!                    =========================
!
!  Interprets command line for "assign" command
!
!  Assigns a given basis set to a given set of spectrum, site indices 
!
!  assign <basis_indx> { to } { spectrum <spec_indx> site <site_indx> } 
!
!  If no spectral or site indices are given, the assignment is made
!  for all currently defined sites/spectra  
!----------------------------------------------------------------------
      subroutine assgnc(line)
!
      use nlsdim
      use eprprm
      use parcom
      use tridag
      use basis
      use stdio
!
      implicit none
      character*80 line
!
      integer i,ibn,isi,isp,ival,ixsi1,ixsi2,ixsm,ixsp1,ixsp2,
     #        ixss(2),lth
      character*30 token
!
      integer isfind
      logical itoken
      external itoken,isfind
!
      ibn=0
!
!----------------------------------------------------------------------
!  Get the index/name name of the basis set file
!----------------------------------------------------------------------
      call gettkn(line,token,lth)
!
!                           *** No index/filename specified
      if (lth.eq.0) then
         write(luttyo,1000)
         return
      endif
!
      ixsm=isfind(token,lth)
      call touppr(token,lth)
      if (ixsm.eq.0) then
!                                             *** Illegal index
         if (.not.itoken(token,lth,ibn)) then
            write(luout,1001) token(:lth)
            if (luout.ne.luttyo) write(luttyo,1001) token(:lth)
            return
         end if
      else
         ibn=abs(ixsm)
      end if
!                                            *** Index out of range      
      if (ibn.lt.0 .or.ibn.gt.nbas) then
         write(luout,1004) ibn
         if (luout.ne.luttyo) write(luttyo,1004) ibn
         return
      end if
!
!----------------------------------------------------------------------
!     Got basis set index/name: now get site/spectrum indices
!----------------------------------------------------------------------
      call getss(line,ixss)
!
!     Set ranges of spectra, site indices
!
      if (ixss(2).le.0) then
         ixsi1=1
         ixsi2=MXSITE
      else
         ixsi1=ixss(2)
         ixsi2=ixss(2)
      end if
!
      if (ixss(1).le.0) then
         ixsp1=1
         ixsp2=MXSPC
      else
         ixsp1=ixss(1)
         ixsp2=ixss(1)
      end if
!
!----------------------------------------------------------------------
!    Now assign basis set 
!----------------------------------------------------------------------
      do isi=ixsi1,ixsi2 
         do isp=ixsp1,ixsp2
            basno(isi,isp)=ibn
            modtd(isi,isp)=1
         end do
      end do
      return
!
! #### Formats #######################################################
 1000 format('*** Basis set number/ID required ***')
 1001 format('*** Illegal index: ''',a,''' ***')
 1004 format('*** basis set ',i2,' is not defined ***')
      end 



!----------------------------------------------------------------------
!                    =========================
!                      subroutine GETSS
!                    =========================
!
!     Interprets command line as follows
!        {TO} { {SITE} <n> {SPECTRUM} <m> }
!
!     Returns a 2-vector with the site and spectrum specified on the
!     line (or zero for an index that wasn't specified)
!
!     If <n> and <m> are given withouth the SITE/SPECTRUM keywords,
!     the routine interprets them in the order site, spectrum
!----------------------------------------------------------------------
!
      subroutine getss(line,ixss)
!
      use stdio
!
      implicit none
      character*80 line
      integer ixss(2)
!
      integer i,ival,ixn,ixsm,lth
      character*30 token
!
      integer isfind
      logical itoken
      external isfind,itoken
!
      integer NKEYWD
      parameter(NKEYWD=3)
!
      character*8 keywrd(NKEYWD)
      data keywrd /'SPECTRUM','SITE','TO'/
!
      ixss(1)=0
      ixss(2)=0
      ixn=1
!
!----------------------------------------------------------------------
!     Look for a keyword or index
!----------------------------------------------------------------------
 5    call gettkn(line,token,lth)
!
      if (lth.ne.0) then
         lth=min(lth,8)
         call touppr(token,lth)
         do i=1,NKEYWD
            if (token(:lth).eq.keywrd(i)(:lth)) go to 7
         end do
!
!----------------------------------------------------------------------
!       Token is not a keyword: check whether it is a symbolic or
!       integer value
!----------------------------------------------------------------------
         ixsm=isfind(token,lth)
         if (ixsm.eq.0) then
            if (.not.itoken(token,lth,ival)) then
!                                             *** Illegal index
               write(luout,1001) token(:lth)
               go to 5
            end if
         else
            ival=abs(ixsm)
         end if
!
!    --- Assign index
!
         ixss(ixn)=ival
         ixn=1+mod(ixn,2)
!
!----------------------------------------------------------------------
!  Keyword found: set index accordingly (ignore "TO" keyword)
!----------------------------------------------------------------------
 7       if (i.lt.3) ixn=i
         go to 5
      end if
!
!----------------------------------------------------------------------
!     No more tokens on the line: return
!----------------------------------------------------------------------
!
      return
 1001 format('*** Illegal index: ''',a,''' ***')
      end
