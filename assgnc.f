c NLSL Version 1.3.2 2/27/94
c----------------------------------------------------------------------
c                    =========================
c                         subroutine ASSGNC
c                    =========================
c
c  Interprets command line for "assign" command
c
c  Assigns a given basis set to a given set of spectrum, site indices 
c
c  assign <basis_indx> { to } { spectrum <spec_indx> site <site_indx> } 
c
c  If no spectral or site indices are given, the assignment is made
c  for all currently defined sites/spectra  
c----------------------------------------------------------------------
      subroutine assgnc(line)
      implicit none
      character*80 line
c
      include 'nlsdim.inc'
      include 'eprprm.inc'
      include 'parcom.inc'
      include 'tridag.inc'
      include 'basis.inc'
      include 'stdio.inc'
c
c
      integer i,ibn,isi,isp,ival,ixsi1,ixsi2,ixsm,ixsp1,ixsp2,
     #        ixss(2),lth
      character*30 token
c
      integer isfind
      logical itoken
      external itoken,isfind
c
      ibn=0
c
c----------------------------------------------------------------------
c  Get the index/name name of the basis set file
c----------------------------------------------------------------------
      call gettkn(line,token,lth)
c
c                           *** No index/filename specified
      if (lth.eq.0) then
         write(luttyo,1000)
         return
      endif
c
      ixsm=isfind(token,lth)
      call touppr(token,lth)
      if (ixsm.eq.0) then
c                                             *** Illegal index
         if (.not.itoken(token,lth,ibn)) then
            write(luout,1001) token(:lth)
            if (luout.ne.luttyo) write(luttyo,1001) token(:lth)
            return
         end if
      else
         ibn=abs(ixsm)
      end if
c                                            *** Index out of range      
      if (ibn.lt.0 .or.ibn.gt.nbas) then
         write(luout,1004) ibn
         if (luout.ne.luttyo) write(luttyo,1004) ibn
         return
      end if
c
c----------------------------------------------------------------------
c     Got basis set index/name: now get site/spectrum indices
c----------------------------------------------------------------------
      call getss(line,ixss)
c
c     Set ranges of spectra, site indices
c
      if (ixss(2).le.0) then
         ixsi1=1
         ixsi2=MXSITE
      else
         ixsi1=ixss(2)
         ixsi2=ixss(2)
      end if
c
      if (ixss(1).le.0) then
         ixsp1=1
         ixsp2=MXSPC
      else
         ixsp1=ixss(1)
         ixsp2=ixss(1)
      end if
c
c----------------------------------------------------------------------
c    Now assign basis set 
c----------------------------------------------------------------------
      do isi=ixsi1,ixsi2 
         do isp=ixsp1,ixsp2
            basno(isi,isp)=ibn
            modtd(isi,isp)=1
         end do
      end do
      return
c
c #### Formats #######################################################
 1000 format('*** Basis set number/ID required ***')
 1001 format('*** Illegal index: ''',a,''' ***')
 1004 format('*** basis set ',i2,' is not defined ***')
      end 



c----------------------------------------------------------------------
c                    =========================
c                      subroutine GETSS
c                    =========================
c
c     Interprets command line as follows
c        {TO} { {SITE} <n> {SPECTRUM} <m> }
c
c     Returns a 2-vector with the site and spectrum specified on the
c     line (or zero for an index that wasn't specified)
c
c     If <n> and <m> are given withouth the SITE/SPECTRUM keywords,
c     the routine interprets them in the order site, spectrum
c----------------------------------------------------------------------
c
      subroutine getss(line,ixss)
      implicit none
      character*80 line
      integer ixss(2)
c
      integer i,ival,ixn,ixsm,lth
      character*30 token
c
      integer isfind
      logical itoken
      external isfind,itoken
c
      integer NKEYWD
      parameter(NKEYWD=3)
c
      include 'stdio.inc'
c
      character*8 keywrd(NKEYWD)
      data keywrd /'SPECTRUM','SITE','TO'/
c
      ixss(1)=0
      ixss(2)=0
      ixn=1
c
c----------------------------------------------------------------------
c     Look for a keyword or index
c----------------------------------------------------------------------
 5    call gettkn(line,token,lth)
c
      if (lth.ne.0) then
         lth=min(lth,8)
         call touppr(token,lth)
         do i=1,NKEYWD
            if (token(:lth).eq.keywrd(i)(:lth)) go to 7
         end do
c
c----------------------------------------------------------------------
c       Token is not a keyword: check whether it is a symbolic or
c       integer value
c----------------------------------------------------------------------
         ixsm=isfind(token,lth)
         if (ixsm.eq.0) then
            if (.not.itoken(token,lth,ival)) then
c                                             *** Illegal index
               write(luout,1001) token(:lth)
               go to 5
            end if
         else
            ival=abs(ixsm)
         end if
c
c    --- Assign index
c
         ixss(ixn)=ival
         ixn=1+mod(ixn,2)
c
c----------------------------------------------------------------------
c  Keyword found: set index accordingly (ignore "TO" keyword)
c----------------------------------------------------------------------
 7       if (i.lt.3) ixn=i
         go to 5
      end if
c
c----------------------------------------------------------------------
c     No more tokens on the line: return
c----------------------------------------------------------------------
c
      return
 1001 format('*** Illegal index: ''',a,''' ***')
      end
