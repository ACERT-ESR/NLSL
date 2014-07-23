c NLSPMC VERSION (VERSION 1.0)  2/5/99
c----------------------------------------------------------------------
c                  =========================
c                      subroutine helpc
c                  =========================
c
c     Main program for a nonlinear least-squares fit using an 
c     EPRP-family slow-motional calculation. The options in running
c     this program are too numerous to detail here. Read the manual...
c     or better yet, wait until the movie comes out.  DB 
c
c----------------------------------------------------------------------
      subroutine helpc(line)
      implicit none
c
      include 'stdio.inc'
c
      character hlptxt*132,line*80, cat1*40, cat2*40, hlpcat*40
      integer ibar,iblk,lth1,lth2,nlines,LINES
      logical found1,found2,kywrd1,kywrd2,match1,noncmd
      parameter(LINES=23)
c
      call gettkn( line, cat1, lth1 )
      call gettkn( line, cat2, lth2 )
      if (lth1.ne.0) call touppr(cat1,lth1) 
      if (lth2.ne.0) call touppr(cat2,lth2) 
      open (ludisk,file='/afs/msc/home/freed/sanghyuk/bin/nlshlp.txt',
     #     status='old',access='sequential',err=10)
c
      found1=.false.
      found2=.false.
      nlines=0
c
c  Search through lines in the help text file
c
 1    read (ludisk,'(a)',end=10,err=11) hlptxt
      ibar=1
 2    if (hlptxt(ibar:ibar).ne.'|'.and.ibar.lt.132) then
         ibar=ibar+1
         go to 2
      end if
c
      kywrd1=hlptxt(1:1).eq.'*'
      kywrd2=hlptxt(1:1).eq.'>'
      noncmd=hlptxt(1:1).eq.' '
      hlpcat=hlptxt(2:ibar-1)
      hlptxt=hlptxt(ibar+1:)
      ibar=ibar-2
c
c     Find last nonblank character in help text line
c
      iblk=132
 3    if (hlptxt(iblk:iblk).eq.' ') then
         iblk=iblk-1
         go to 3
      end if
c
c     If help text entry represents a major category, check the
c     first keyword specified in the help command (if any) against it 
c
      if (kywrd1) then
         match1=.false.
         if (lth1.eq.0) then
            call linchk(nlines)
            write(luttyo,1000) hlpcat(:ibar),hlptxt(:iblk)
         else if (cat1(:lth1).eq.hlpcat(:lth1)) then 
            call linchk(nlines)
            write(luttyo,1000) hlpcat(:ibar),hlptxt(:iblk)
            match1=.true.
            found1=.true.
         end if
c
c     If help text entry represents a subcategory, check the
c     second keyword specified in the help command (if any) against it 

      else if (kywrd2) then
         if (match1.and.lth2.eq.0) then
            call linchk(nlines)
            write(luttyo,1004) hlpcat(:ibar),hlptxt(:iblk)
         else if (match1.and.cat2(:lth2).eq.hlpcat(:lth2)) then
            call linchk(nlines)
            write(luttyo,1004) hlpcat(:ibar),hlptxt(:iblk)
            found2=.true.
         end if
c
      else if (noncmd .and. lth1.ne.0) then
         if (cat1(:lth1).eq.hlpcat(:lth1)) then
            call linchk(nlines)
            write(luttyo,1000) hlpcat(:ibar),hlptxt(:iblk)
            found1=.true.
         end if
      end if
      go to 1
c
 10   if ((lth1.ne.0 .and. .not.found1) .or. 
     #    (lth2.ne.0 .and. .not.found2)) 
     #     write (luttyo,1001) cat1(:lth1),cat2(:lth2)
c
      close(ludisk)
      write (luttyo,*)
      return
c
 11   write (luttyo,1003) hlpcat(:ibar),hlptxt(:iblk)
      close(ludisk)
      write (luttyo,*)
      return
c
 1000 format(a,t22,a)
 1001 format('*** No help available for ''',a,' ',a,''' ***')
 1002 format('*** File ''nlshlp.txt'' not available ***')
 1003 format('*** Error reading file ''nlshlp.txt'' ***')
 1004 format(2x,a,t24,a)
      end


      subroutine linchk( nlines )
      implicit none
      include 'stdio.inc'
c
      integer nlines,MXLINES
      character dummy*1
      parameter(MXLINES=20)
c
      if (nlines.eq.0) write (luttyo,1001)
      nlines=nlines+1
      if (nlines.gt.MXLINES) then
         write (luttyo,1000)
         read (luttyo,'(a)') dummy
         nlines=1
      end if
      return
 1000 format('...press <RETURN> to continue...') 
 1001 format(/15x,' *** NLSPMC on-line help ***')
      end
