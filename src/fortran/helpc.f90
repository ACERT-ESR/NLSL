! Version 1.2
!----------------------------------------------------------------------
!                  =========================
!                      subroutine HELPC
!                  =========================
!----------------------------------------------------------------------
      subroutine helpc(line)
!
      use stdio
!
      implicit none
      character hlptxt*132,line*80, cat1*30, cat2*30, hlpcat*30
      integer ioerr,ibar,iblk,lth1,lth2,nlines,LINES
      logical found1,found2,kywrd1,kywrd2,match1,noncmd
      parameter(LINES=23)
!
      ibar=1
      iblk=1
      hlpcat=''
      match1=.false.

      call gettkn( line, cat1, lth1 )
      call gettkn( line, cat2, lth2 )
      if (lth1.ne.0) call touppr(cat1,lth1) 
      if (lth2.ne.0) call touppr(cat2,lth2) 
      open (ludisk,file='nlshlp.txt',
     #     status='old',access='sequential',iostat=ioerr)
      if (ioerr.ne.0) then
         write (luttyo,1002)
         return
      end if
!
      found1=.false.
      found2=.false.
      nlines=0
!
!  Search through lines in the help text file
!
 1    read (ludisk,'(a)',end=4,iostat=ioerr) hlptxt
      if (ioerr.ne.0) then
         write (luttyo,1003) hlpcat(:ibar),hlptxt(:iblk)
         close(ludisk)
         return      
      end if
      ibar=1
 2    if (hlptxt(ibar:ibar).ne.'|'.and.ibar.lt.132) then
         ibar=ibar+1
         go to 2
      end if
!
      kywrd1=hlptxt(1:1).eq.'*'
      kywrd2=hlptxt(1:1).eq.'>'
      noncmd=hlptxt(1:1).eq.' '
      hlpcat=hlptxt(2:ibar-1)
      hlptxt=hlptxt(ibar+1:)
      ibar=ibar-2
!
!     Find last nonblank character in help text line
!
      iblk=132
 3    if (hlptxt(iblk:iblk).eq.' ') then
         iblk=iblk-1
         go to 3
      end if
!
!     If help text entry represents a major category, check the
!     first keyword specified in the help command (if any) against it 
!
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
!
!     If help text entry represents a subcategory, check the
!     second keyword specified in the help command (if any) against it 

      else if (kywrd2) then
         if (match1.and.lth2.eq.0) then
            call linchk(nlines)
            write(luttyo,1004) hlpcat(:ibar),hlptxt(:iblk)
         else if (match1.and.cat2(:lth2).eq.hlpcat(:lth2)) then
            call linchk(nlines)
            write(luttyo,1004) hlpcat(:ibar),hlptxt(:iblk)
            found2=.true.
         end if
!
      else if (noncmd .and. lth1.ne.0) then
         if (cat1(:lth1).eq.hlpcat(:lth1)) then
            call linchk(nlines)
            write(luttyo,1000) hlpcat(:ibar),hlptxt(:iblk)
            found1=.true.
         end if
      end if
      go to 1
!
 4    if ((lth1.ne.0 .and. .not.found1) .or. 
     #    (lth2.ne.0 .and. .not.found2)) 
     #     write (luttyo,1001) cat1(:lth1),cat2(:lth2)
!
      close(ludisk)
      return
!
 1000 format(a,t22,a)
 1001 format('*** No help available for ''',a,' ',a,''' ***')
 1002 format('*** File ''nlshlp.txt'' not available ***')
 1003 format('*** Error reading file ''nlshlp.txt'' ***')
 1004 format(2x,a,t24,a)
      end


      subroutine linchk( nlines )
!
      use stdio
!
      implicit none
      integer nlines,MXLINES
      character dummy*1
      parameter(MXLINES=20)
!
      if (luttyo.ne.luttyo) return
      if (nlines.eq.0) write (luttyo,1001)
      nlines=nlines+1
      if (nlines.gt.MXLINES) then
         write (luttyo,1000)
         read (luttyo,'(a)') dummy
         nlines=1
      end if
      return
 1000 format('...press <RETURN> to continue...') 
 1001 format(/15x,' *** NLSL on-line help ***'/)
      end
