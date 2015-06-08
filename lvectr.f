c Version 1.6  5/2/94
c*********************************************************************
c
c                       ================
c                       PROGRAM : LVECTR
c                       ================
c
c       This program will calculate the starting vector for EPRLL
c
c       written by DJS 29-AUG-87 (originally program STVT)
c       modified by DEB 15-MAR-92 to use new filename convention
c
c       Uses :
c               rddat.f
c               stvect.f
c               getids.f
c               setnam.f
c
c*********************************************************************
c
      program lvectr
c
      implicit none
c
      include 'stddim.inc'
      include 'stdio.inc'
      include 'eprdat.inc'
      include 'indexl.inc'
      include 'fnames.inc'
      include 'rndoff.inc'
c
      integer i,ierr,ifile,nfiles,lr,kr,krmx,mr,mrmx,
     #     ipnr,ipnrmx,ipnrmn,iqnr,iqnrmx,iqnrmn
      logical fexist
c
      double precision x
      dimension x(2,mxdim)
c
      integer ipar
      external ipar
c
      include 'version.inc'
c
c#####################################################################
c
      write (luttyo,1000)
      write (luttyo,1010) version,vdate
c
      call getids(nfiles,mxcalc)
      if (nfiles.eq.0) goto 9999
c
      do 200 ifile=1,nfiles
c
      fileid=files(ifile)
      call setnam
c
c----------------------------------------------------------------------
c     read in parameters from file
c----------------------------------------------------------------------
c
        inquire(file=prname(:namlth),exist=fexist)
c
        if (fexist) then
          call rddat(prname,namlth,ierr)
          if (ierr.eq.-100) then
            write (luttyo,1201) prname(:namlth),version
            write (lulog,1201) prname(:namlth),version
            go to 200
          else if (ierr.ne.0) then
            write (luttyo,1200) prname(:namlth)
            go to 200
          end if
        else
          write (luttyo,1100) prname(:namlth)
          go to 200
        endif
c
      call lbasix( ixname,namlth,0 )
c
c----------------------------------------------------------------------
c     starting vector calculation
c----------------------------------------------------------------------
c
      call stvect( x )
c
      write (luttyo,1070) nelv
      write (luttyo,1050)

      do 100 i=1,ndim
       if(dabs(x(1,i)).gt.rndoff) write(luttyo,1060) i,x(1,i),x(2,i),
     #                           l1(i),k1(i),m1(i),pi1(i),qi1(i)
 100   continue
c
      write(luttyo,1080) ndim
c
 200  continue
c
 9999 write(luttyo,1000)
      stop
c
c======================================================================
c     format statements
c======================================================================
c
 1000 format(2x,70('#'))
 1010 format(25x,'program LVECTR'/22x,'Version ',a,1x,a/
     #       22x,'----------------')
 1050 format (/,3x,'Row',3x,'Real Part',6x,'Imaginary Part',6x,
     #        'L  K   M   pI  qI',/)
 1060 format (1x,i5,' (',g14.7,',',g14.7,')   ',
     # '|',4(i3,','),i3,'>')
 1070 format(2x,'Non-zero elements in starting vector : ',i5)
 1080 format(2x,'Matrix dimension : ',i5)
 1100 format(20x,'*** Cannot find file ',a,' ****')
 1200 format(10x,'*** Error reading file ',a,' ***')
 1201 format(15x,'*** ''',a,''' is not an EPRLL Version ',a,' file',
     #' ***')
c
c=====================================================================
c
      end
