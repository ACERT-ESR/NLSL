c Version 1.6  5/2/94
c**********************************************************************
c
c                       =================
c                       PROGRAM : PRUNEL
c                       =================
c
c       This program reads in the basis set truncation file <name>.bss 
c       written by EPRBL, and creates an output file containing a
c       list of basis vector indices determined using a prespecified
c       cutoff tolerance.
c
c       Uses:
c               ipar.f
c               getids.f
c               setnam.f
c               setflg.f
c               lbasix.f
c               rddat.f
c
c**********************************************************************
c
      program prunel
c
      implicit none
c
      include 'stdio.inc'
      include 'eprdat.inc'
      include 'stddim.inc'
      include 'fnames.inc'
      include 'indexl.inc'
      include 'rndoff.inc'
c
      integer i,ierr,ifile,ipar,j,nbss,ndim1,nfiles
      double precision bsswt(MXDIM),bmax,dmy,ftoken
      logical fexist,fprmpt,prmpt
      character line*80
c
      external ipar,prmpt,ftoken
c
      include 'version.inc'
c
c######################################################################
c
      write(luttyo,1000)
      write(luttyo,1010) version,vdate
c
      ifile=0
      call getids(nfiles,MXCALC)
      if (nfiles.eq.0) goto 9999
c
 1    ifile=ifile+1
      if (ifile.gt.nfiles) goto 9999
      fileid=files(ifile)
      call setnam
      call setflg( flags(i) )
c
c----------------------------------------------------------------------
c     read input parameters from input files
c----------------------------------------------------------------------
c
        inquire(file=prname,exist=fexist)
c
        if (fexist) then
          call rddat(prname,namlth,ierr)
          if (ierr.eq.-100) then
            write (luttyo,1201) prname(:namlth),version
            write (lulog,1201) prname(:namlth),version
            go to 1
          else if (ierr.ne.0) then
            write (luttyo,1200) prname(:namlth)
            write (lulog,1200) prname(:namlth)
            go to 1
          end if
        else
          write (luttyo,1100) prname(:namlth)
          write (lulog,1100) prname(:namlth)
          go to 1
        endif
c
c----------------------------------------------------------------------
c       Read in basis weights from .bss file
c----------------------------------------------------------------------
        write (luttyo,5000) ifile,bsname(:namlth)
c
      open (unit=ludisk,file=bsname,status='old',
     #     access='sequential',form='unformatted',iostat=ierr)
      if (ierr.ne.0) then
         write (luttyo,1100) bsname
         write (lulog,1100) bsname
	 goto 1
      end if
      read (ludisk,iostat=ierr) ndim,lemx,lomx,kmn,kmx,mmn,mmx,ipnmx,
     #                          in2,ldelta,kdelta,ipsi0,jmmn,jkmn
      if (ierr.eq.0) read (ludisk,iostat=ierr) (dmy,bsswt(i),i=1,ndim)
      close (unit=ludisk)
      if (ierr.ne.0) then
         write (luttyo,1200) bsname
         write (lulog,1200) bsname
	 goto 1
      end if
c
c----------------------------------------------------------------------
c     Build the full list of basis indices corresponding to the MTS
c----------------------------------------------------------------------
c
      ndim1=ndim
      call lbasix( ixname,namlth,1 )
      if (ndim.ne.ndim1) write(luttyo,6000) ndim,bsname(:namlth),ndim1
c
      j=0
c
      write (luttyo,5030) bsname(:namlth-4),btol
      fprmpt=prmpt( line )
      if (fprmpt) btol=ftoken( line )
      if (btol.lt.0.0D0) btol=0.0D0
      if (btol.gt.0.5D0) btol=0.5D0
c
c----------------------------------------------------------------------
c     Prune the basis set using the new tolerance
c----------------------------------------------------------------------
c
      lemx=0
      lomx=0
      kmx=0
      kmn=0
      mmx=0
      mmn=0
      ipnmx=0
      do i=1,ndim
         if (bsswt(i).ge.btol) then
            j=j+1
            l1(j)=l1(i)
            if (ipar(l1(j)).eq.1) then
               if (l1(j).gt.lemx) lemx=l1(j)
            else
               if (l1(j).gt.lomx) lomx=l1(j)
            endif 
c
            k1(j)=k1(i)
            if (k1(j).gt.kmx) kmx=k1(j)
            if (k1(j).lt.kmn) kmn=k1(j)
c
            m1(j)=m1(i)
            if (m1(j).gt.mmx) mmx=m1(j)
            if (m1(j).lt.mmn) mmn=m1(j)
c
            pi1(j)=pi1(i)
            if (abs(pi1(j)).gt.ipnmx) ipnmx=abs(pi1(j))
            qi1(j)=qi1(i)
         end if
      end do
c 
      write (luttyo,5020) ndim,j,btol,lemx,lomx,kmn,kmx,mmn,mmx,ipnmx
      ndim=j
c
c----------------------------------------------------------------------     
c   Write out basis set indices for pruned set
c----------------------------------------------------------------------
c
      write (luttyo,5010) ixname(:namlth)
      open (unit=ludisk,file=ixname(:namlth),status='unknown',
     #     access='sequential',form='unformatted')
      write (ludisk) ndim,lemx,lomx,kmn,kmx,mmn,mmx,ipnmx
      write (ludisk) (l1(i),k1(i),m1(i),pi1(i),qi1(i),i=1,ndim)
      close (unit=ludisk)
c
      go to 1
c
c----------------------------------------------------------------------
c     exit from program
c----------------------------------------------------------------------
c
 9999 write(luttyo,1000)
      stop
c
c======================================================================
c    format statements
c======================================================================
c
 1000 format(2x,70('#'))
 1010 format(25x,'program PRUNEL'/22x,'Version ',a,1x,a/
     #       22x,'----------------')
 1100 format(20x,'*** Cannot find file ',a,' ****')
 1200 format(10x,'*** Error reading file ',a,' ***')
 1201 format(10x,'*** ''',a,''' is not an EPRLL Version ',a,' file',
     #' ***')
 5000 format(15x,38('='),/,
     #     15x,'Pruning procedure ',i2,' using file ',a,/,
     #     15x,38('='))
 5010 format(20x,'*** writing file ',a,' ***')
 5020 format(20x,'MTS basis dimension: ',i4/
     #       20x,'pruned dimension:    ',i4/
     #       20x,'pruning tolerance :  ',f8.4/
     #       20x,'Trunc. indices for pruned set:'/
     #       20x,'lemx=',i3,' lomx=',i3,' kmn=',i3,' kmx=',i3/
     #       20x,' mmn=',i3,' mmx=',i3,' ipnmx=',i2/ )
 5030 format(2x,'Pruning tolerance for file ''',a,''' [btol] ',g9.3,
     #     ' : ')
 6000 format(10x,'*** Dimension discrepancy: ndim=',i6,'; Dim in file ',
     #a,'=',i6)
c

c======================================================================
c
      end
