c Version 1.6 8/12/94
c**********************************************************************
c
c                       =================
c                         PROGRAM : BSSL
c                       =================
c
c       This program outputs a formatted record of the basis set 
c       truncation file <name>.bss written by EPRBL.
c
c       Uses:
c               getids.f
c               setnam.f
c               setflg.f
c
c**********************************************************************
c
      program bssl
      implicit none
c
      include 'stdio.inc'
      include 'stddim.inc'
      include 'eprdat.inc'
      include 'fnames.inc'
      include 'indexl.inc'
      include 'baswt.inc'
      include 'rndoff.inc'
c
      integer i,ierr,ifile,ndim1,nfiles
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
c       Read in basis weights from .bss file
c----------------------------------------------------------------------
      open (unit=ludisk,file=bsname,status='old',
     #     access='sequential',form='unformatted',iostat=ierr)
      if (ierr.ne.0) then
         write (luttyo,1100) bsname
         write (lulog,1100) bsname
	 goto 1
      end if
c
      read (ludisk,iostat=ierr) ndim,lemx,lomx,kmn,kmx,mmn,mmx,ipnmx,
     #                          in2,ldelta,kdelta,ipsi0,jmmn,jkmn
      if (ierr.eq.0) read (ludisk,iostat=ierr) 
     #               (totwt(i),bsswt(i),i=1,ndim)
      close (unit=ludisk)
      if (ierr.ne.0) then
         write (luttyo,1200) bsname
         write (lulog,1200) bsname
	 goto 1
      end if
      ndim1=ndim
c
c
c----------------------------------------------------------------------
c     Build the full list of basis indices specified by the truncation
c     indices. It should correspond to the number of elements in the bss 
c     file. Report any discrepancy.
c----------------------------------------------------------------------
c
      call lbasix( ixname,namlth,1 )
c
      if (ndim.ne.ndim1) write(luttyo,5000) ndim,bsname(:namlth),ndim1
c
c----------------------------------------------------------------------     
c   Write out basis set projections
c----------------------------------------------------------------------
c
      write (luttyo,5010) bsname(:namlth),ndim,lemx,lomx,kmn,kmx,
     #                    mmn,mmx,ipnmx
      do i=1,ndim
         write (luttyo,5011) i,l1(i),k1(i),m1(i),pi1(i),qi1(i),
     #                       totwt(i),bsswt(i)
      end do
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
 1010 format(25x,'program BSSL'/22x,'Version ',a,1x,a/
     #       22x,'----------------')
 1100 format(20x,'*** Cannot find file ',a,' ****')
 1200 format(10x,'*** Error reading file ',a,' ***')
 5000 format(10x,'*** Dimension discrepancy: ndim=',i6,'; Dim in file ',
     #a,'=',i6)
 5010 format(/20x,'File: ',a/20x,'BSS Dimension:',i6/20x,
     #'Trunc. indices:',' lemx=',i3,' lomx=',i3,' kmx=',i3,' mmx=',i3,
     #' ipnmx=',i2//9x,'Basis vector',11x,'Field(max)',5x,'Norm. Proj.'
     #/66('-'))
 5011 format(i5,2x,'|',4(i3,','),i3,'>',3x,2(2x,g12.5))
c
      end
