c Version 1.6  5/2/94
c**********************************************************************
c
c                       ================
c                        program LMATRX
c                       ================
c
c       This program will read in the matrix data file written by
c       the main calculation programs EPRCGL or EPRL and then 
c       write it to the screen in an understandable form.
c
c       written by DEB 19-MAR-92 based on MATLST by DJS 29-SEP-87
c       Bug fixed by SL on 9-DEC-92: handles case of an empty row
c
c       Includes:
c               stddim.inc
c               stdio.inc
c               eprmat.inc
c               eprdat.inc
c               indexl.inc
c               fnames.inc
c
c       Uses:
c               rddat.f
c               setnam.f
c               lbasix.f
c               getids.f
c
c**********************************************************************
c
      program lmatrx
c
      implicit none
c
      include 'stddim.inc'
      include 'stdio.inc'
      include 'eprmat.inc'
      include 'eprdat.inc'
      include 'indexl.inc'
      include 'fnames.inc'
      include 'rndoff.inc'
c
      integer i,ierr,j,jmax,k,kmax,m,n,ifile,nfiles
      double precision zr,zi
      logical fexist
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
c---------------------------------------------------------------------
c     get file identifier and construct file names
c---------------------------------------------------------------------
c
         fileid=files(ifile)
         call setnam
c
c---------------------------------------------------------------------
c     read in data from parameter file
c     and build list of basis indices
c---------------------------------------------------------------------
c
        inquire(file=prname,exist=fexist)
c
        if (fexist) then
          call rddat(prname,namlth,ierr)
          if (ierr.eq.-100) then
            write (luttyo,1201) prname(:namlth),version
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
c
      call lbasix( ixname,namlth,0 )
c       
c---------------------------------------------------------------------
c     read in matrix
c---------------------------------------------------------------------
c
      write (luttyo,1110) mtname(:namlth)
c
      open (unit=ludisk,file=mtname(:namlth),status='old',
     #     access='sequential',form='unformatted')
c
      read (ludisk) (zdiag(1,i),zdiag(2,i),i=1,ndim)
      read (ludisk) (jzmat(i),i=1,ndim+1) 
      read (ludisk) (izmat(i),zmat(i),i=1,nelim)
      read (ludisk) (kzmat(i),i=1,ndim+1)
      read (ludisk) (izmat(MXEL-i+1),zmat(MXEL-i+1),i=1,nelre)
c
      close (unit=ludisk)
c
      write (luttyo,2002) ndim,nelim,nelre,nelim+nelre
c
c---------------------------------------------------------------------
c     write out upper half of matrix
c---------------------------------------------------------------------
c
      open (unit=lulog,file=mfname(:namlth),status='unknown',
     #     access='sequential',form='formatted')
      write (luttyo,1120) mfname(:namlth)
c
      write (lulog,1000)
      write (lulog,1040) mfname(:namlth)
      write (lulog,2002) ndim,nelim,nelre,nelim+nelre
c
      j=1
      k=1
c
c   --- loop over rows
c
      do 100 m=1,ndim
c
c    --- write out diagonal element
c
        write (lulog,1050) m,l1(m),k1(m),m1(m),pi1(m),qi1(m)
        if ((abs(zdiag(1,m)).gt.RNDOFF).or.
     #      (abs(zdiag(2,m)).gt.RNDOFF)) then
          write (lulog,1060) m,zdiag(1,m),zdiag(2,m),
     #          l1(m),k1(m),m1(m),pi1(m),qi1(m)
        endif
c
c   --- If there are any real or imaginary elements in this row,
c       find indices of last real and imaginary elements in this row
c
        if (jzmat(m).lt.jzmat(m+1) .or. kzmat(m).lt.kzmat(m+1) ) then
           jmax=izmat(jzmat(m+1)-1)
           kmax=izmat(MXEL-kzmat(m+1)+2)
c
c   --- loop over columns
c
           do 110 n=m+1,max(jmax,kmax)
c
c   ---  Save any imaginary element in this column
c
              if (n.eq.izmat(j) .and. n.le.jmax) then
                 zi=zmat(j)
                 j=j+1
              else
                 zi=0.0D0
              endif
c
c   --- Save any real element in this column
c
              if (n.eq.izmat(MXEL-k+1) .and. n.le.kmax) then
                 zr=zmat(MXEL-k+1)     
                 k=k+1
              else
                 zr=0.0D0
              endif
c
c  --- Output the entry if there was a real or imaginary part
c
              if ((abs(zr).gt.RNDOFF).or.(abs(zi).gt.RNDOFF)) then
                 write (lulog,1060) n,zr,zi,l1(n),k1(n),m1(n),
     #                              pi1(n),qi1(n)
              end if
c
 110       continue
c                   --- end if jzmat(m).ne.jzm(m+1) 
c                         .or. kzmat(m).ne.kzmat(m+1)
        end if
c
      if ((j.ne.jzmat(m+1)).or.(k.ne.kzmat(m+1))) then
        write (luttyo,2000) mtname(:namlth)
        write (lulog,2000) mtname(:namlth)
        go to 200
      endif
c
 100  continue
c
      write (lulog,1000) 
      close (unit=lulog)        
c
 200  continue
c---------------------------------------------------------------------
c                  exit program
c---------------------------------------------------------------------
c
 9999 write (luttyo,1000)
c
      stop
c
c=====================================================================
c     format statements
c=====================================================================
c
 1000 format (70('#'))
 1010 format(25x,'program LMATRX'/22x,'Version ',a,1x,a/
     #       22x,'----------------')
 1040 format (/10x,'matrix elements from file : ',a)
 1050 format (/,11x,'*** Row : ',i5,' ***',11x,'<',4(i3,','),i3,'|',
     #        /,3x,'Col',3x,'Real Part',6x,'Imaginary Part',6x,
     #        'L   K   M  pI  qI',/)
 1060 format (1x,i5,' (',g14.7,',',g14.7,')   ',
     # '|',4(i3,','),i3,'>')
c
 1100 format(/,15x,'*** File ''',a,''' not found ***',/)
 1110 format(/15x,' *** Reading file ''',a,''' ***')
 1120 format(/15x,' *** Writing file ''',a,''' ***')
 1200 format(/,15x,'*** Error reading file ''',a,''' ***',/)
 1201 format(/,15x,'*** ''',a,''' is not an EPRLL Version ',a,' file',
     #' ***')
 2000 format (//,8x,'*** ERROR: File ''',a,''' is corrupted ***',//)
 2002 format(/10x,'Matrix dimension: ',i6/
     #        10x,'Off-diagonal elements: ',i7,' imaginary'/
     #        33x,i7,' real'/33x,i7,' total'/)
c
c=====================================================================
c
      end 
