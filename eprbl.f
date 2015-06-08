c Version 1.6  8/12/94
c----------------------------------------------------------------------
c
c                            ===============
c                            PROGRAM : EPRBL
c                            ===============
c
c       This program does the EPR "field sweep" conjugate gradients
c       calculation in order to determine the maximum projection of each 
c       basis vector on the solution vector across the EPR spectrum.
c       The program stores this information in a basis set file named
c       <filename>.bss. Each projection is normalized to the largest
c       projection found in the basis set.
c
c       The basis set projections are then used to determine the minimum 
c       truncation set (MTS) indices as originally done by the program
c       TNLL in Versions 1.0 through 1.3. This produces an MTS database
c       file ( <filename>.mts ) containing the MTS parameters as a function
c       of "significance factor", i.e., the minimum projection for which
c       a basis element is considered to have a "significant" overlap with
c       the solution vector. Sixteen different MTS sets are specified for
c       significance factors of 2**-(i+1) for i from 1 to 16. 
c
c       The basis set projections are also used to determine a "pruned"
c       basis set that is obtained by removing all basis vectors whose
c       maximum projection on the solution vector across the spectrum is
c       less than a specified cutoff value, btol. This procedure may often 
c       produce basis sets that are considerably smaller than those obtained
c       using the simple MTS truncation rules. If btol is zero, then no
c       pruning will take place. 
c
c       In the case that the input parameter file specifies a MOMD calculation,
c       the MTS database file is not produced; however the MTS and pruned 
c       basis set are determined for the entire series of tilt angles.
c
c       Written by SHL based on original EPRBL, TNLL programs from DJS
c       (through version 1.3). Modified by DEB 16-APR-93 for EPRLL calculation.
c
c       Uses :
c               getids.f
c               setnam.f
c               rddat.f
c               wrdat.f
c               lbasix.f
c               matrll.f
c               stvect.f
c               cspccg.f
c               zdotu.f
c               zaxpy.f
c               scmvm.f
c               dtime.c
c               mtsfmt (coded in this file)
c
c----------------------------------------------------------------------
c
      program eprbl
c
      implicit none
c
      include 'stddim.inc'
      include 'stdio.inc'
      include 'eprdat.inc'
      include 'eprmat.inc'
      include 'spectr.inc'
      include 'baswt.inc'
      include 'pidef.inc'
      include 'fnames.inc'
      include 'indexl.inc'
      include 'maxl.inc'
      include 'timer.inc'
c
      integer itfld
      common /fldtab/ itfld(MXDIM)
c
      double precision b
      common /bvec/ b(2,MXDIM)
c
      integer MXCUT
      parameter (MXCUT=16)
c
      double precision ZERO
      parameter (ZERO=0.0D0)
c
      logical fexist,fmomd,report
c
      integer i,ierr,iort,maxle,maxlo,maxk,maxm,maxpi,minle,minlo,
     #     mink,minm,minpi,ncalc,nfiles,nbas,nmts,nprune,nstpmx,
     #     ntotal
c
      real cpufs,cpumat,cpuvec,cputot
      double precision bcut,bmax,cspsi,dfield,terrmx,tmp
      character mtsln*76
c
      character*9 mtsfmt
      external mtsfmt
c
      include 'version.inc'
c
c#####################################################################
c
      write (luttyo,1000)
      write (luttyo,1010) version,vdate
      write (luttyo,1000)
c
      call getids(nfiles,MXCALC)
c
c======================================================================
c   Loop over calculations
c======================================================================
c
      cputot=dtime( tarray )
      do 100 ncalc=1,nfiles
c
c---------------------------------------------------------------------
c     construct file names
c---------------------------------------------------------------------
c
         fileid = files(ncalc)
c
         call setnam
         call setflg( flags(ncalc) )
c
c----------------------------------------------------------------------
c     open log file and write message
c----------------------------------------------------------------------
c
         open (unit=lulog,file=rlname(:namlth),status='unknown',
     #      access='sequential',form='formatted')
         rewind (unit=lulog)
c
         write (lulog,1000)
         write (luttyo,1030) prname(:namlth)
         write (lulog,1030) prname(:namlth)
         write (lulog,1000)
c
c----------------------------------------------------------------------
c     read parameters from file
c----------------------------------------------------------------------
c
         inquire(file=prname(:namlth),exist=fexist)
c
         if (fexist) then
            call rddat(prname,namlth,ierr)
           if (ierr.eq.-100) then
              write (luttyo,1201) prname(:namlth),version
              write (lulog,1201) prname(:namlth),version
              go to 100
           else if (ierr.ne.0) then
              write (luttyo,1200) prname(:namlth)
              write (lulog,1200) prname(:namlth)
              go to 100
           end if
         else
            write (luttyo,1100) prname(:namlth)
            write (lulog,1100) prname(:namlth)
            go to 100
         endif
c
c----------------------------------------------------------------------
c     check if a field-swept calculation is specified
c----------------------------------------------------------------------
c
         if (itype.ne.2) then
            write (luttyo,1020) prname(:namlth)
            write (lulog,1020) prname(:namlth)
            go to 999
         end if
c
c      ----------------------------------------------------------
c       Set up calculation by determining initial basis indices, 
c       and setting the starting vector
c      ----------------------------------------------------------
c
         call lbasix(ixname,namlth,1)
c
         if (ndim.gt.MXDIM) then
            write (luttyo,1060) ndim,MXDIM
            write (lulog,1060) ndim,MXDIM
            go to 999
         end if          
c
         write (luttyo,2000) 
         write (lulog,2000)
c
         call stvect( b )
         cpuvec=dtime( tarray )
c
         write (luttyo,2010) nelv,cpuvec
         write (lulog,2010) nelv,cpuvec
c
c
c    -------------------------
c     Zero out spectrum array
c    -------------------------
c     
         do i=1,MXPT
            spect(1,i)=ZERO
            spect(2,i)=ZERO
         end do
c
c    -----------------------------
c     Zero basis projection array 
c    -----------------------------
c
         nprune=0
         do i=1,MXDIM
            itfld(i)=0
            keep(i)=0
            totwt(i)=ZERO
         end do
c
c----------------------------------------------------------------------
c   Determine basis projections for calculation (or first orientation
c   in a MOMD calculation) using swept-field conjugate gradients
c----------------------------------------------------------------------
c
         fmomd=nort.gt.1
         report=.not.fmomd
c
c       --------------------------------------------------
c       Set director tilt angle psi in MOMD calculations
c       --------------------------------------------------
         do iort=1,max(nort,1)
            if (fmomd) then
               if (iort.eq.1) then
                  write (luttyo,2100)
                  write (lulog,2100)
               end if
               cspsi=dfloat(iort-1)/dfloat(nort-1)
               psi=dacos(cspsi)*1.8D2/PI
            else
               write(luttyo,1040)
               write(lulog,1040) 
            end if
c
            call matrll( ierr )
            cpumat=dtime( tarray )
c
            if (ierr.ne.0) then
               write (luttyo,1070) cpumat
               write (lulog,1070) cpumat
               if (ierr.eq.1) then
                  write (luttyo,1080) MXEL
                  write (lulog,1080) MXEL
               end if
               goto 999
            end if
c
            bmax=ZERO
            call eprfsl( b,ixfld,bsswt,spect,bmax,terrmx,nstpmx,
     #                   cpufs,ntotal,report )
            cputot=cputot+cpufs
c
c            -----------------------------------------
c            Determine MTS for the current calculation
c            -----------------------------------------
c
            call mtsl( bsswt,ndim,btol*bmax,maxle,minle,maxlo,minlo,
     #                 maxk,mink,maxm,minm,maxpi,minpi,nbas,nmts )
c
c             ---------------------------------------
c              Add elements to the pruned basis set 
c             ---------------------------------------
c
            do i=1,ndim
               if (bsswt(i).ge.btol*bmax) then
                  if (keep(i).eq.0) nprune=nprune+1
                  keep(i)=1
               end if
               tmp=bsswt(i)/bmax
               if (tmp.gt.totwt(i)) then
                  totwt(i)=tmp
                  itfld(i)=ixfld(i)
               end if
            end do
c
c           ------------------------------------------------------- 
c           For MOMD calculations, write out current angle, unpruned
c           dimension, pruned dimension, and dimension of combined
c           pruned sets so far
c           ------------------------------------------------------- 
c
            if (fmomd) then
               write (luttyo,2110) iort,psi,neltot,ndim,nbas,
     #                             nprune,cpufs
               write (lulog,2110) iort,psi,neltot,ndim,nbas,
     #                            nprune,cpufs
            end if
c
         end do
c                    *** End of loop over fields
c
c       ----------------------------------------------------
c        Write out basis set projections in file <name>.bss
c       ----------------------------------------------------
c
         dfield=(fieldf-fieldi)/dble(nfield-1)
         write (luttyo,5000) bsname(:namlth)
         open (unit=ludisk,file=bsname(:namlth),status='unknown',
     #        access='sequential',form='unformatted')
         write (ludisk) ndim,lemx,lomx,kmn,kmx,mmn,mmx,ipnmx,in2,
     #                  ldelta,kdelta,ipsi0,jkmn,jmmn
         write (ludisk) (b0+fieldi+dfield*(itfld(i)-1),totwt(i),
     #                  i=1,ndim)
         close (unit=ludisk)
c
c       ----------------------------------------
c        Write out spectrum in file <name>.spc
c       ----------------------------------------
         write (luttyo,5000) spname(:namlth)
         open (unit=ludisk,file=spname(:namlth),status='unknown',
     #        access='sequential',form='formatted')
         write (ludisk,3200) (b0+fieldi+(i-1)*dfield,spect(1,i),
     #                        spect(2,i),i=1,nfield)
         close (unit=ludisk)
c
c       -------------------------------------------------------
c        Determine the MTS as a function of significance factor
c        and write the MTS database file <name>.mts 
c       -------------------------------------------------------
c
         write (luttyo,5000) tsname(:namlth)
         open (unit=ludisk,file=tsname(:namlth),status='unknown',
     #         access='sequential',form='formatted')
         write(ludisk,4200)
         write(ludisk,4208)
c
         do i=1,MXCUT
            bcut = 0.5d0**(i-1)
            nmts=-1
            call mtsl( totwt,ndim,bcut*bmax,maxle,minle,maxlo,minlo,
     #                 maxk,mink,maxm,minm,maxpi,minpi,nbas,nmts )
            write(mtsln,4209) bcut,nmts,nbas
            if (nmts.gt.0) then
               mtsln(13:21)=mtsfmt(minle,maxle,MXLVAL+1)
               mtsln(23:31)=mtsfmt(minlo,maxlo,MXLVAL+1)
               mtsln(33:41)=mtsfmt(mink,maxk,MXLVAL+1)
               mtsln(43:51)=mtsfmt(minm,maxm,MXLVAL+1)
               mtsln(53:61)=mtsfmt(minpi,maxpi,in2+1)
            end if
            write (ludisk,'(a)') mtsln
         end do
         write (ludisk,4210)
         close (ludisk)
c
c       -----------------------------------------------------------
c        Condense basis projection and index arrays to retain only
c        the pruned basis 
c       -----------------------------------------------------------
c
         nprune=0
         do i=1,ndim
            if (keep(i).ne.0) then
               nprune=nprune+1
               totwt(nprune)=totwt(i)
               l1(nprune)=l1(i)
               k1(nprune)=k1(i)
               m1(nprune)=m1(i)
               pi1(nprune)=pi1(i)
               qi1(nprune)=qi1(i)
            end if
         end do
c
         write (luttyo,4000) ndim,nprune
         write (lulog,4000) ndim,nprune
c
c       --------------------------------------------
c        Find MTS corresponding to pruned basis set 
c       --------------------------------------------
c
         ndim=nprune
         nmts=-1
         call mtsl( totwt,ndim,ZERO,maxle,minle,maxlo,minlo,maxk,mink,
     #              maxm,minm,maxpi,minpi,nbas,nmts )
c
         write(luttyo,4200)
	 write(luttyo,4208)
         write(mtsln,4209) btol,nmts,nbas
	 if (nmts.gt.0) then
            mtsln(13:21)=mtsfmt(minle,maxle,MXLVAL+1)
            mtsln(23:31)=mtsfmt(minlo,maxlo,MXLVAL+1)
            mtsln(33:41)=mtsfmt(mink,maxk,MXLVAL+1)
            mtsln(43:51)=mtsfmt(minm,maxm,MXLVAL+1)
            mtsln(53:61)=mtsfmt(minpi,maxpi,in2+1)
         end if
         write(luttyo,'(a)') mtsln
         write(luttyo,4210)
c
         write(lulog,4200)
	 write(lulog,4208)
	 write(mtsln,4209) btol,nmts,nbas
         if (nmts.gt.0) then
            mtsln(13:21)=mtsfmt(minle,maxle,MXLVAL)
            mtsln(23:31)=mtsfmt(minlo,maxlo,MXLVAL)
            mtsln(33:41)=mtsfmt(mink,maxk,MXLVAL)
            mtsln(43:51)=mtsfmt(minm,maxm,MXLVAL)
            mtsln(53:61)=mtsfmt(minpi,maxpi,in2+1)
         end if
         write(lulog,'(a)') mtsln
         write(lulog,4210)
c
c       ----------------------------------------------------
c        Output updated parameter file and basis index file
c       ----------------------------------------------------
c
c
         lemx=maxle
         lomx=max0(maxlo,0)
         kmx=maxk
         mmx=maxm
c -----------------------------------
c Two lines needed to reset kmn and mmn if they change.  -JPB April 20 1998
c -----------------------------------
         kmn=mink
         mmn=minm
c -----------------------------------
         ipnmx=maxpi
c
         write (luttyo,5000) prname(:namlth)
         call wrdat(prname,namlth)
c     
c        -----------------------------------------------
c        Write out basis set indices in <name>.ind file
c        -----------------------------------------------
c
         write (luttyo,5000) ixname(:namlth)
         open ( unit=ludisk,file=ixname(:namlth),status='unknown',
     #          access='sequential',form='unformatted' )
         write (ludisk) ndim,lemx,lomx,kmn,kmx,mmn,mmx,ipnmx
         write (ludisk) (l1(i),k1(i),m1(i),pi1(i),qi1(i),i=1,ndim)
         close (unit=ludisk)
c
         write (luttyo,5100) cputot
         write (lulog,5100) cputot
c     
         write(luttyo,1000)
         write(lulog,1000)	
c
c----------------------------------------------------------------------
c       End of loop over input files
c----------------------------------------------------------------------
 100  continue
c
 999  close (unit=lulog)
      stop
c
c=====================================================================
c
 1000 format(2x,70('#'))
 1010 format(25x,'program EPRBL'/22x,'Version ',a,1x,a/
     #       22x,'----------------')
c
 1020 format(10x,
     #     '*** ''',a,''' is not set for EPRBL (itype=2) ***')
 1030 format(25x,'file : ',a)
 1040 format(/15x,'*** Calculating matrix elements ***')
 1050 format(2x,'Dimension of the matrix : ',i6,
     #     /,2x,'Total number of non-zero entries : ',i8,
     #     /,2x,'CPU time (seconds) : ',f10.2)
 1060 format(20x,'*** Basis set is too large ***',
     #       20x, 'Specified dimension=',i7,'; maximum =',i7)
 1070 format(10x,
     #     '*** Serious error encountered in matrix generation ***'
     #     /20x,'CPU time (seconds) : ',f10.2)
 1080 format(10x,
     #  '*** Number of off-diagonal matrix elements exceeds ',i8,' ***')
 1090 format(20x,'*** Matrix dimension exceeds ',i6,' ***')
c
 1100 format(10x,'*** file ',a,' not found ***')
 1200 format(10x,'*** Error reading file ''',a,''' ***')
 1201 format(10x,'*** ''',a,''' is not an EPRLL Version ',a,' file',
     #' ***')
 2000 format(/15x,'*** Calculating starting vector ***')
 2010 format(2x,'Number of non-zero elements in ',
     #     'starting vector : ',i5,
     #     /,2x,'CPU time (seconds) : ',f10.2)
c
 2100 format(/,2x,70('-'),/,4x,'i',3x,'angle',4x,'neltot',3x,
     #      'updim',3x,'pdim',3x,'totdim',4x,'cpu'/,2x,70('-'))
 2110 format(2x,i3,3x,f5.2,3x,i7,3x,i4,3x,i4,5x,i4,3x,f6.2)
 3200 format(1x,3(g14.7,','))
 4000 format(/,2x,'unpruned dimension = ',i4,
     #           5x,'pruned dimension = ',i5)
 4200 format(/,2x,74('-')/
     #       3x,'error',7x,'Le',8x,'Lo',8x,'K',9x,'M',9x,'pI',6x,'MTS',
     #       4x,'pruned'/,12x,5('min  max  '),' dim',5x,' dim')
 4208 format(2x,6('---------+'),'-------+','------')
 4209 format(1x,g10.3,5('|',9x),'|',i6,' |',i6)
 4210 format(2x,6('---------+'),'-------+','------'/)
 5000 format(13x,'*** writing file ',a,' ***')
 5100 format(2x,'Total CPU time for EPRBL calculation (seconds) : ',
     #     f10.2)
c
c=====================================================================
c
      end

      function mtsfmt(ixmin,ixmax,maxabs)
      implicit none
      character*9 mtsfmt
      integer ixmin,ixmax,maxabs
c
      mtsfmt=' '
      if (ixmin.ge.abs(maxabs) .or. ixmax.le.-abs(maxabs) .or.
     #    ixmin.gt.ixmax) then
         mtsfmt='  --   --'
      else
         write(mtsfmt(1:4),'(i4)') ixmin
         write(mtsfmt(5:9),'(i4)') ixmax
      endif
      return
      end

