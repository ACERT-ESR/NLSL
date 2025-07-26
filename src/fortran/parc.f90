! NLSL Version 1.5 beta 11/25/95
!----------------------------------------------------------------------
!                    =========================
!                       subroutine PARC 
!                    =========================
!
!  Prints out a list of parameter values on the given logical unit
!  number.
!
!----------------------------------------------------------------------
      subroutine parc( line )
!
      use nlsdim
      use eprprm
      use expdat
      use parcom
      use lpnam
      use stdio
      use rnddbl
      use symdef
      use mtsdef
!
      implicit none
      character line*80,hdr*4,shdr*4,fileID*30
!
      integer i,ioerr,isym,jx,lth,lu,ntmp,npotn
!
      integer itrim
      external itrim
!
      character*2 nonbr(3)
      data nonbr/ 'l', 'xy', 'zz'/
!
!......................................................................
!
      call gettkn(line,fileID,lth)
!
!----------------------------------------------------------------------
!     No name specified: output to terminal 
!----------------------------------------------------------------------
      if (lth.eq.0) then
         lu=luttyo
         hdr=' '
         shdr=' '
      else
!
!     ------------------------------------------------
!      Set output of "let" commands to specified file 
!     ------------------------------------------------
         open(ludisk,file=fileID(:lth),status='unknown',
     #        access='sequential',form='formatted',iostat=ioerr)
         if (ioerr.ne.0) then
            write (luout,3000) fileID(:lth)
            return
         end if
         lu=ludisk
         hdr='let '
         shdr='c '
      end if
!
!     --------------------
!     Series command
!     --------------------
      if (iser.ne.0) then 
         write(lu,1023) parnam(iser)(:itrim(parnam(iser))),
     #        (serval(i),i=1,nser)
      else
         write(lu,2023) shdr
      end if
!
      ntmp=max0(nser,1)
      write(lu,1024) (hdr,jx,sb0(jx),jx=1,ntmp)
      write(lu,1021) (hdr,jx,slb(jx),jx=1,ntmp)
      write(lu,1025) (hdr,jx,spsi(jx),jx=1,ntmp)
      write(lu,1026) (hdr,jx,sphs(jx),jx=1,ntmp)
!
      write (lu,1018) nsite
!

! ========================================
!     Loop over site parameters
! ========================================
      do 10 jx=1,nsite
         write (lu,1019) shdr,jx
!
!     -----------
!      g-tensor
!     -----------
         isym = iparm(IIGFLG,jx)
         if (lu.ne.luttyo) write (lu,2002) symstr(isym), 'g'
         if (isym.eq.AXIAL) then
            write (lu,2001) hdr,(fparm(IGXX+i,jx),i=0,2,2)
         else if (isym.eq.SPHERICAL) then
            write (lu,1001) hdr,(fparm(IGXX+i,jx),i=0,2)
         else
            write (lu,1000) hdr,(fparm(IGXX+i,jx),i=0,2)
         end if
!
!    -----------------------------------
!     Nuclear spin and hyperfine tensor
!    -----------------------------------
      if (iparm(IIN2,jx).ne.0) then
         isym=iparm(IIAFLG,jx) 
         if (lu.ne.luttyo) write (lu,2002) symstr(isym), 'A' 
         if (isym.eq.AXIAL) then
            write(lu,2004) hdr,iparm(IIN2,jx),(fparm(IAXX+i,jx),i=0,2,2)
         else if (isym.eq.SPHERICAL) then
            write(lu,1004) hdr,iparm(IIN2,jx),(fparm(IAXX+i,jx),i=0,2)
         else
            write(lu,1003) hdr,iparm(IIN2,jx),(fparm(IAXX+i,jx),i=0,2)
         end if 
      else
         write (lu,1002) hdr,iparm(IIN2,jx)
      end if
!
!   -------------------
!     Linewidth tensor 
!   -------------------
      isym=iparm(IIWFLG,jx) 
      if (lu.ne.luttyo) write (lu,2002) symstr(isym), 'W'
      if (isym.eq.AXIAL) then
         write (lu,3003) hdr,(fparm(IWXX+i,jx),i=0,2,2)
      else if (isym.eq.SPHERICAL) then
         write (lu,3002) hdr,(fparm(IWXX+i,jx),i=0,2)
      else
        write (lu,3001) hdr,(fparm(IWXX+i,jx),i=0,2)
      end if
!
!     -------------------
!       Diffusion tensor 
!     -------------------
      isym=iparm(IIRFLG,jx)
      if (lu.ne.luttyo) write (lu,2002) symstr(isym), 'R'
      if (isym.eq.AXIAL) then
        write(lu,2006) hdr,(fparm(IDX+i,jx),i=0,2,2)
      else if (isym.eq.SPHERICAL) then
         write(lu,1006) hdr,(fparm(IDX+i,jx),i=0,2)
      else
         write(lu,1005) hdr,(fparm(IDX+i,jx),i=0,2)
      end if
!
!    ------------------------------------------
!      Gaussian inhomogeneous broadening
!    ------------------------------------------
      write (lu,1022) hdr, fparm(IGIB0,jx), fparm(IGIB2,jx)
!
!    -------------------------------------------
!      Non-Brownian rotational model parameters
!     ------------------------------------------
      write (lu,1007) hdr,iparm(IIPDF,jx)
      if (iparm(IIPDF,jx).eq.1) then
        do i=0,2
            write(lu,1008) hdr,nonbr(i+1),nonbr(i+1),iparm(IML+i,jx),
     #                     fparm(IPML+i,jx)
         end do
!
!   -- Anisotropic viscosity
!
      else if (iparm(IIPDF,jx).eq.2) then
        write (lu,1009) hdr,fparm(IDJF,jx),fparm(IDJFPRP,jx)
      end if
!
!   --- Discrete jump motion
!
      npotn=0
      if (iparm(IIPDF,jx).ne.2 .and. iparm(IIST,jx).ne.0 .and.
     #    npotn.gt.0) write (lu,1010) hdr,iparm(IIST,jx),fparm(IDJF,jx)
!
!    ----------------------
!      Orienting potential 
!    ----------------------
      npotn=0
      do i=0,4
        if (dabs(fparm(IC20+i,jx)).gt.RNDOFF) npotn=i+1
      end do
      if (npotn.gt.0) then
        write (lu,1011) hdr,(fparm(IC20+i,jx),i=0,4)
        if (iparm(INORT,jx).gt.0) write (lu,1012) hdr,iparm(INORT,jx)
      end if
!
!    --------------------
!    Heisenberg exchange
!    --------------------
      if (dabs(fparm(IOSS,jx)).gt.RNDOFF) 
     #    write (lu,1013) hdr,fparm(IOSS,jx)
!
!    -----------------
!      Diffusion tilt 
!    -----------------
      if (dabs(fparm(IBED,jx)).gt.RNDOFF) 
     #    write (lu,1014) hdr,fparm(IBED,jx)
!
!    -----------------------------
!     Basis set and CG parameters
!    -----------------------------
      write (lu,1016) hdr,(iparm(ILEMX+i,jx),i=0,NTRC-1)
      write (lu,1017) hdr,iparm(INSTEP,jx),fparm(ICGTOL,jx),
     #                fparm(ISHIFT,jx),fparm(ISHIFT+1,jx)
!
   10 continue
!
!
      return
!
!######################################################################
!   
 1000 format(a,'gxx,gyy,gzz = ',2(f9.6,','),f9.6)
 1001 format(a,'g1,g2,g3 = ',2(f9.6,','),f9.6)
 2001 format(a,'gprp,gpll = ',f9.6,',',f9.6)
 2002 format(a,' ',a)
 1002 format(a,'in2 =',i2)
 1003 format(a,'in2,Axx,Ayy,Azz = ',i2,',',2(f7.2,','),f7.2)
 1004 format(a,'in2,A1,A2,A3 = ',i2,',',2(f7.2,','),f7.2)
 2004 format(a,'in2,Aprp,Apll = ',i2,',',f7.2,',',f7.2)
 1005 format(a,'Rx,Ry,Rz = ',2(f8.4,','),f8.4)
 1006 format(a,'Rbar,N,Nxy = ',2(f8.4,','),f8.4)
 2006 format(a,'Rprp,Rpll = ',f8.4,',',f8.4)
 3001 format(a,'Wxx,Wyy,Wzz = ',2(f7.3,','),f7.3)
 3002 format(a,'W1,W2,W3 = ',2(f7.3,','),f7.3)
 3003 format(a,'Wprp,Wpll = ',f7.3,',',f7.3)
 1007 format(a,'ipdf = ',i1)
 1008 format(a,'m',a2,', pm',a2,' = ',i2,',',g10.3)
 1009 format(a,'djf,djfprp = ',f7.4,',',f7.4)
 1010 format(a,'ist,djf =',i3,',',f7.4)
 1011 format(a,'c20,c22,c40,c42,c44 = ',4(f7.3,','),f7.3)
 1012 format(a,'nort = ',i3)
 1013 format(a,'oss = ',f7.4)
 1014 format(a,'bed =',f7.2)
 1016 format(a,'lemx,lomx,kmn,kmx,mmn,mmx,ipnmx = ',6(i3,','),i3)
 1017 format(a,'nstep,cgtol,shiftr,shifti = ',i4,',',2(g10.3,','),g10.3)
 1018 format(/'sites =',i3)
 1019 format(/a,' ','*** Parameters for site',i2,' ***'/)
 1021 format(a,'lb(',i2,')=',f7.3)
 1022 format(a,'gib0,gib2 =',f7.3,',',f7.3)
 1023 format('series ',a,' = ',8(g10.3,' '))
 2023 format(/a,' ','*** No series defined ***'/)
 1024 format(a,'B0(',i2,')=',f10.3)
 1025 format(a,'psi(',i2,')=',f8.3)
 1026 format(a,'phase(',i2,')=',f7.2)
 3000 format('*** Unable to open file ',a,' for parameter output ***')
      end
