c Version 1.6  5/2/94
c**********************************************************************
c       
c                       ==============
c                       PROGRAM : LBLL
c                       ==============
c
c       This program will prompt the user to enter the parameters 
c       required as input to the EPR spectral calculation program.
c
c       Version 1.5 includes an orientation-dependent linewidth
c       Version 1.6 includes all three Euler angles in specifying
c                   diffusion tilt, three angles specifying magnetic
c                   tilt, and allows for a nuclear Zeeman interaction
c
c       written by DJS 15-NOV-87
c       Modified by DEB 15-MAR-92 and later
c
c       Includes:
c               stdio.inc
c               rndoff.inc
c               eprdat.inc
c               maxl.inc
c               stddim.inc
c               fnames.inc
c               indexl.inc
c
c       Uses:
c               ipar.f
c               cd2km.f
c               rddat.f
c               wrdat.f
c               lbasix.f
c               wrlbas.f
c
c**********************************************************************
c
      program lbll
c
      implicit none
c
      include 'stddim.inc'
      include 'maxl.inc'
      include 'stdio.inc'
      include 'eprdat.inc'
      include 'fnames.inc'
      include 'indexl.inc'
      include 'rndoff.inc'
c
      logical fprmpt,fexist,flagav,axiala,axialg,axialr,axialw
      integer i,ierr,ifile,j,k,l,listb,mxfld,nfiles,pndim,
     #        inlemx,inlomx,inkmn,inkmx,inmmn,inmmx,inpnmx,momd
      double precision d2km(2,5,5)
      character*80 line
c
      double precision ZERO,HALF,DSQ23
      parameter (ZERO=0.0D0,HALF=0.5D0)
c
      integer ipar,itoken
      double precision ftoken
      logical prmpt
      external ipar,prmpt,itoken,ftoken
c
      include 'version.inc'
c
c######################################################################
c
      DSQ23=dsqrt(2.0D0/3.0D0)
      flagav=.false.
c
      write (luttyo,1000)
      write (luttyo,1010) version,vdate
c
      call getids(nfiles,mxcalc)
      if (nfiles.eq.0) goto 9999
      ifile=0
c
    1 ifile=ifile+1
      if (ifile.gt.nfiles) goto 9999
      fileid=files(ifile)
      call setnam
c
c----------------------------------------------------------------------
c     check for existence and open parameter file
c----------------------------------------------------------------------
c
      inquire(file=prname,exist=fexist)
      if (fexist) then
         write (luttyo,1040) prname(:namlth)
         call rddat(prname,namlth,ierr)
         if (ierr.eq.-100) then
            write (luttyo,1041) prname(:namlth)
         elseif (ierr.lt.0) then
            write (luttyo,1042) prname(:namlth)
         end if
      else
         write (luttyo,1050) prname(:namlth)
      end if
c
c
c----------------------------------------------------------------------
c     initialize spherical tensors
c----------------------------------------------------------------------
c
        do j=1,5
           faa(j)=ZERO
           fgm(j)=ZERO
           do i=1,2
              fam(i,j)=ZERO
              fgd(i,j)=ZERO
              fad(i,j)=ZERO
           end do
        end do
c
c----------------------------------------------------------------------
c     get g-tensor
c----------------------------------------------------------------------
c
 11   write (luttyo,2000) gxx,gyy,gzz
      fprmpt=prmpt( line )
      if (fprmpt) then
         gxx=ftoken( line )
         gyy=ftoken( line )
         gzz=ftoken( line )
      end if

      if ( abs(gxx) .lt. RNDOFF  .or.
     #     abs(gyy) .lt. RNDOFF  .or.
     #     abs(gzz) .lt. RNDOFF) go to 11
c
c----------------------------------------------------------------------
c     get nuclear spin and A-tensor (in gauss)
c----------------------------------------------------------------------
c
      write (luttyo,2010) in2
      fprmpt=prmpt( line )
      if (fprmpt) in2=itoken( line )
      if (in2.le.0) then
         in2=0
      else
         write (luttyo,2020) axx,ayy,azz
         fprmpt=prmpt( line )
         if (fprmpt) then
            axx=ftoken( line )
            ayy=ftoken( line )
            azz=ftoken( line )
         end if
      end if
c
c----------------------------------------------------------------------
c     Get nuclear gyromagnetic ratio (in radians / sec*gauss)
c----------------------------------------------------------------------
      write (luttyo,2021) gamman
      fprmpt=prmpt( line )
      if (fprmpt) gamman=ftoken( line )
c
c----------------------------------------------------------------------
c     get linewidth tensor (peak-peak derivative Lorentzian, in gauss)
c----------------------------------------------------------------------
c
      write (luttyo,2280) wxx,wyy,wzz
      fprmpt=prmpt( line )
      if (fprmpt) then
         wxx=ftoken( line )
         wyy=ftoken( line )
         wzz=ftoken( line )
      end if
      if (wxx.lt.0) wxx=ZERO
      if (wyy.lt.0) wyy=ZERO
      if (wzz.lt.0) wzz=ZERO
c
c----------------------------------------------------------------------
c     get center field (in gauss)
c----------------------------------------------------------------------
c
 20   write (luttyo,2030) b0
      fprmpt=prmpt( line )
      if (fprmpt) b0=ftoken( line )
c
      if (b0 .lt. ZERO) b0=abs(b0)
      if (b0.lt.RNDOFF) go to 20
      if ( (in2 .le. 0) .or.
     #     (dabs(axx) .lt. RNDOFF .and. 
     #     dabs(ayy) .lt. RNDOFF .and.
     #     dabs(azz) .lt. RNDOFF)     ) then
         in2=0
         axx=ZERO
         ayy=ZERO
         azz=ZERO
      endif
c
c----------------------------------------------------------------------
c     calculated spherical components of F_g, F_A and F_w
c----------------------------------------------------------------------
c
c                                        *** g-matrix ***
      g0=(gxx+gyy+gzz)/3.0D0
      fgm(1)=HALF*(gxx-gyy)*(b0/g0)
      fgm(3)=DSQ23*(gzz-HALF*(gxx+gyy))*(b0/g0)
      fgm(5)=fgm(1)
      axialg=dabs(fgm(1)).lt.RNDOFF
c                                        *** Hyperfine tensor ***
      a0=(axx+ayy+azz)/3.0D0
      faa(1)=HALF*(axx-ayy)
      faa(3)=DSQ23*(azz-HALF*(axx+ayy))
      faa(5)=faa(1)
      axiala=dabs(faa(1)).lt.RNDOFF
c                                        *** Linewidth tensor ***
      w0=(wxx+wyy+wzz)/3.0D0
      fwm(1)=HALF*(wxx-wyy)
      fwm(3)=DSQ23*(wzz-HALF*(wxx+wyy))
      fwm(5)=fwm(1)
      axialw=dabs(faa(1)).lt.RNDOFF
c
c----------------------------------------------------------------------
c     get diffusion model index (0 to 2)
c                         0: Brownian motion
c                         1: Nonbrownian motion
c                         2: Anisotropic viscosity
c----------------------------------------------------------------------
c
      write (luttyo,2040) ipdf
      fprmpt=prmpt( line )
      if (fprmpt) ipdf=itoken( line )
      if (ipdf.gt.2) ipdf=2
      if (ipdf.lt.0) ipdf=0
c
c----------------------------------------------------------------------
c     get diffusion tensor (per second)
c----------------------------------------------------------------------
c
      write (luttyo,2050) dx,dy,dz
      fprmpt=prmpt( line )
      if (fprmpt) then
         dx=ftoken( line )
         dy=ftoken( line )
         dz=ftoken( line )
      end if
      axialr=dabs(dx-dy).le.RNDOFF
c
c----------------------------------------------------------------------
c     get more motional parameters if needed
c----------------------------------------------------------------------
c
      if (ipdf.eq.0) go to 200
c
      if (ipdf.eq.2) then
         write (luttyo,2060) tl,tkxy,tkzz
      else
         write (luttyo,2070)  tl,tkxy,tkzz
      end if
      fprmpt=prmpt( line )
      if (fprmpt) then
         tl=ftoken( line )
         tkxy=ftoken( line )
         tkzz=ftoken( line )
      endif
      if ( (ipdf.eq.2) .and.
     #     (abs(tl).lt.RNDOFF) .and.
     #     (abs(tkxy).lt.RNDOFF) .and.
     #     (abs(tkzz).lt.RNDOFF)       ) then
         flagav=.true.
         go to 200
      else
         flagav=.false.
         djfprp=ZERO 
      end if
c
 210  write (luttyo,2080) pl,pkxy,pkzz
      fprmpt=prmpt( line )
c
      if (fprmpt) then
         pl=ftoken( line ) 
         pkxy=ftoken( line )
         pkzz=ftoken( line )
      end if
c
 200  continue
c
c----------------------------------------------------------------------
c     get number of sites and jump rate for site jump model
c----------------------------------------------------------------------
c
      write (luttyo,2090) ist,djf
      fprmpt=prmpt( line )
      if (fprmpt) then
         ist=itoken( line )
         djf=ftoken( line )
      end if
c
c----------------------------------------------------------------------
c     get perpendicular rate for anisotropic viscosity term
c----------------------------------------------------------------------
c
      if (flagav.and.(ist.eq.0)) then
         write (luttyo,2073) djfprp
         fprmpt=prmpt( line )
         if (fprmpt) djfprp=ftoken( line )
      else
         djfprp=ZERO
      end if            
c
c----------------------------------------------------------------------
c     get Heisenberg spin exchange rate (should be positive)
c----------------------------------------------------------------------
c
      write(luttyo,2100) oss
      fprmpt=prmpt( line )
      if (fprmpt) oss=abs(ftoken( line ))
      if (oss.lt.RNDOFF) oss=ZERO
c
c----------------------------------------------------------------------
c     get potential parameters
c----------------------------------------------------------------------
c
      write (luttyo,2110) ipt
      fprmpt=prmpt( line )
      if (fprmpt) ipt=itoken( line )
      if (ipt.lt.0) ipt=0
      if (ipt.gt.6) ipt=6
      if (ipt.gt.0) then
         write(luttyo,2120)
         do j=1,ipt
            write(luttyo,2130) ipt,(ixp(j,i),i=1,2),cxp(j)
            fprmpt=prmpt( line )
            if (fprmpt) then
               ixp(j,1)=itoken( line )
               ixp(j,2)=itoken( line )
               cxp(j)=ftoken( line )
            end if
         end do
c
         do j=1,ipt
            i=ixp(j,1)
            if(i.gt.4) i=4
            if(i.lt.2) i=2
            if(i.eq.3) i=2
            ixp(j,1)=i
            i=ixp(j,2)
            if(i.gt.ixp(j,1)) i=ixp(j,1)
            if(i.lt.0) i=0
            if(ipar(i).eq.-1) i=i-1
            ixp(j,2)=i
         end do
      end if
c
      do j=1,5
         do i=1,5
            cpot(i,j)=ZERO
         end do
      end do
c
      lptmx=0
      kptmx=0
c
      do i=1,ipt
         l=ixp(i,1)
         k=ixp(i,2)
         if(l.gt.lptmx) lptmx=l
         if(k.gt.kptmx) kptmx=k
         l=1+l/2
         k=1+k/2
         cpot(l,k)=cxp(i)
      end do
c
      if (lptmx.ge.2) then
         lband=lptmx*2
      else
         lband=2
      end if
c
      kband=kptmx*2
      if (.not.axialr) kband=kband+2
c
c----------------------------------------------------------------------
c Check specified motional model parameters for consistency with
c the presence of a potential
c----------------------------------------------------------------------
c
c     *** Non-Brownian motion: must have at least one nonzero residence
c         time, and no potential specified
c
      if (ipdf .eq. 1) then
         tl=abs(tl)
         tkxy=abs(tkxy)
         tkzz=abs(tkzz)
c
         if ((tl .lt. RNDOFF) .and. (tkxy .lt. RNDOFF) .and.
     #       (tkzz .lt. RNDOFF) )  then
            write (luttyo,5004)
            ipdf=0
         endif
c
         if (ipt .gt. 0) then
            write (luttyo,5005)
            ipdf=0
         endif
      endif
c
c     *** Anisotropic viscosity model: must have potential and
c         no discrete jump motion specified
c
      if (ipdf .eq. 2) then
         if (ipt .eq. 0) then
            write (luttyo,5006)
            ipdf=0
            flagav=.false.
         endif
c
        if (ist .gt. 0) then
           write (luttyo,5007)
           ipdf=0
           flagav=.false.
        endif
      endif
c
c     *** Brownian motion model: set residence times to ZERO
c
      if (ipdf .eq. 0) then
         tl  =ZERO
         tkxy=ZERO
         tkzz=ZERO
      endif
c
c--------------------------------------------------------------------------
c     get MOMD flag or director tilt if there is a non-vanishing potential
c--------------------------------------------------------------------------
c
      if (ipt.gt.0) then
         write (luttyo,2135) nort
         fprmpt=prmpt( line )
         if (fprmpt) nort=itoken(line)
         if (nort.gt.MXORT) nort=MXORT
         if (nort.lt.2) then
            nort=1
	    momd=0
         else
            momd=1
         endif
c
         if (momd.eq.0) then
            nort=1
            write(luttyo,2140) psi
            fprmpt=prmpt( line )
            if (fprmpt) psi=ftoken( line )
         end if
c
      else
         psi=ZERO
         momd=0
         nort=1
      end if
c
      if( ( (dabs(psi).lt.RNDOFF).or.(abs(psi-180.0d0).lt.RNDOFF) )
     #     .and. momd.eq.0) then
         ipsi0=0
      else
         ipsi0=1
      end if
c
c----------------------------------------------------------------------
c     get MOMD parameters
c----------------------------------------------------------------------
c
c      *** Orientation-dependent Gaussian inhomogeneous broadening
c
      if (momd.ne.0) then
         write (luttyo,2240) gib0, gib2
         fprmpt=prmpt( line )
         if (fprmpt) then
            gib0=ftoken( line )
            gib2=ftoken( line )
         end if
c
c      *** Derivative-mode flag
c
         write (luttyo,2260) ideriv
         fprmpt=prmpt( line )
         if (fprmpt) ideriv=itoken( line )
         if (ideriv .ne. 0) ideriv=1
      end if
c
c----------------------------------------------------------------------
c     Get diffusion tilt angles
c----------------------------------------------------------------------
c
      write(luttyo,2160) ald,bed,gad
      fprmpt=prmpt( line )
      if (fprmpt) then
         ald=ftoken( line )
         bed=ftoken( line )
         gad=ftoken( line )
      end if
c
c     *** diff. tilt angle alpha not needed if diff. tensor is axial
c
      if (axialr .and. dabs(ald).gt.RNDOFF) then
        ald=ZERO
        write(luttyo,5015)
      end if
c
      bed=dabs( bed )
      itd=0
      if ( dabs(ald).gt.RNDOFF .or.
     #     dabs(gad).gt.RNDOFF .or. 
     #     bed.gt.RNDOFF) itd=1
c
c----------------------------------------------------------------------
c     Get magnetic tilt angles
c----------------------------------------------------------------------
c
      write(luttyo,2162) alm,bem,gam
      fprmpt=prmpt( line )
      if (fprmpt) then
        alm=ftoken( line )
        bem=ftoken( line )
        gam=ftoken( line )
      end if
c
c     *** mag. tilt angle alpha not needed if g-matrix is axial
c
      if (axialg .and. axialw .and. dabs(alm).gt.RNDOFF) then
        alm=ZERO
        write(luttyo,5014)
      end if
c
      if (dabs(alm).gt.RNDOFF .or. bem.gt.RNDOFF 
     #    .or. dabs(gam).gt.RNDOFF) itm=1
c
c----------------------------------------------------------------------
c     Set A tensor in magnetic (g-tensor) frame
c----------------------------------------------------------------------
      if (itm .eq. 0) then
c                           *** no magnetic tilt: copy A tensor directly
         do j=1,5
            fam(1,j)=faa(j)
         end do
      else
c                          *** transform A tensor into g axis system
c
         call cd2km(d2km,alm,bem,gam)
         do i=1,5
            do j=1,5
               fam(1,i)=fam(1,i)+d2km(1,i,j)*faa(j)
               fam(2,i)=fam(2,i)+d2km(2,i,j)*faa(j)
            end do
         end do
      end if
c
c----------------------------------------------------------------------
c     Set F_g and F_A tensors in the diffusion frame
c----------------------------------------------------------------------
      if (itd.eq.0) then
c                           *** no diffusion tilt: copy tensors directly
         do j=1,5
            fgd(1,j)=fgm(j)
            fwd(1,j)=fwm(j)
            fad(1,j)=fam(1,j)
            fad(2,j)=fam(2,j)
         end do
c
      else
c                    *** Transform A, g, W tensors into the diffusion frame
c
         call cd2km(d2km,ald,bed,gad)
         do i=1,5
            do j=1,5
               fgd(1,i)=fgd(1,i)+d2km(1,i,j)*fgm(j)
               fgd(2,i)=fgd(2,i)+d2km(2,i,j)*fgm(j)
               fwd(1,i)=fwd(1,i)+d2km(1,i,j)*fwm(j)
               fwd(2,i)=fwd(2,i)+d2km(2,i,j)*fwm(j)
               fad(1,i)=fad(1,i)+d2km(1,i,j)*fam(1,j)
     #                          -d2km(2,i,j)*fam(2,j)
               fad(2,i)=fad(2,i)+d2km(1,i,j)*fam(2,j)
     #                          +d2km(2,i,j)*fam(1,j)
            end do
         end do
c
      end if
c
c----------------------------------------------------------------------
c     Get basis set truncation parameters
c----------------------------------------------------------------------
c
 272  write(luttyo,2170) lemx,lomx,kmx,mmx,ipnmx
      fprmpt=prmpt( line )
      if (fprmpt) then
         lemx=itoken( line )
         lomx=itoken( line )
         kmx=itoken( line )
         mmx=itoken( line )
         ipnmx=itoken( line )
      end if
c
c----------------------------------------------------------------------
c Use antisymmetric K combinations (jmmn=-1) if the magnetic 
c tensors in the diffusion frame are complex-valued. This can result
c from nonzero alpha or gamma magnetic or diffusion tilt angles.
c----------------------------------------------------------------------
      jkmn=1
      do i=1,5
         if (     abs(fgd(2,i)).ge.RNDOFF
     #       .or. abs(fad(2,i)).ge.RNDOFF
     #       .or. abs(fwd(2,i)).ge.RNDOFF ) jkmn=-1
      end do
c
c     --- Get truncation value for antisymmetric K-space if needed
c
      if (jkmn.eq.-1) then
         kmn=-kmx
         write(luttyo,2171) kmn
         fprmpt=prmpt( line )
         if (fprmpt) kmn=itoken( line )
      else
        kmn=0
      endif
c
c----------------------------------------------------------------------
c Use antisymmetric M/pI combinations with a nonzero nuclear Zeeman 
c----------------------------------------------------------------------
      jmmn=1
      if (dabs(gamman).ge.RNDOFF) jmmn=-1
      if (jmmn.eq.1) mmn=0
c
c     Get truncation value for antisymmetric M-space if needed
c
      if (jmmn.eq.-1) then
         mmn=-mmx
         write(luttyo,2172) mmn
         fprmpt=prmpt( line )
         if (fprmpt) mmn=itoken( line )
      else
         mmn=0
      endif
c
c----------------------------------------------------------------------
c     Check basis set parameters
c----------------------------------------------------------------------
      inlemx=lemx
      inlomx=lomx
      inkmn=kmn
      inkmx =kmx
      inmmn =mmn
      inmmx =mmx
      inpnmx=ipnmx
c
      if(lemx.lt.0) lemx=0
      if(lomx.lt.0) lomx=0
      if(ipar(lemx).ne.1) lemx=lemx-1
      if(lomx.gt.lemx) lomx=lemx-1
c
      if(kmx.lt.0) kmx=0
      if(kmn.gt.0) kmn=0
      if(kmx.gt.lemx) kmx=lemx
      if(ipar(kmx).ne.1) kmx=kmx-1
      if(ipar(kmn).ne.1) kmn=kmn+1
      if(-kmn.gt.kmx) kmn=-kmx
c
      if(mmx.gt.lemx) mmx=lemx
      if(mmn.gt.0) mmn=0
      if(-mmn.gt.mmx) mmn=-mmx
c
      if((lemx.gt.MXLVAL)) then
         write (luttyo,5008) MXLVAL
         lemx=MXLVAL
      endif
c
      if(mmx.lt.0) mmx=0
      if(ipnmx.lt.0) ipnmx=0
      if(ipnmx.gt.in2) ipnmx=in2
      if((ipsi0.eq.0).and.nort.le.1.and.(mmx.gt.ipnmx)) mmx=ipnmx
c
c----------------------------------------------------------------------
c  Determine basis set using rules in Meirovitch,Igner,Igner,Moro & Freed
c  (J. Phys. Chem. 77 (1982) p. 3935)
c----------------------------------------------------------------------
c
c----------------------------------------------------------------------
c Use only even K values (kdelta=2) if there is no magnetic or
c diffusion tilt.
c----------------------------------------------------------------------
      if (itm.eq.0 .and.itd.eq.0) then
        kdelta=2
      else
        kdelta=1
      end if
c
c----------------------------------------------------------------------
c Use only even L values (ldelta=2) and no K values (kmx=0) in the case
c   (1) there is no magnetic or diffusion tilt (kdelta.eq.2)
c   (2) magnetic tensors are axial
c   (3) potential is axial (no K terms in potential)
c----------------------------------------------------------------------
      if (axiala .and. axialg .and. axialw .and. kdelta.eq.2 .and. 
     #   kptmx.eq.0) then
         ldelta=2
         kmx=0
      else
         ldelta=1
      end if
c
c----------------------------------------------------------------------
c     Report any changes in basis truncation parameters
c----------------------------------------------------------------------
      if (inlemx .ne. lemx .or.
     #    inlomx .ne. lomx .or.
     #    inkmn  .ne. kmn  .or.
     #    inkmx  .ne. kmx  .or.
     #    inmmn  .ne. mmn  .or.
     #    inmmx  .ne. mmx  .or.
     #    inpnmx .ne. ipnmx)
     # write (luttyo,5009) lemx,lomx,kmn,kmx,mmn,mmx,ipnmx
c
      call lbasix( ixname, namlth, 1 )
      write(luttyo,5010) ndim
c
      if (ndim.gt.MXDIM) then
         write (luttyo,5011) MXDIM
         goto 272
      end if
c
      inquire(file=ixname(:namlth),exist=fexist)
      if (fexist) then
         open (unit=ludisk,file=ixname(:namlth),status='old',
     #        access='sequential',form='unformatted')
            read (ludisk) pndim
            close (unit=ludisk)
            write (luttyo,5012) ixname(:namlth),pndim
      end if	
c

c----------------------------------------------------------------------
c     get number of Lanczos/CG steps to execute and
c     the maximum allowed error for CG
c----------------------------------------------------------------------
c
      if (nstep.gt.ndim) nstep=ndim
      write(luttyo,2180) nstep
      fprmpt=prmpt( line )
      if (fprmpt) nstep=itoken( line )
      if (nstep.gt.mxstep) nstep=mxstep
c
      write (luttyo,2190) itype
      fprmpt=prmpt( line )
      if (fprmpt) itype=itoken( line )
      if (itype.lt.0) itype=0
      if (momd.ne.0 .and. itype.eq.0) itype=1
c
c --- Specify CG tolerance for CG (or MOMD) calculations
c
      if (itype.ge.1 .or. momd.ne.0) then
         write (luttyo,2200) cgtol
         fprmpt=prmpt( line )
         if (fprmpt) cgtol=ftoken( line )
         if (cgtol.lt.10.0D0*RNDOFF) cgtol=10.0D0*RNDOFF
c
         write (luttyo,2210) shiftr,shifti
         fprmpt=prmpt( line )
         if (fprmpt) then
            shiftr=ftoken( line )
            shifti=ftoken( line )
         end if
         if (shiftr.lt.10.0D0*RNDOFF) shiftr=10.0D0*RNDOFF
      else
         shiftr=ZERO 
         shifti=ZERO 
      end if
c
c   Choose a field range and number of fields for either 
c   field-swept or MOMD calculation
c
      if (itype.ge.2 .or. momd.ne.0) then
         write (luttyo,2220) b0,fieldi,fieldf
         fprmpt=prmpt( line )
         if (fprmpt) then
            fieldi=ftoken( line )
            fieldf=ftoken( line )
        end if
c
        mxfld=mxcalc
        if (momd.ne.0) mxfld=mxpt
        write (luttyo,2230) mxfld+1,nfield
        fprmpt=prmpt( line )
        if (fprmpt) nfield=itoken( line )
        if (nfield.gt.mxfld) nfield=mxfld
        if (nfield.lt.2) nfield=2
      end if
c
c   Set pruning tolerance for EPRBL calculation
c
      if (itype.ge.2) then
         write (luttyo,2270) btol
         fprmpt=prmpt( line )
         if (fprmpt) btol=ftoken( line )
         if (btol.lt.ZERO) btol=ZERO
         if (btol.gt.HALF) btol=HALF
      end if
c         
c----------------------------------------------------------------------
c     print out new values of parameters
c----------------------------------------------------------------------
c
      call wrparm(prname,namlth,luttyo) 
c
      listb=0
      write(luttyo,4000) fmname(:namlth)
      read (luttyi,4001) line
      if (line(1:1).eq.'y' .or. line(1:1).eq.'Y') listb=1
c     
c
c----------------------------------------------------------------------
c     write parameters to formatted disk file
c----------------------------------------------------------------------
      open (unit=lulog,file=fmname,status='unknown',
     #     access='sequential',form='formatted')
      rewind (unit=lulog)
c
      call wrparm(prname,namlth,lulog)
      if (listb.eq.1) call wrlbas(prname,namlth,lulog)
c
      close (unit=lulog)
c
c----------------------------------------------------------------------
c     write parameters to binary disk file 
c----------------------------------------------------------------------
c
      call wrdat(prname,namlth)
      go to 1
c
 9999 write(luttyo,1000)
      stop
c
c======================================================================
c               format statements
c======================================================================
c
 1000 format(//,2x,70('#'),//)
 1010 format(25x,'program LBLL'/22x,'Version ',a,1x,a/
     #       22x,'----------------')
 1040 format(2x,'Data file ',a,' already exists on the disk')
 1041 format(2x,'*** ''',a,''' is not an EPRLL Version 1.5 .PAR file',
     #' ***')
 1042 format(2x,'*** Error reading file ''',a,''' ***')
 1050 format(2x,'Data file ',a,' not found')
c
 2000 format(2x,'g-tensor [gxx,gyy,gzz] : ',3(f10.7,2x))
 2010 format(2x,'twice the nuclear spin [in2] : ',i2)
 2020 format(2x,'A-tensor [axx,ayy,azz] (gauss) : ',3(f9.3,2x))
 2021 format(2x,'Nuclear gyromagnetic ratio [gamman] ',
     #       'radian/(sec*gauss) : ',g13.6)
 2030 format(2x,'static field [B0] (gauss) : ',f10.3)
c
 2040 format(2x,'diffusion parameter [ipdf] : ',i1,/,10x,
     #     'ipdf=0 --> brownian diffusion',/,10x,
     #     'ipdf=1 --> free or jump diffusion',/,10x,
     #     'ipdf=2 --> anisotropic viscosity (if ist=0)')
 2050 format(2x,'diffusion tensor [dx,dy,dz] (1/sec) =',3(g11.4,2x))
 2060 format(2x,'residence times [tl,tkxy,tkzz] (sec)',/,
     #     5x,'(brownian model-->tl,tkxy,tkzz=0) ',/,2x,3(1x,g11.4))
 2070 format(2x,'residence times [tl,tkxy,tkzz]  (sec)',/,2x,
     #     8x,3(1x,g11.4))
 2073 format(2x,'perpendicular diffusion coefficient in the director',
     #     /,5x,'frame [djfprp] (1/sec) : ',g11.4)
 2080 format(2x,'non-Brownian exponents [pl,pkxy,pkzz] : ',/,
     #     10x,'0.5 --> free',/,10x,'1.0 --> jump',/,
     #     2x,3(f6.2,2x))
 2090 format(2x,'discrete jumps motion [ist,djf] : ',
     #     i2,2x,g13.6,/,4x,
     #     '(if ipdf=2 then ist=0 and djf is the parallel',/,4x,
     #     'diffusion coefficient in the director frame)')
c
 2100 format(2x,'Heisenberg spin exchange frequency [oss] = ',g11.4)
c
 2110 format(2x,'number of terms in the potential [ipt] = ',i1)
 2120 format(2x,'coefficients of the potential : ')
 2130 format(10x,'ipt = ',i1,2x,'[l,k,coef.] = ',2(i1,2x),f8.4)
 2135 format(2x,'number of orientations (>1 for MOMD) [nort] :',i5)
 2140 format(2x,'angle between B0 and local director [psi] ',
     #     '(degrees) : ',g13.6)
 2160 format(2x,'diffusion tilt [ald,bed,gad] (degrees): ',
     # 2(f6.1,','),f6.1)
 2162 format(2x,'magnetic tilt [alm,bem,gam] (degrees): ',
     # 2(f6.1,','),f6.1)
 2170 format(2x,'truncation values  [lemx,lomx,kmx,'
     #     ,'mmx,ipnmx] : ',4(i3,','),i3)
 2171 format(2x,'truncation for antisymmetric K-space [kmn] : ',i3 )
 2172 format(2x,'truncation for antisymmetric M-space [mmn] : ',i3 )
 2180 format(2x,'number of Lanczos/CG steps [nstep] : ',i5)
 2190 format(2x,'calculation type (0=L,1=CG,2=FS) [itype] : ',i1)
 2200 format(2x,'error tolerance for residual [cgtol] : ',g11.4)
 2210 format(2x,'origin shift [shiftr,shifti] :',2(2x,g11.4))
 2220 format(2x,'initial and final fields relative to B0=', f10.3,
     #     /,2x,'[fieldi,fieldf] (gauss) : ',f10.3,2x,f10.3)
 2230 format(2x,'number of field positions [nfield<',i3,'] : ',i3)
 2240 format(2x,'Gaussian inhomog. broadening [gib0,gib2] (gauss) :',
     #     2x,f7.3,2x,f7.3)
 2260 format(2x,'derivative mode [ideriv] : ',i2)
 2270 format(2x,'pruning tolerance [btol] : ',g9.3)
 2280 format(2x,'Linewidth tensor [wxx,wyy,wzz] (gauss) : ',
     #     3(f6.2,2x))
c
 4000 format(2x,'Write basis set indices to ',a,'? [N]: ')
 4001 format(a)
 5004 format(' Times for jump/free model are all zero:',
     #       ' ipdf=0 assumed')
 5005 format(' Nonzero potential with jump/free model:',
     #       ' ipdf=0 assumed') 
 5006 format(' Zero potential with anisotropic viscosity:',
     #       ' ipdf=0 assumed')
 5007 format(' Discrete jumps with anisotropic viscosity:',
     #       ' ipdf=0 assumed')
 5008 format(' lemx must be ',i4,' or less')
 5009 format(' ** Basis set parameters adjusted to [lemx,lomx,',
     #       'kmn,kmx,mmn,mmx,ipnmx]  = ',/,10x,6(i3,','), i3,' **')
 5010 format(/,2x,'*** Dimension of basis = ',i6,' ***',/)
 5011 format(/,18x,'*** WARNING: basis dimension exceeds ',i7,' ***'/
     #     2x,'Please respecify basis truncation parameters'/)
 5012 format(/' File ''',a,''' exists: pruned dimension =',i6/)
 5014 format(' Angle alm not needed with axial g-matrix')
 5015 format(' Angle ald not needed with axial diffusion tensor')
c
c======================================================================
c
      end
