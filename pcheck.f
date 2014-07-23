c  VERSION 1.0  (NLSPMC version)   2/5/99
c**********************************************************************
c       
c                        =================
c                        subroutine PCHECK
c                        =================
c
c       This procedure replaces the function of program LBLF & TDIP for
c       performing initial computations based on the input parameters
c       for the EPRP 2D calculation programs.  It assumes the floating
c       point and integer variables in the fparm and iparm arrays have 
c       been copied into common block /eprprm/ and checks to make sure
c       they are valid.
c
c       Returns ierr=0 if parameters are valid;  a negative number
c       identifying a fatal error, and a positive number issuing warning
c       messages. (it is reset to zero before returning to the calling
c       routine.)
c
c       lumsg specifies a logical unit for output of any error messages.
c       None are produced if lumsg=0.
c
c       Fatal parameter errors
c          1) zero gxx, gyy, or gzz
c          2) zero B0
c      
c       Includes:
c               nlsdim.inc
c               nlsnam.inc
c               maxl.inc
c               eprprm.inc
c               rndoff.inc
c       Uses:
c               ipar.f
c               cd2km.f
c               tocart.f (coded in tensym.f)
c               fbasis.f
c
c**********************************************************************
c
      subroutine pcheck(lumsg,ierr)
c
      implicit none
c
      include 'limits.inc'
      include 'names.inc'
      include 'maxl.inc'
      include 'simparm.inc'
      include 'rndoff.inc'
c
      integer lumsg,ierr
c
      integer i,j,itmp,inlemx,inlomx,inkmx,inmmx,inpnmx
      double precision d2km(2,5,5),tmp,gmax,gmin,hmax
c
      double precision zero,half,dsq23,one,ten
      parameter (zero=0.0D0,half=0.5D0,one=1.0D0,ten=1.0D1)
      parameter (dsq23=0.816496580927726d0)
c
      logical axiala,axialg 
c
      integer ipar
      external ipar
c
c######################################################################
c
      ierr=0
c
c----------------------------------------------------------------------
c     initialize spherical tensors
c----------------------------------------------------------------------
        do 10 j=1,5
          faa(j)=zero
          fgm(j)=zero
          do 9 i=1,2
            fam(i,j)=zero
            fgd(i,j)=zero
            fad(i,j)=zero
    9       continue
 10       continue
c
c----------------------------------------------------------------------
c     Check for zero field
c----------------------------------------------------------------------
      if (b0 .lt. zero) then
         if (lumsg.ne.0) write (lumsg,1002)
         ierr=2
         b0=abs(b0)
      elseif (b0 .le. rndoff) then
         if (lumsg.ne.0) write (lumsg,2001)
         ierr=-1
         return
      endif
c
c----------------------------------------------------------------------
c     Convert g, A, and R tensors to Cartesian form
c----------------------------------------------------------------------
      call tocart( gxx, igflg )
      call tocart( axx, iaflg )
      call tocart( dx,  irflg )
c
c  *** NOTE: in NLS version, diffusion rates are specified as
c            a power of 10
c
c     set criterion for weighting factor to store eigenvector
c
      tmp=7.99
      do 30 i=1,6
        if ((dx.gt.tmp).and.(dy.gt.tmp).and.(dz.gt.tmp)) then
          go to 40
        else
          tmp=tmp-1.0d0
        end if
 30   continue
c
c  Wtol eliminated.  RC 9/5/96
c
 40   wtol=1.0D-9
c
      dx=ten**dx
      dy=ten**dy
      dz=ten**dz
c
c----------------------------------------------------------------------
c     Check for zero values in g-tensor and zero A tensor
c----------------------------------------------------------------------
      if ( abs(gxx) .lt. rndoff  .or.
     #     abs(gyy) .lt. rndoff  .or.
     #     abs(gzz) .lt. rndoff)   then
        if (lumsg.ne.0) write (lumsg,2002)
        ierr=-2
        return
      endif

      if ( (in2 .le. 0) .or.
     #    (dabs(axx) .lt. rndoff .and. 
     #     dabs(ayy) .lt. rndoff .and.
     #     dabs(azz) .lt. rndoff)     ) then
        in2=0
        axx=zero
        ayy=zero
        azz=zero
        if (lumsg.ne.0) write (lumsg,1001)
        ierr=1
      endif
c
      gamman=zero
c
c----------------------------------------------------------------------
c     Calculate spherical components of tensors in magnetic frames
c----------------------------------------------------------------------
c                                        *** Electronic Zeeman ***
      g0=(gxx+gyy+gzz)/3.0D0
      fgm(1)=half*(gxx-gyy)*(b0/g0)
      fgm(3)=dsq23*(gzz-half*(gxx+gyy))*(b0/g0)
      fgm(5)=fgm(1)
      axialg = dabs(fgm(1)).lt.rndoff
c                                        *** Hyperfine ***
      a0=(axx+ayy+azz)/3.0D0
      faa(1)=half*(axx-ayy)

      faa(3)=dsq23*(azz-half*(axx+ayy))
      faa(5)=faa(1)
      axiala = dabs(faa(1)).lt.rndoff
c                                        *** Nuclear Zeeman ***
      zeen=zero
c
      gmin=dmin1( gxx, gyy, gzz )
      gmax=dmax1( gxx, gyy, gzz )
      hmax=dmax1( dabs(axx), dabs(ayy), dabs(azz) )
c
c----------------------------------------------------------------------
c Issue warning if high-field approximation has been violated
c (this is a very rough criterion for the high-field approx.!)
c----------------------------------------------------------------------
      if (b0 .lt. 10.0D0*dmax1( hmax, (gmax-gmin)*b0/g0) ) then
         if (lumsg.ne.0) write (lumsg,1003)
         ierr=3
      endif
c----------------------------------------------------------------------
c  Set ipt and cpot array according to the potential coefficients given
c----------------------------------------------------------------------
      ipt=0
      do 240 j=1,5
        do 245 i=1,5
          cpot(i,j)=zero
 245    continue
 240  continue

      if (dabs(c20) .gt. rndoff) then
        ipt=ipt+1
        cpot(2,1)=c20
        lptmx=2
      endif
      if (dabs(c22) .gt. rndoff) then
        ipt=ipt+1
        cpot(2,2)=c22
        lptmx=2
        kptmx=2
      endif
      if (dabs(c40) .gt. rndoff) then
        ipt=ipt+1
        cpot(3,1)=c40
        lptmx=4
      endif
      if (dabs(c42) .gt. rndoff) then
        ipt=ipt+1
	cpot(3,2)=c42
        lptmx=4
        kptmx=2
      endif
      if (dabs(c44) .gt. rndoff) then
        ipt=ipt+1
        cpot(3,3)=c44
        lptmx=4
        kptmx=4
      endif
c
      if (lptmx.ge.2) then
        lband=lptmx*2
      else
        lband=2
      end if
c
      kband=kptmx*2
      if (dabs(dx-dy).gt.rndoff) kband=kband+2
c
c----------------------------------------------------------------------
c Check for consistency of specified motion parameters with given model
c----------------------------------------------------------------------
      if (ipdf.gt.2) ipdf=2
      if (ipdf.lt.0) ipdf=0
      if (ist.lt.0) ist=0
c 
c     Set exponents according to non-Brownian model
c
      pl=one
      pkxy=one
      pkzz=one
      if (ml.eq.1)  pl=half
      if (mxy.eq.1) pkxy=half
      if (mzz.eq.1) pkzz=half
c
c     *** Non-Brownian motion: must have at least one nonzero residence
c         time, and no potential specified
c         
c NOTE: in NLS version, R*residence time products are specified as 
c       powers of ten
c
      if (ipdf .ne. 0) then
         pml = ten**pml
         pmxy = ten**pmxy
         pmzz = ten**pmzz
c 
         if ((ml.eq.0) .and.(mxy.eq.0) .and. (mzz.eq.0) ) then
            if (lumsg.ne.0) write (lumsg,1004)
            ierr=4
            ipdf=0
         end if
c
c  *** NOTE: in NLS version, djf, djfprp specified as powers of ten
c
         djf=ten**djf
         djfprp=ten**djfprp
c
         if (ipt .gt. 0) then
            if (lumsg.ne.0) write (lumsg,1005)
            ierr=5
            ipdf=0
         end if
      end if
c
c     *** Anisotropic viscosity model: must have potential and
c         no discrete jump motion specified
c
      if (ipdf .eq. 2) then
         if (ipt .eq. 0) then
            if (lumsg.ne.0) write (lumsg,1006)
            ierr=6
            ipdf=0
         end if
c
         if (ist .gt. 0) then
            if (lumsg.ne.0) write (lumsg,1007)
            ierr=7
            ipdf=0
         end if
      end if
c
c     *** Brownian motion model: set residence times to zero
c
      if (ipdf .eq. 1) then
         if (ml.eq.0)  pml =zero
         if (mxy.eq.0) pmxy=zero
         if (mzz.eq.0) pmzz=zero
      end if
c
c     *** MOMD calculation: must have potential
c
      if (nort.gt.1) then
        if (nort.gt.MXORT) then
          if (lumsg.ne.0) write (lumsg,*)'nort too big, reset to ',
     #			MXORT
          nort=MXORT
        end if
        if (ipt.lt.1) then
          if (lumsg.ne.0) write (lumsg,1015)
          ierr=15
        end if
      end if
c
c----------------------------------------------------------------------
c     Set director tilt flag
c----------------------------------------------------------------------
      ipsi0=0
      if (nort.gt.1) ipsi0=1
      if((abs(psi).gt.rndoff).and.(abs(psi-180.0d0).gt.rndoff)) ipsi0=1
c
c----------------------------------------------------------------------
c     Set diffusion tilt flag
c----------------------------------------------------------------------
      bed=abs( bed )
      itd=0
      if ( dabs(ald).gt.rndoff .or. dabs(gad).gt.rndoff
     #    .or. bed.gt.rndoff) itd=1
c
c----------------------------------------------------------------------
c     Get magnetic tilt angles
c----------------------------------------------------------------------
      bem = dabs( bem )
      if (axialg .and. dabs(alm).gt.rndoff) then
         if (lumsg.ne.0) write (lumsg,1013)
         ierr=13
         alm=zero
      end if
c
c     *** Only beta tilt angles allowed if all tensors are axial
c
      if (axialg.and.axiala.and.(dabs(alm).gt.rndoff.or.
     #    dabs(gam).gt.rndoff.or.dabs(gad).gt.rndoff) ) then
         if (lumsg.ne.0) write (lumsg,1014)
         ierr=14
         alm=zero
         gam=zero
         gad=zero
      end if
c
      itm=0
      if (dabs(alm).gt.rndoff .or. bem.gt.rndoff 
     #    .or. dabs(gam).gt.rndoff)     itm=1
c
c----------------------------------------------------------------------
c     Set A tensor in magnetic (g-tensor) frame
c----------------------------------------------------------------------
      if (itm .eq. 0) then
c                           *** no tilt: copy tensor directly
         do 252 j=1,5
            fam(1,j)=faa(j)
  252       continue
      else
c                          *** transform A tensor into g axis system
c 
         call cd2km(d2km,alm,bem,gam)
         do 255 i=1,5
            do 254 j=1,5
               fam(1,i)=fam(1,i)+d2km(1,i,j)*faa(j)
               fam(2,i)=fam(2,i)+d2km(2,i,j)*faa(j)
 254        continue
 255     continue
      end if
c
c----------------------------------------------------------------------
c     Set F_g and F_A tensors in the diffusion frame
c----------------------------------------------------------------------
      if (itd .eq. 0) then
c                           *** no tilt: copy tensors directly
         do 260 j=1,5
         fgd(1,j)=fgm(j)
           do 259 i=1,2
            fad(i,j)=fam(i,j)
  259       continue
  260    continue
c
      else
c                    *** Transform A and g tensors into the diffusion frame
c
         call cd2km(d2km,ald,bed,gad)
         do 265 i=1,5
            do 264 j=1,5
               fgd(1,i)=fgd(1,i)+d2km(1,i,j)*fgm(j)
               fgd(2,i)=fgd(2,i)+d2km(2,i,j)*fgm(j)
               fad(1,i)=fad(1,i)+d2km(1,i,j)*fam(1,j)
     #                          -d2km(2,i,j)*fam(2,j)
               fad(2,i)=fad(2,i)+d2km(1,i,j)*fam(2,j)
     #                          +d2km(2,i,j)*fam(1,j)
 264        continue
 265     continue
      end if
c
c----------------------------------------------------------------------
c     Relaxation parameters & Heisenberg exchange parameter
c----------------------------------------------------------------------
c         
c NOTE: in NLS version, these rates are given in powers of ten
c
      if (t2edi.gt.rndoff) t2edi=ten**t2edi
      if (t2ndi.gt.rndoff) t2ndi=ten**t2ndi
      if (t2efi.gt.rndoff) t2efi=ten**t2efi
      if (t1edi.gt.rndoff) t1edi=ten**t1edi
      if (t1ndi.gt.rndoff) t1ndi=ten**t1ndi
      if (oss.gt.rndoff) oss=ten**oss
c
c----------------------------------------------------------------------
c     Check basis set parameters
c----------------------------------------------------------------------
c                             ** basis already set **
      if (setbas) go to 500
c
      inlemx = lemx
      inlomx = lomx
      inkmx  = kmx
      inmmx  = mmx
      inpnmx = ipnmx
c
      if((lemx.gt.mxlval)) lemx=mxlval
      if(ipar(lemx).ne.1) lemx=lemx-1
      if(ipar(lomx).ne.-1) lomx=lomx-1
      if(lomx.gt.lemx) lomx=lemx-1
      if(kmx.gt.lemx) kmx=lemx
      if(mmx.gt.lemx) mmx=lemx
      if(ipar(kmx).ne.1) kmx=kmx-1
      if(ipnmx.gt.in2) ipnmx=in2
      if((ipsi0.eq.0).and.(mmx.gt.ipnmx+1)) mmx=ipnmx
c
      if(lemx.lt.0)  lemx=0
      if(lomx.lt.0)  lomx=0
      if(kmx.lt.0)   kmx=0
      if(mmx.lt.0)   mmx=0
      if(ipnmx.lt.0) ipnmx=0
c
      if (inlemx .ne. lemx .or.
     #    inlomx .ne. lomx .or.
     #    inkmx  .ne. kmx  .or.
     #    inmmx  .ne. mmx  .or.
     #    inpnmx .ne. ipnmx ) then
c
         if (lumsg.ne.0) write (lumsg,1009) lemx,lomx,kmx,mmx,ipnmx
         ierr=9
      endif
c
c----------------------------------------------------------------------
c  Determine basis set using rules in M,I,I,M & Freed
c  (1) jm=1 (use only symmetric M combinations) for no nuclear Zeeman
c  (2) jkmn=1 (use only symmetric K combinations) if the magnetic 
c          tensors in the diffusion frame are real-valued. This is the
c          case if alpha_m, gamma_m (magnetic tilt), and gamma_d 
c          (diffusion tilt) are all zero.
c  (3) only even K values if there is no magnetic and diffusion tilt.
c  (4) only even L values no K values (kmx=0) in case of axial magnetic
c          tensors, axial potential, and no magnetic/diffusion tilt.
c----------------------------------------------------------------------
      jmmn=1
      if (abs(zeen).gt.rndoff) jmmn=-1
c
      jkmn=1
      do 270 i=1,5
         if (    dabs(fgd(2,i)).ge.rndoff 
     #      .or. dabs(fad(2,i)).ge.rndoff ) jkmn=-1
 270  continue
c
      if (itm.eq.0 .and.itd.eq.0) then
         kdelta=2
      else
         kdelta=1
      end if
c
      if (axiala .and. axialg .and. kdelta.eq.2 .and. kptmx.eq.0)
     # then
         ldelta=2
         lomx=0
         kmx=0
      else
         ldelta=1
      end if
c
c----------------------------------------------------------------------
c     check the calculation type and related information
c----------------------------------------------------------------------
 500  continue
c
      if (nstep.gt.MXSTEP-2) then
        if (lumsg.ne.0) write (lumsg,*) 
     #		'adjusting nstep to ',MXSTEP-2
	write(*,*)'adjusting nstep to ',MXSTEP-2
	nstep=MXSTEP-2
      end if
c
      if ((itype.ne.1).and.(itype.ne.4)) then
         if (lumsg.ne.0) write (lumsg,1016)
         ierr=16
         itype=4
      end if
c                                   *** eprbf calculation ***
      if (itype.eq.2) then
         if (nstep.eq.0) then
           if(ndimo.gt.0) then
             if (lumsg.ne.0) write (lumsg,1017)ndimo
             ierr=17
             nstep=ndimo
           else
             if (lumsg.ne.0) write (lumsg,1170)
             ierr=-17
             return
           end if  
         end if
c
         if (ptol.lt.rndoff) then
            if (lumsg.ne.0) write (lumsg,1018)
            ierr=18
            ptol=1.0d-2
         end if
c
         if ((abs(fieldi).lt.1.0d0).and.(abs(fieldi).lt.1.0d0)) then
            if (lumsg.ne.0) write (lumsg,1019)
            ierr=19
            fieldi=-50.0d0
            fieldf=-50.0d0
         end if
c
         if (nfield.lt.2) then
            if (lumsg.ne.0) write (lumsg,1020)
            ierr=20
            nfield=10
         end if
c
      end if
c                                   *** reset ierr flag ***
      if (ierr.gt.0) ierr=0
      return
c
c======================================================================
c
c Formats for warning and error messages
c----------------------------------------------------------------------
c Warning messages
c
 1001 format(' Nuclear spin in2=0 assumed')
 1002 format(' Positive B0 value assumed')
 1003 format(' Warning: high-field approximation may not apply for',
     #       ' given parameters')
 1004 format(' No non-Brownian motion specified for l, xy, zz:',
     #       ' ipdf=0 assumed')
 1005 format(' Nonzero potential with jump/free model:',
     #       ' ipdf=0 assumed') 
 1006 format(' Zero potential with anisotropic viscosity:',
     #       ' ipdf=0 assumed')
 1007 format(' Discrete jumps with anisotropic viscosity:',
     #       ' ipdf=0 assumed')
 1008 format(' lemx must be 48 or less with potential;',
     #       ' lemx=48 assumed')
 1009 format(' Basis set parameters adjusted to ',4(i3,','), i3)
 1010 format(' Too many Lanczos/CG steps specified: nstep=',i3,
     #       ' assumed ')
 1011 format(' Specified cgtol too small: cgtol=',g9.3,' assumed')
 1012 format(' Specified shiftr too small: shiftr=',g9.3,' assumed')
 1013 format(' Axial g-tensor: alm=0 assumed')
 1014 format(' Axial g, hf-tensors: alm=gam=gad=0 assumed')
 1015 format(' Zero potential with nort > 1 attempted')
 1016 format(' Calculation type: Rutishauser assumed')
 1017 format(' Zero CG steps: nstep=ndimo assumed',I8)
 1170 format(' Zero CG steps: fatel error in pcheck')
 1018 format(' Zero pruning tolerance for EPRBF: ptol=0.01 assumed') 
 1019 format(' Zero field range for EPRBF: -50 -- +50 assumed')
 1020 format(' Zero field points for EPRBF: nfield=10 assumed')
c
c----------------------------------------------------------------------
c Fatal error messages
c
 2001 format(' Fatal parameter error: zero B0')
 2002 format(' Fatal parameter error: zero values in g-tensor')
c
c======================================================================
c
      end


c----------------------------------------------------------------------
c                      =========================
c                        subroutine PRMSOK
c                      =========================
c
c  Checks all the parameter blocks up to the currently defined 
c  number of spectra and sites to see that:
c
c   (1) B0 is nonzero
c   (2) gxx,gyy,gzz are all nonzero (Cartesian tensor) or
c       at least one of g1,g2,g3 is nonzero (spherical tensor)
c   (3) Dimension of the matrix is valid
c   (4) nort is valid
c
c   Returns .false. if any of these conditions is violated for any
c   of the given sites.
c----------------------------------------------------------------------
      function prmsOK( lu )
      implicit none
      logical prmsOK
      integer lu
c
      include 'limits.inc'
      include 'simparm.inc'
      include 'datas.inc'
      include 'parms.inc'
      include 'miscel.inc'
      include 'basis.inc'
      include 'rndoff.inc'
c
      integer j
      logical allOK,gOK,dimOK,dimdOK,ortOK
c
      integer CARTESIAN,SPHERICAL,AXIAL,isite,ispec
      parameter(CARTESIAN=1,SPHERICAL=2,AXIAL=3)
c
      gOK=.true.
      allOK=.true.
      dimOK=.true.
      dimdOK=.true.
      ortOK=.true.
c total basis set storage check:
      dimOK = dimOK .and. ((pidptr(nbasis)+ndimoss(nbasis)-1.gt.0)
     #  .and. (pidptr(nbasis)+ndimoss(nbasis)-1.le.MMXDIM))
      dimdOK = dimdOK .and. ((dpidptr(nbasis)+ndimdss(nbasis)-1.gt.0)
     #  .and. (dpidptr(nbasis)+ndimdss(nbasis)-1.le.MMXDIM))
c
c  for all sites/spectra:
c
      do 10 ispec=1,nspectra
 	do 10 isite=1,ncomps
          allOK = allOK .and. (dabs(fparm(IB0,ispec,isite)).gt.
     #		rndoff)
c
c   Check g-tensors. For Cartesian tensor, all elements must be nonzero.
c   For spherical tensor, at least one element must be nonzero.
c
          if (iparm(IIGFLG,ispec,isite).eq.CARTESIAN) then
            gOK = gOK .and. (dabs(fparm(IGXX,ispec,isite)).gt.
     #        rndoff .and. dabs(fparm(IGXX+1,ispec,isite)).gt.
     #        rndoff .and. dabs(fparm(IGXX+2,ispec,isite)).gt.rndoff)
          else if (iparm(IIGFLG,ispec,isite).eq.SPHERICAL) then
            gOK = gOK .and. (dabs(fparm(IGXX,ispec,isite)).gt.
     #        rndoff .or. dabs(fparm(IGXX+1,ispec,isite)).gt.
     #        rndoff .or. dabs(fparm(IGXX+2,ispec,isite)).gt.rndoff)
          else if (iparm(IIGFLG,ispec,isite).eq.AXIAL) then
            gOK = gOK .and. (dabs(fparm(IGXX,ispec,isite)).gt.
     #        rndoff .and. dabs(fparm(IGXX+1,ispec,isite)).gt.rndoff)
          end if
c
          allOK = allOK .and. gOK
c  Individual basis size check:
          dimOK = dimOK .and. ((ndimoss(basinfo(1,ispec,isite)).gt.0)
     #	    .and. (ndimoss(basinfo(1,ispec,isite)).le.mxdim))
          dimdOK = dimdOK.and.((ndimdss(basinfo(1,ispec,isite)).gt.0)
     #      .and. (ndimoss(basinfo(1,ispec,isite)).le.mxdim))
          allOK = allOK .and. dimOK .and. dimdOK
c
          ortOK = ortOK .and. ((iparm(INORT,ispec,isite).gt.0).and.
     #		(iparm(INORT,ispec,isite).le.mxort))
          allOK = allOK .and. ortOK
 10   continue
c
      if (.not.allOK) write (lu,1000)
      if (.not.gOK) write (lu,1001)
      if (diaflg) dimOK = dimOK .and. dimdOK
      if (.not.dimOK) write (lu,1002)
      if (.not.ortOK) write (lu,1003)
c
      prmsOK=allOK
      return
c
 1000 format('*** B0 is zero for spectrum ***')
 1001 format('*** Zero g-tensor value set ***')
 1002 format('*** Dimension of the matrix is not valid ***')
 1003 format('*** Number of orientations for MOMD is not valid ***')
      end


c----------------------------------------------------------------------
c                      =========================
c                        subroutine EVALOK
c                      =========================
c
c  Checks if real part of all the eigenvalues is positive.
c
c----------------------------------------------------------------------
      function evalOK( ips,eval,nvar,lu )
c
      implicit none
      logical evalOK
      integer ips,nvar,lu,luold,i,j,nev,ii
c
      include 'limits.inc'
      include 'stdio.inc'
      include 'simparm.inc'
      include 'parms.inc'
      include 'rndoff.inc'
c
      double precision tol
      complex*16 eval(mxdim)      
c
      evalOK=.true.
      tol=rndoff*1.0d9
      luold=lu
      nev=nevo
      if (ips.eq.0) nev=nevd
c
c RC - Modify code to set neg eigenvalues to zero and continue:
c
c      do 100 i=1,nev
c         if ( dreal(eval(i)).lt.-tol ) then
c            if (ips.eq.0) then
c               write (lu,1001) psi
c               if (lu.ne.luttyo) write (luttyo,1001) psi
c            else
c               write (lu,1002) psi
c               if (lu.ne.luttyo) write (luttyo,1002) psi
c            end if
c            write (lu,1003) (j,eval(j),j=1,nev)
c            if (nvar.gt.0) then
c               write (lu,1004)
c               write (lu,1005) (tag(j),fparm(ixpr(j)),j=1,nvar)
c            end if
c            go to 20
c         end if
c 100  continue
c      return
c
c 20   negegv=negegv+1
c      evalOK=.false.
c      lu=luold
c      return
c
c             write(*,*)(ii,eval(ii),ii=1,nev)
      do 100 i=1,nev
         if ( dreal(eval(i)).lt.-tol ) then
           if (ips.eq.0) then
             write (lu,2001) i,psi,dreal(eval(i)),
     #          dimag(eval(i))
             if (lu.ne.luttyo) write (luttyo,2001) i,psi,dreal(eval(i)),
     #		dimag(eval(i))
c             write(*,*)(ii,eval(ii),ii=1,nev)
             eval(i)=(0.0d0,0.0d0)
c             eval(i)=(0.0d0,dimag(eval(i)))
           else
             write (lu,2002) i,psi,dreal(eval(i)),
     #          dimag(eval(i))
             if (lu.ne.luttyo) write (luttyo,2002) i,psi,dreal(eval(i)),
     #          dimag(eval(i))
c             write(*,*)(ii,eval(ii),ii=1,nev)
             eval(i)=(0.0d0,0.0d0)
c             eval(i)=(0.0d0,dimag(eval(i)))
           end if
         end if
 100   continue
      return
c
c## format statements ###############################################
c
 1001 format(/5x,'*** Negative real part in psi = ',f5.2,
     #' of DIAGONAL eigenvalues ***')
 1002 format(/5x,'*** Negative real part in psi = ',f5.2,
     #' of OFFDIAGONAL eigenvalues ***')
 2001 format(/5x,'*** Neg. real at psi = ',i7,f6.2,
     #' of DIAG. egvs. *** ',2g15.5,/,'Forced to zero')
 2002 format(/5x,'*** Neg. real at psi = ',i7,f6.2,
     #' of OFFDIAG. egvs. *** ',2g15.5,/,'Forced to zero')
 1003 format(8x,i4,4x,g14.7,2x,g14.7)
 1004 format(/8x,'EPR parameters')
 1005 format(11x,a,4x,g16.9)
c
      end
