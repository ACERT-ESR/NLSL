c Version 1.5.1 abeta 03/07/96
c----------------------------------------------------------------------
c                    =========================
c                        subroutine LCHECK
c                    =========================
c
c      This procedure replaces the function of program LBLL for 
c      performing initial computations based on the input parameters
c      for the EPR lineshape programs. It loads the floating point and
c      integer variables in the fparm and iparm arrays into common block
c      /eprprm/ and checks to make sure they are valid.
c
c      Returns ierr=0 if parameters are valid;  a negative number
c      identifying a fatal error, and a positive number indicating 
c      other conditions.
c
c      lumsg specifies a logical unit for output of any error messages.
c      None are produced if lumsg=0.
c
c     Note: Check of MTS indices that was originally coded here
c           is now in file mtschk
c
c      Fatal parameter errors
c       1) zero gxx, gyy, or gzz
c       2) zero B0
c      
c       Includes:
c               nlsdim.inc
c               maxl.inc
c               rndoff.inc
c               eprprm.inc
c
c       Uses:
c               ipar
c               cd2km
c               sp2car 
c
c----------------------------------------------------------------------
      subroutine lcheck(ierr)
c
      use nlsdim
      use maxl
      use eprprm
c      use prmeqv
      use errmsg
      use rnddbl
c
      implicit none
      integer ierr
c
      integer i,itmp,j 
      integer inlemx,inlomx,inkmn,inkmx,inmmn,inmmx,inpnmx
      double precision gmax,gmin,hmax,dn
      double precision d2km(2,5,5)
      logical axiala,axialg,axialw,axialr
c
      double precision ZERO,HALF,DSQ23,ONE,TEN
      parameter (ZERO=0.0D0,HALF=0.5D0,ONE=1.0D0,TEN=1.0D1)
      parameter (DSQ23=0.816496580927726d0)
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
      do i=1,5
         faa(i)=ZERO
         fgm(i)=ZERO
         fwm(i)=ZERO
         do j=1,2
            fam(j,i)=ZERO
            fgd(j,i)=ZERO
            fad(j,i)=ZERO
            fwd(j,i)=ZERO
         end do
      end do

c----------------------------------------------------------------------
c     Check for zero field
c----------------------------------------------------------------------
      if (b0 .lt. ZERO) then
         ierr=POSB0
         b0=abs(b0)
      elseif (b0 .eq. ZERO) then
         ierr=ZEROB0
         return
      endif
c----------------------------------------------------------------------
c     Check for nonzero number of points and field range
c----------------------------------------------------------------------
      if (nfld.le.0) then
         ierr=ZERONFLD
         return
      end if
c
c----------------------------------------------------------------------
c     Convert g, A, and R tensors to Cartesian form
c----------------------------------------------------------------------
      call tocart( fepr(IGXX), igflg )
      call tocart( fepr(IAXX), iaflg )
      call tocart( fepr(IWXX), iwflg )
      call tocart( fepr(IDX),  irflg )
c
c  *** NOTE: in NLS version, diffusion rates are specified as
c            a power of 10
c
      dx=TEN**dx
      dy=TEN**dy
      dz=TEN**dz
      axialr=abs(dx-dy).le.RNDOFF
c
c----------------------------------------------------------------------
c     Check for zero values in g-tensor and zero A tensor
c----------------------------------------------------------------------
      if ( abs(gxx) .lt. RNDOFF  .or.
     #     abs(gyy) .lt. RNDOFF  .or.
     #     abs(gzz) .lt. RNDOFF)   then
        ierr=ZEROG
        return
      endif

      if ( (in2 .le. 0) .or.
     #    (abs(axx) .lt. RNDOFF .and. 
     #     abs(ayy) .lt. RNDOFF .and.
     #     abs(azz) .lt. RNDOFF)     ) then
        in2=0
        axx=ZERO
        ayy=ZERO
        azz=ZERO
        ierr=ZEROIN2
      endif

c----------------------------------------------------------------------
c     Calculate spherical components of F_g and F_A
c----------------------------------------------------------------------
c                                        *** Electronic Zeeman ***
      g0=(gxx+gyy+gzz)/3.0D0
      fgm(1)=HALF*(gxx-gyy)*(b0/g0)
      fgm(3)=DSQ23*(gzz-HALF*(gxx+gyy))*(b0/g0)
      fgm(5)=fgm(1)
      axialg=abs(fgm(1)).lt.RNDOFF
c                                        *** Hyperfine ***
      a0=(axx+ayy+azz)/3.0D0
      faa(1)=HALF*(axx-ayy)
      faa(3)=DSQ23*(azz-HALF*(axx+ayy))
      faa(5)=faa(1)
      axiala=abs(faa(1)).lt.RNDOFF
c                                        *** Line-broadening ***
      w0=(wxx+wyy+wzz)/3.0D0
      fwm(1)=HALF*(wxx-wyy)
      fwm(3)=DSQ23*(wzz-HALF*(wxx+wyy))
      fwm(5)=fwm(1)
      axialw=abs(fwm(1)).lt.RNDOFF
c
      gmin=dmin1( gxx, gyy, gzz )
      gmax=dmax1( gxx, gyy, gzz )
      hmax=dmax1( abs(axx), abs(ayy), abs(azz) )
c
c----------------------------------------------------------------------
c Issue warning if high-field approximation has been violated
c (this is a very rough criterion for the high-field approx.!)
c----------------------------------------------------------------------

      if (b0 .lt. 10.0D0*dmax1( hmax, (gmax-gmin)*b0/g0) ) then
         ierr=BADHIFLD
      endif

c----------------------------------------------------------------------
c  Set ipt and cpot array according to the potential coefficients given
c----------------------------------------------------------------------
      ipt=0
      do j=1,5
        do i=1,5
          cpot(i,j)=ZERO
       end do
      end do

      if (abs(c20) .gt. RNDOFF) then
         ipt=ipt+1
         cpot(2,1)=c20
         lptmx=2
      endif
      if (abs(c22) .gt. RNDOFF) then
         ipt=ipt+1
         cpot(2,2)=c22
         lptmx=2
         kptmx=2
      endif
      if (abs(c40) .gt. RNDOFF) then
         ipt=ipt+1
         cpot(3,1)=c40
         lptmx=4
      endif
      if (abs(c42) .gt. RNDOFF) then
         ipt=ipt+1
         cpot(3,2)=c42
         lptmx=4
         kptmx=2
      endif
      if (abs(c44) .gt. RNDOFF) then
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
      if (axialr) kband=kband+2

c----------------------------------------------------------------------
c Check for consistency of specified motion parameters with given model
c----------------------------------------------------------------------
      if (ipdf.gt.2) ipdf=2
      if (ipdf.lt.0) ipdf=0
      if (ist.lt.0) ist=0
c 
c  Set exponents according to non-Brownian model
c
      expl=ONE
      expkxy=ONE
      expkzz=ONE
      if (ml.eq.1)  expl=HALF
      if (mxy.eq.1) expkxy=HALF
      if (mzz.eq.1) expkzz=HALF
c
c
c     *** Non-Brownian motion: must have at least one nonzero residence
c         time, and no potential specified
c         
c
c    Set residence times to zero for Brownian model
c
      if (ipdf .eq. 0) then
         pml = ZERO
         pmxy = ZERO
         pmzz = ZERO
c
c NOTE: in NLS version, R*residence time products are specified as 
c       powers of ten
c
      else
         pml = TEN**pml
         pmxy = TEN**pmxy
         pmzz = TEN**pmzz
c 
c  *** NOTE: in NLS version, djf, djfprp specified as powers of ten
c
         djf=TEN**djf
         djfprp=TEN**djfprp
c
         if (ipt .gt. 0) then
            ierr=BADJMP
            ipdf=0
         endif
      endif
c
c     *** Anisotropic viscosity model: must have potential and
c         no discrete jump motion specified
c
      if (ipdf .eq. 2) then
         if (ipt .eq. 0) then
            ierr=BADAV
            ipdf=0
         endif
c
         if (ist .gt. 0) then
            ierr=BADDJ
            ipdf=0
         endif
      endif
c
      if (oss.ne.ZERO) oss=TEN**oss
c     
      ipsi0=0
      if( ((abs(psi).gt.RNDOFF).and.(abs(psi-180.0d0).gt.RNDOFF))
     #   .or. nort.gt.1) ipsi0=1
c
c----------------------------------------------------------------------
c     Set magnetic tilt flag and A tensor in g frame
c----------------------------------------------------------------------
c
c     --- mag. tilt angle gamma not needed if g,W-matrices are axial
c
c      if (axialg .and. axialw .and. abs(gam).gt.RNDOFF) then
c        gam=ZERO
cc        ierr=NOALPHAM
c      end if
c
      if (abs(alm).gt.RNDOFF .or. bem.gt.RNDOFF 
     #    .or. abs(gam).gt.RNDOFF) itm=1
c
c
      if (itm .eq. 0) then
c                           *** no tilt: copy A tensor directly
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
c     Set diffusion tilt flag and g,A,W tensors in diffusion frame
c----------------------------------------------------------------------
c
c      if (axialr .and. abs(gad).gt.RNDOFF) then
c        gad=ZERO
cc        ierr=NOALPHAD
c      end if
c
      if ( abs(ald).gt.RNDOFF .or.
     #     abs(gad).gt.RNDOFF .or. 
     #     bed.gt.RNDOFF) itd=1
c
c                           *** no tilt: copy tensors directly
      if (itd .eq. 0) then
         do i=1,5
            fgd(1,i)=fgm(i)
            fwd(1,i)=fwm(i)
            fad(1,i)=fam(1,i)
            fad(2,i)=fam(2,i)
         end do
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
     #                       +d2km(2,i,j)*fam(1,j)
            end do
         end do
      end if
c
c    --- Allow antisymmetric K combinations if any tensor has
c        imaginary elements (nonzero alm, gam, ald, or gad)
c
      jkmn=1
      do i=1,5
         if (     abs(fgd(2,i)).ge.RNDOFF 
     #       .or. abs(fad(2,i)).ge.RNDOFF
     #       .or. abs(fwd(2,i)).ge.RNDOFF ) jkmn=-1
      end do
c
c    --- Allow antisymmetric M combinations if there is a nonzero
c        nuclear Zeeman interaction
c
      if (abs(gamman).ge.RNDOFF .and. in2.gt.0) then
         jmmn=-1
      else
         jmmn=1
      end if

c----------------------------------------------------------------------
c     Check basis set parameters
c----------------------------------------------------------------------
      inlemx=lemx
      inlomx=lomx
      inkmn=kmn
      inkmx =kmx
      inmmn=mmn
      inmmx =mmx
      inpnmx=ipnmx
c
      if((lemx.gt.MXLVAL).and.(ipt.ne.0)) then
         ierr=LEMXHI
         lemx=MXLVAL
      endif
c
      if(ipar(lemx).ne.1) lemx=lemx-1
      if(ipar(lomx).ne.-1) lomx=lomx-1
      if(lomx.gt.lemx) lomx=lemx-1
      if(kmx.gt.lemx) kmx=lemx
      if(mmx.gt.lemx) mmx=lemx
      if(ipar(kmx).ne.1) kmx=kmx-1
      if(ipnmx.gt.in2) ipnmx=in2
      if((ipsi0.eq.0).and.(mmx.gt.ipnmx)) mmx=ipnmx
      if(-kmn.gt.kmx) kmn=-kmx
      if(-mmn.gt.mmx) mmn=-mmx
      if(jkmn.eq.1) kmn=max(kmn,0)
      if(jmmn.eq.1) mmn=max(kmn,0)
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
     #    inkmn  .ne. kmn  .or.
     #    inmmx  .ne. mmx  .or.
     #    inmmn  .ne. mmn  .or.
     #    inpnmx .ne. ipnmx ) then
c
         ierr=BSSADJ
         write (eprerr(BSSADJ)(24:),'(6(i3,'',''),i2)') 
     #         lemx,lomx,kmn,kmx,mmn,mmx,ipnmx
      endif
c
c----------------------------------------------------------------------
c  Determine basis set using rules in M,I,I,M & F
c  (1)   If there is no diffusion or magnetic tilt, only even K values
c        are needed
c----------------------------------------------------------------------
      if(itm.eq.0 .and. itd.eq.0) then
        kdelta=2
      else
        kdelta=1
      end if
c
c----------------------------------------------------------------------
c  (2)   In addition, if the magnetic and linewidth tensors are axial
c        (in the diffusion frame) only even L values and no K values
c         are needed
c----------------------------------------------------------------------
      if (axiala .and. axialg .and. axialw .and. 
     #    kdelta.eq.2 .and. kptmx.eq.0) then
         ldelta=2
         kmx=0
      else
         ldelta=1
      end if
c
c----------------------------------------------------------------------
c     check number of Lanczos/CG steps to execute 
c----------------------------------------------------------------------
c
      if (nstep.gt.MXSTEP) then
         ierr=NSTEPHI
         write (eprerr(NSTEPHI)(30:),'(i5)') MXSTEP
         nstep=MXSTEP
      endif
c
c----------------------------------------------------------------------
c     Check specified tolerances and shifts for CG calculation
c----------------------------------------------------------------------
c
      if (abs(cgtol) .lt. TEN*RNDOFF) then
         cgtol=TEN*RNDOFF
         ierr=CGTOLHI
         write (eprerr(CGTOLHI)(30:),'(g9.3)') shiftr
      end if
      if (abs(shiftr) .lt. TEN*RNDOFF) then
         shiftr=10.0D0*RNDOFF
         ierr=SHIFTHI
         write (eprerr(SHIFTHI)(31:),'(g9.3)') shiftr
      end if        
c
      return
c
      end

