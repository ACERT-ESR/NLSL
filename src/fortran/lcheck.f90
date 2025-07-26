! Version 1.5.1 abeta 03/07/96
!----------------------------------------------------------------------
!                    =========================
!                        subroutine LCHECK
!                    =========================
!
!      This procedure replaces the function of program LBLL for 
!      performing initial computations based on the input parameters
!      for the EPR lineshape programs. It loads the floating point and
!      integer variables in the fparm and iparm arrays into common block
!      /eprprm/ and checks to make sure they are valid.
!
!      Returns ierr=0 if parameters are valid;  a negative number
!      identifying a fatal error, and a positive number indicating 
!      other conditions.
!
!      lumsg specifies a logical unit for output of any error messages.
!      None are produced if lumsg=0.
!
!     Note: Check of MTS indices that was originally coded here
!           is now in file mtschk
!
!      Fatal parameter errors
!       1) zero gxx, gyy, or gzz
!       2) zero B0
!      
!       Includes:
!               nlsdim.inc
!               maxl.inc
!               rndoff.inc
!               eprprm.inc
!
!       Uses:
!               ipar
!               cd2km
!               sp2car 
!
!----------------------------------------------------------------------
      subroutine lcheck(ierr)
!
      use nlsdim
      use maxl
      use eprprm
!      use prmeqv
      use errmsg
      use rnddbl
!
      implicit none
      integer ierr
!
      integer i,j
      integer inlemx,inlomx,inkmn,inkmx,inmmn,inmmx,inpnmx
      double precision gmax,gmin,hmax
      double precision d2km(2,5,5)
      logical axiala,axialg,axialw,axialr
!
      double precision ZERO,HALF,DSQ23,ONE,TEN
      parameter (ZERO=0.0D0,HALF=0.5D0,ONE=1.0D0,TEN=1.0D1)
      parameter (DSQ23=0.816496580927726d0)
!
      integer ipar
      external ipar
!

!######################################################################
!
      ierr=0
!
!----------------------------------------------------------------------
!     initialize spherical tensors
!----------------------------------------------------------------------
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

!----------------------------------------------------------------------
!     Check for zero field
!----------------------------------------------------------------------
      if (b0 .lt. ZERO) then
         ierr=POSB0
         b0=abs(b0)
      elseif (b0 .eq. ZERO) then
         ierr=ZEROB0
         return
      endif
!----------------------------------------------------------------------
!     Check for nonzero number of points and field range
!----------------------------------------------------------------------
      if (nfld.le.0) then
         ierr=ZERONFLD
         return
      end if
!
!----------------------------------------------------------------------
!     Convert g, A, and R tensors to Cartesian form
!----------------------------------------------------------------------
      call tocart( fepr(IGXX), igflg )
      call tocart( fepr(IAXX), iaflg )
      call tocart( fepr(IWXX), iwflg )
      call tocart( fepr(IDX),  irflg )
!
!  *** NOTE: in NLS version, diffusion rates are specified as
!            a power of 10
!
      dx=TEN**dx
      dy=TEN**dy
      dz=TEN**dz
      axialr=abs(dx-dy).le.RNDOFF
!
!----------------------------------------------------------------------
!     Check for zero values in g-tensor and zero A tensor
!----------------------------------------------------------------------
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

!----------------------------------------------------------------------
!     Calculate spherical components of F_g and F_A
!----------------------------------------------------------------------
!                                        *** Electronic Zeeman ***
      g0=(gxx+gyy+gzz)/3.0D0
      fgm(1)=HALF*(gxx-gyy)*(b0/g0)
      fgm(3)=DSQ23*(gzz-HALF*(gxx+gyy))*(b0/g0)
      fgm(5)=fgm(1)
      axialg=abs(fgm(1)).lt.RNDOFF
!                                        *** Hyperfine ***
      a0=(axx+ayy+azz)/3.0D0
      faa(1)=HALF*(axx-ayy)
      faa(3)=DSQ23*(azz-HALF*(axx+ayy))
      faa(5)=faa(1)
      axiala=abs(faa(1)).lt.RNDOFF
!                                        *** Line-broadening ***
      w0=(wxx+wyy+wzz)/3.0D0
      fwm(1)=HALF*(wxx-wyy)
      fwm(3)=DSQ23*(wzz-HALF*(wxx+wyy))
      fwm(5)=fwm(1)
      axialw=abs(fwm(1)).lt.RNDOFF
!
      gmin=dmin1( gxx, gyy, gzz )
      gmax=dmax1( gxx, gyy, gzz )
      hmax=dmax1( abs(axx), abs(ayy), abs(azz) )
!
!----------------------------------------------------------------------
! Issue warning if high-field approximation has been violated
! (this is a very rough criterion for the high-field approx.!)
!----------------------------------------------------------------------

      if (b0 .lt. 10.0D0*dmax1( hmax, (gmax-gmin)*b0/g0) ) then
         ierr=BADHIFLD
      endif

!----------------------------------------------------------------------
!  Set ipt and cpot array according to the potential coefficients given
!----------------------------------------------------------------------
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
!
      if (lptmx.ge.2) then
         lband=lptmx*2
      else
         lband=2
      end if
!
      kband=kptmx*2
      if (axialr) kband=kband+2

!----------------------------------------------------------------------
! Check for consistency of specified motion parameters with given model
!----------------------------------------------------------------------
      if (ipdf.gt.2) ipdf=2
      if (ipdf.lt.0) ipdf=0
      if (ist.lt.0) ist=0
! 
!  Set exponents according to non-Brownian model
!
      expl=ONE
      expkxy=ONE
      expkzz=ONE
      if (ml.eq.1)  expl=HALF
      if (mxy.eq.1) expkxy=HALF
      if (mzz.eq.1) expkzz=HALF
!
!
!     *** Non-Brownian motion: must have at least one nonzero residence
!         time, and no potential specified
!         
!
!    Set residence times to zero for Brownian model
!
      if (ipdf .eq. 0) then
         pml = ZERO
         pmxy = ZERO
         pmzz = ZERO
!
! NOTE: in NLS version, R*residence time products are specified as 
!       powers of ten
!
      else
         pml = TEN**pml
         pmxy = TEN**pmxy
         pmzz = TEN**pmzz
! 
!  *** NOTE: in NLS version, djf, djfprp specified as powers of ten
!
         djf=TEN**djf
         djfprp=TEN**djfprp
!
         if (ipt .gt. 0) then
            ierr=BADJMP
            ipdf=0
         endif
      endif
!
!     *** Anisotropic viscosity model: must have potential and
!         no discrete jump motion specified
!
      if (ipdf .eq. 2) then
         if (ipt .eq. 0) then
            ierr=BADAV
            ipdf=0
         endif
!
         if (ist .gt. 0) then
            ierr=BADDJ
            ipdf=0
         endif
      endif
!
      if (oss.ne.ZERO) oss=TEN**oss
!     
      ipsi0=0
      if( ((abs(psi).gt.RNDOFF).and.(abs(psi-180.0d0).gt.RNDOFF))
     #   .or. nort.gt.1) ipsi0=1
!
!----------------------------------------------------------------------
!     Set magnetic tilt flag and A tensor in g frame
!----------------------------------------------------------------------
!
!     --- mag. tilt angle gamma not needed if g,W-matrices are axial
!
!      if (axialg .and. axialw .and. abs(gam).gt.RNDOFF) then
!        gam=ZERO
!c        ierr=NOALPHAM
!      end if
!
      if (abs(alm).gt.RNDOFF .or. bem.gt.RNDOFF 
     #    .or. abs(gam).gt.RNDOFF) itm=1
!
!
      if (itm .eq. 0) then
!                           *** no tilt: copy A tensor directly
         do j=1,5
            fam(1,j)=faa(j)
         end do
      else
!                          *** transform A tensor into g axis system
!
         call cd2km(d2km,alm,bem,gam)
         do i=1,5
            do j=1,5
               fam(1,i)=fam(1,i)+d2km(1,i,j)*faa(j)
               fam(2,i)=fam(2,i)+d2km(2,i,j)*faa(j)
            end do
         end do
      end if
!
!----------------------------------------------------------------------
!     Set diffusion tilt flag and g,A,W tensors in diffusion frame
!----------------------------------------------------------------------
!
!      if (axialr .and. abs(gad).gt.RNDOFF) then
!        gad=ZERO
!c        ierr=NOALPHAD
!      end if
!
      if ( abs(ald).gt.RNDOFF .or.
     #     abs(gad).gt.RNDOFF .or. 
     #     bed.gt.RNDOFF) itd=1
!
!                           *** no tilt: copy tensors directly
      if (itd .eq. 0) then
         do i=1,5
            fgd(1,i)=fgm(i)
            fwd(1,i)=fwm(i)
            fad(1,i)=fam(1,i)
            fad(2,i)=fam(2,i)
         end do
      else
!                    *** Transform A, g, W tensors into the diffusion frame
!
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
!
!    --- Allow antisymmetric K combinations if any tensor has
!        imaginary elements (nonzero alm, gam, ald, or gad)
!
      jkmn=1
      do i=1,5
         if (     abs(fgd(2,i)).ge.RNDOFF 
     #       .or. abs(fad(2,i)).ge.RNDOFF
     #       .or. abs(fwd(2,i)).ge.RNDOFF ) jkmn=-1
      end do
!
!    --- Allow antisymmetric M combinations if there is a nonzero
!        nuclear Zeeman interaction
!
      if (abs(gamman).ge.RNDOFF .and. in2.gt.0) then
         jmmn=-1
      else
         jmmn=1
      end if

!----------------------------------------------------------------------
!     Check basis set parameters
!----------------------------------------------------------------------
      inlemx=lemx
      inlomx=lomx
      inkmn=kmn
      inkmx =kmx
      inmmn=mmn
      inmmx =mmx
      inpnmx=ipnmx
!
      if((lemx.gt.MXLVAL).and.(ipt.ne.0)) then
         ierr=LEMXHI
         lemx=MXLVAL
      endif
!
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
!
      if(lemx.lt.0)  lemx=0
      if(lomx.lt.0)  lomx=0
      if(kmx.lt.0)   kmx=0
      if(mmx.lt.0)   mmx=0
      if(ipnmx.lt.0) ipnmx=0
!
      if (inlemx .ne. lemx .or.
     #    inlomx .ne. lomx .or.
     #    inkmx  .ne. kmx  .or.
     #    inkmn  .ne. kmn  .or.
     #    inmmx  .ne. mmx  .or.
     #    inmmn  .ne. mmn  .or.
     #    inpnmx .ne. ipnmx ) then
!
         ierr=BSSADJ
         write (eprerr(BSSADJ)(24:),'(6(i3,'',''),i2)') 
     #         lemx,lomx,kmn,kmx,mmn,mmx,ipnmx
      endif
!
!----------------------------------------------------------------------
!  Determine basis set using rules in M,I,I,M & F
!  (1)   If there is no diffusion or magnetic tilt, only even K values
!        are needed
!----------------------------------------------------------------------
      if(itm.eq.0 .and. itd.eq.0) then
        kdelta=2
      else
        kdelta=1
      end if
!
!----------------------------------------------------------------------
!  (2)   In addition, if the magnetic and linewidth tensors are axial
!        (in the diffusion frame) only even L values and no K values
!         are needed
!----------------------------------------------------------------------
      if (axiala .and. axialg .and. axialw .and. 
     #    kdelta.eq.2 .and. kptmx.eq.0) then
         ldelta=2
         kmx=0
      else
         ldelta=1
      end if
!
!----------------------------------------------------------------------
!     check number of Lanczos/CG steps to execute 
!----------------------------------------------------------------------
!
      if (nstep.gt.MXSTEP) then
         ierr=NSTEPHI
         write (eprerr(NSTEPHI)(30:),'(i5)') MXSTEP
         nstep=MXSTEP
      endif
!
!----------------------------------------------------------------------
!     Check specified tolerances and shifts for CG calculation
!----------------------------------------------------------------------
!
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
!
      return
!
      end

