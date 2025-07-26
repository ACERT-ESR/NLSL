! Version 1.5.1 beta 3/21/96
!----------------------------------------------------------------------
!                    =========================
!                        function SETMTS   
!                    =========================
!
!  Given a set of slow-motional EPR parameters, check the MTS indices
!  and adjust them accordingly.
!
!  Inputs
!    fparmi   Array of floating-point EPR parameters appearing in the
!            order defined in /eprprm/ 
!    iparmi   Array of integer EPR parameters appearing in the
!            order defined in /eprprm/ 
!
!  Output
!    mts     Array to receive the corrected MTS parameters
!              Information is packed as documented in mtsdef.inc
!
!  Returns 0 for successful completion
!          8 if lemx exceeded limit for calculations w/potential
!          9 if mts parameters were changed from those in iparmi array
!
!----------------------------------------------------------------------
      function setmts(fparmi,iparmi,mts)
!
      use nlsdim
      use eprprm
      use errmsg
      use mtsdef
      use maxl
      use rnddbl
!
      implicit none
      double precision fparmi(NFPRM)
      integer setmts,iparmi(NIPRM),mts(MXMTS)
!
      integer i
      logical knon0p,changed
!
      integer ipar
      logical isaxial
      external ipar,isaxial
!
!----------------------------------------------------------------------
!     If there is no associated list of basis set indices,
!     check basis set parameters
!----------------------------------------------------------------------
!
      setmts=0
      knon0p = abs(fparmi(IC20+1)).gt.RNDOFF .or.
     #         abs(fparmi(IC20+2)).gt.RNDOFF .or.
     #         abs(fparmi(IC20+3)).gt.RNDOFF .or.
     #         abs(fparmi(IC20+4)).gt.RNDOFF
!
      do i=1,NTRC
         mts(i)=iparmi(ILEMX+i-1)
      end do
!
      mts(NIPSI0)=0 
      if ( (abs(fparmi(IPSI)).gt.RNDOFF .and.
     #     abs(fparmi(IPSI)-1.8D2).gt.RNDOFF) 
     #    .or. iparmi(INORT).gt.1) mts(NIPSI0)=1
      mts(NIN2)=iparmi(IIN2)
!
!    --- Allow antisymmetric K combinations if any tensor has
!        imaginary elements (nonzero alm, gam, ald, or gad)
!
      mts(NJKMN)=1
      if (abs(fparmi(IALM)).ge.RNDOFF .or.
     #    abs(fparmi(IGAM)).ge.RNDOFF .or.
     #    abs(fparmi(IALD)).ge.RNDOFF .or.
     #    abs(fparmi(IGAD)).ge.RNDOFF) mts(NJKMN)=-1
!
!    --- Allow antisymmetric M combinations if there is a nonzero
!        nuclear Zeeman interaction
!
      if (abs(fparmi(IGAMAN)).ge.RNDOFF .and.
     #    iparmi(IIN2).gt.0) then
         mts(NJMMN)=-1
      else
         mts(NJMMN)=1
      end if

!
      if((iparmi(ILEMX).gt.MXLVAL).and.(ipt.ne.0)) then
         setmts=LEMXHI
         iparmi(ILEMX)=MXLVAL
      endif
!
      if (mts(NLEMX).gt.MXLVAL) mts(NLEMX)=MXLVAL
      if (ipar(mts(NLEMX)).ne.1) mts(NLEMX)=mts(NLEMX)-1
      if (ipar(mts(NLOMX)).ne.-1) mts(NLOMX)=mts(NLOMX)-1
      if (mts(NLOMX).gt.mts(NLEMX)) mts(NLOMX)=mts(NLEMX)-1
      if (mts(NKMX).gt.mts(NLEMX)) mts(NKMX)=mts(NLEMX)
      if (mts(NMMX).gt.mts(NLEMX)) mts(NMMX)=mts(NLEMX)
      if (ipar(mts(NKMX)).ne.1) mts(NKMX)=mts(NKMX)-1
      if (mts(NIPNMX).gt.in2) mts(NIPNMX)=iparmi(IIN2)
      if (mts(NIPSI0).eq.0.and.mts(NMMX).gt.mts(NIPNMX) )
     #    mts(NMMX)=mts(NIPNMX)
      if (mts(NJKMN).eq.1) mts(NKMN)=max(mts(NKMN),0)
      if (mts(NJMMN).eq.1) mts(NMMN)=max(mts(NMMN),0)
!     
      if(mts(NLEMX).lt.0)  mts(NLEMX)=0
      if(mts(NLOMX).lt.0)  mts(NLOMX)=0
      if(mts(NKMX).lt.0)   mts(NKMX)=0
      if(mts(NMMX).lt.0)   mts(NMMX)=0
      if(mts(NIPNMX).lt.0) mts(NIPNMX)=0
!
      changed=.false.
      do i=1,NTRC
         changed = changed .or.(mts(i).ne.iparmi(ILEMX+i-1)) 
      end do
!
      if (changed) then
         setmts=BSSADJ
         write (eprerr(BSSADJ)(24:),'(6(i3,'',''),i2)') 
     #         lemx,lomx,kmn,kmx,mmn,mmx,ipnmx
      end if
!
!----------------------------------------------------------------------
!  Determine basis set using rules in M,I,I,M & F
!  (1)   If there is no diffusion tilt, only even K values are needed
!----------------------------------------------------------------------
      if(abs(ald).gt.RNDOFF .or. abs(bed).gt.RNDOFF .or. 
     #  (abs(bed)-1.8D2).gt.RNDOFF .or.abs(gad).gt.RNDOFF) then
         mts(NKDEL)=1
      else
         mts(NKDEL)=2
      end if
!
!----------------------------------------------------------------------
!  (2)   In addition, if the magnetic and linewidth tensors are axial, 
!        and there is no diffusion tilt or K.ne.0 potential terms,
!        only even L values and no K values are needed
!----------------------------------------------------------------------
      if( isaxial(fparmi(IGXX),iparmi(IIGFLG)) .and.
     #    isaxial(fparmi(IAXX),iparmi(IIAFLG)) .and.
     #    isaxial(fparmi(IWXX),iparmi(IIWFLG)) .and.
     #     (mts(NKDEL).eq.2).and.(.not.knon0p)) then
         mts(NLDEL)=2
         mts(NKMX)=0
      else
         mts(NLDEL)=1
      end if
!
      return
      end
