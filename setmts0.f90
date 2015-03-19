c Version 1.5.1 beta 3/21/96
c----------------------------------------------------------------------
c                    =========================
c                        function SETMTS   
c                    =========================
c
c  Given a set of slow-motional EPR parameters, check the MTS indices
c  and adjust them accordingly.
c
c  Inputs
c    fparm   Array of floating-point EPR parameters appearing in the
c            order defined in /eprprm/ 
c    iparm   Array of integer EPR parameters appearing in the
c            order defined in /eprprm/ 
c
c  Output
c    mts     Array to receive the corrected MTS parameters
c              Information is packed as documented in mtsdef.inc
c
c  Returns 0 for successful completion
c          8 if lemx exceeded limit for calculations w/potential
c          9 if mts parameters were changed from those in iparm array
c
c----------------------------------------------------------------------
      function setmts(fparm,iparm,mts)
      implicit none
c
      use nlsdim
      use eprprm
      use errmsg
      use mtsdef
      use maxl
      use rnddbl
c
      double precision fparm(NFPRM)
      integer setmts,iparm(NIPRM),mts(MXMTS)
c
      integer i
      logical knon0p,changed
c
      integer ipar
      logical isaxial
      external ipar,isaxial
c
c----------------------------------------------------------------------
c     If there is no associated list of basis set indices,
c     check basis set parameters
c----------------------------------------------------------------------
c
      setmts=0
      knon0p = abs(fparm(IC20+1)).gt.RNDOFF .or.
     #         abs(fparm(IC20+2)).gt.RNDOFF .or.
     #         abs(fparm(IC20+3)).gt.RNDOFF .or.
     #         abs(fparm(IC20+4)).gt.RNDOFF
c
      do i=1,NTRC
         mts(i)=iparm(ILEMX+i-1)
      end do
c
      mts(NIPSI0)=0 
      if ( (abs(fparm(IPSI)).gt.RNDOFF .and.
     #     abs(fparm(IPSI)-1.8D1).gt.RNDOFF) 
     #    .or. iparm(INORT).gt.1) mts(NIPSI0)=1
      mts(NIN2)=iparm(IIN2)
c
c    --- Allow antisymmetric K combinations if any tensor has
c        imaginary elements (nonzero alm, gam, ald, or gad)
c
      mts(NJKMN)=1
      if (abs(fparm(IALM)).ge.RNDOFF .or.
     #    abs(fparm(IGAM)).ge.RNDOFF .or.
     #    abs(fparm(IALD)).ge.RNDOFF .or.
     #    abs(fparm(IGAD)).ge.RNDOFF) mts(NJKMN)=-1
c
c    --- Allow antisymmetric M combinations if there is a nonzero
c        nuclear Zeeman interaction
c
      if (abs(fparm(IGAMAN)).ge.RNDOFF .and. iparm(IIN2).gt.0) then
         mts(NJMMN)=-1
      else
         mts(NJMMN)=1
      end if

c
      if((iparm(ILEMX).gt.MXLVAL).and.(ipt.ne.0)) then
         setmts=LEMXHI
         iparm(ILEMX)=MXLVAL
      endif
c
      if (mts(NLEMX).gt.MXLVAL) mts(NLEMX)=MXLVAL
      if (ipar(mts(NLEMX)).ne.1) mts(NLEMX)=mts(NLEMX)-1
      if (ipar(mts(NLOMX)).ne.-1) mts(NLOMX)=mts(NLOMX)-1
      if (mts(NLOMX).gt.mts(NLEMX)) mts(NLOMX)=mts(NLEMX)-1
      if (mts(NKMX).gt.mts(NLEMX)) mts(NKMX)=mts(NLEMX)
      if (mts(NMMX).gt.mts(NLEMX)) mts(NMMX)=mts(NLEMX)
      if (ipar(mts(NKMX)).ne.1) mts(NKMX)=mts(NKMX)-1
      if (mts(NIPNMX).gt.in2) mts(NIPNMX)=iparm(IIN2)
      if (mts(NIPSI0).eq.0.and.mts(NMMX).gt.mts(NIPNMX) )
     #    mts(NMMX)=mts(NIPNMX)
      if (mts(NJKMN).eq.1) mts(NKMN)=max(mts(NKMN),0)
      if (mts(NJMMN).eq.1) mts(NMMN)=max(mts(NMMN),0)
c     
      if(mts(NLEMX).lt.0)  mts(NLEMX)=0
      if(mts(NLOMX).lt.0)  mts(NLOMX)=0
      if(mts(NKMX).lt.0)   mts(NKMX)=0
      if(mts(NMMX).lt.0)   mts(NMMX)=0
      if(mts(NIPNMX).lt.0) mts(NIPNMX)=0
c
      changed=.false.
      do i=1,NTRC
         changed = changed .or.(mts(i).ne.iparm(ILEMX+i-1)) 
      end do
c
      if (changed) then
         setmts=BSSADJ
         write (eprerr(BSSADJ)(24:),'(6(i3,'',''),i2)') 
     #         lemx,lomx,kmn,kmx,mmn,mmx,ipnmx
      end if
c
c----------------------------------------------------------------------
c  Determine basis set using rules in M,I,I,M & F
c  (1)   If there is no diffusion tilt, only even K values are needed
c----------------------------------------------------------------------
      if(abs(ald).gt.RNDOFF .or. abs(bed).gt.RNDOFF .or. 
     #  (abs(bed)-1.8D1).gt.RNDOFF .or.abs(gad).gt.RNDOFF) then
         mts(NKDEL)=2
      else
         mts(NKDEL)=1
      end if
c
c----------------------------------------------------------------------
c  (2)   In addition, if the magnetic and linewidth tensors are axial, 
c        and there is no diffusion tilt or K.ne.0 potential terms,
c        only even L values and no K values are needed
c----------------------------------------------------------------------
      if( isaxial(fparm(IGXX),iparm(IIGFLG)) .and.
     #    isaxial(fparm(IAXX),iparm(IIAFLG)) .and.
     #    isaxial(fparm(IWXX),iparm(IIWFLG)) .and.
     #     (mts(NKDEL).eq.2).and.(.not.knon0p)) then
         mts(NLDEL)=2
         mts(NKMX)=0
      else
         mts(NLDEL)=1
      end if
c
      return
      end
