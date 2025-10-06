! NLSL Version 1.5.1 beta 2/3/96
!----------------------------------------------------------------------
!                    =========================
!                       subroutine SETSPC
!                    =========================
!
!  This routine sets all the parameters in the fparm and iparm array
!  (kept in common /parcom/) for the specified site and spectrum
!  (arguments isite and ispc).
!
!  Certain parameter values are not kept in the fparm and iparm
!  arrays. In particular, if a series of spectra are to be calculated,
!  the series variable appropriate for the current spectrum is loaded
!  into fparm. This routine also sets the spectrum-dependent parameters
!  from the arrays in common /expdat/.  These parameters are not allowed
!  to vary within a single spectrum, and specifically include the following:
!
!    Floating:
!       B0
!       LB
!       PHASE
!       PSI
!       FLDI
!       DFLD
!
!    Integer:
!       NFLD
!       IDERIV
!----------------------------------------------------------------------
      subroutine setspc(isite,ise)
!
      use nlsdim
      use eprprm
      use expdat
      use parcom
      use basis
      use errmsg
      use iterat
      use rnddbl
      use stdio
      use pidef
!
      implicit none
      integer isite,ise
!
      integer i
!
!######################################################################
!
      if (ise.le.0 .or. ise.gt.MXSPC
     #   .or. isite.le.0 .or. isite.gt.MXSITE) return
!
!     --------------------------
!     Set spectral parameters
!     --------------------------
!
      fparm(ILB,   isite) = slb(ise)
      fparm(IPHASE,isite) = sphs(ise)
      fparm(IPSI,  isite) = spsi(ise)
      fparm(IB0,   isite) = sb0(ise)
      fparm(IFLDI, isite) = sbi(ise)-shft(ise)-tmpshft(ise)
      fparm(IDFLD, isite) = sdb(ise)
      iparm(INFLD, isite) = npts(ise)
      iparm(IIDERV,isite) = idrv(ise)
!
!     -------------------------------------- 
!     set the SERIES variable if necessary
!     -------------------------------------- 
      if (iser.gt.0) fparm(iser,isite)=serval(ise)
!
!     ------------------------------------------------------
!     Set the matrix dimension if using a pruned basis set
!     ------------------------------------------------------
      if (basno(isite,ise).gt.0) then
         iparm(INDIM,isite) = ltbas( basno(isite,ise) ) 
         do i=0,6
            iparm(ILEMX+i,isite) = mts( i+1, basno(isite,ise) )
         end do
      end if
!
!     ------------------------------------------------------------
!     If B0 is zero, 
!       find B0 from the center field of the spectrum of the data,
!       resetting series value if necessary
!       If data are not available, report an error.
!     ------------------------------------------------------------
      if (sb0(ise).lt.RNDOFF) then
         if (ise.le.nspc) then
            sb0(ise)=sbi(ise)+sdb(ise)*(npts(ise)+1)/2.0D0
            fparm(IB0,isite)=sb0(ise)
            if (iser.eq.IB0) serval(ise)=sb0(ise)
!
            write (luout,1000) ise,ise,sb0(ise)
            if (luttyo.ne.luout) write (luttyo,1000)ise,ise,sb0(ise)
         else
            write (luout,1001) ise
            if (luttyo.ne.luout) write (luttyo,1000) ise
         end if
      end if
!
!     ------------------------------------------------------------
!     Reset negative GIB parameters to absolute values if needed
!     ------------------------------------------------------------
      if (fparm(IGIB0,isite).lt.0.0d0 .or. 
     #    ( (iser.eq.IPSI .or.iparm(INORT,isite).gt.1).and.
     #     fparm(IGIB2,isite).lt.0.0d0 ) ) then
         fparm(IGIB2,isite)=abs(fparm(IGIB0,isite)+fparm(IGIB2,isite))
     #                      -abs(fparm(IGIB0,isite))
         fparm(IGIB0,isite)=abs(fparm(IGIB0,isite))
         xreset=.true.
      end if


      return
!
 1000 format('*** B0(',i1,') outside range of spectrum ',i1,
     #': reset to',f8.1)
 1001 format('*** B0(',i1,') has not been specified ***')
!
      end
