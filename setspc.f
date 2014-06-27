c NLSL Version 1.5.1 beta 2/3/96
c----------------------------------------------------------------------
c                    =========================
c                       subroutine SETSPC
c                    =========================
c
c  This routine sets all the parameters in the fparm and iparm array
c  (kept in common /parcom/) for the specified site and spectrum
c  (arguments isite and ispc).
c
c  Certain parameter values are not kept in the fparm and iparm
c  arrays. In particular, if a series of spectra are to be calculated,
c  the series variable appropriate for the current spectrum is loaded
c  into fparm. This routine also sets the spectrum-dependent parameters
c  from the arrays in common /expdat/.  These parameters are not allowed
c  to vary within a single spectrum, and specifically include the following:
c
c    Floating:
c       B0
c       LB
c       PHASE
c       PSI
c       FLDI
c       DFLD
c
c    Integer:
c       NFLD
c       IDERIV
c----------------------------------------------------------------------
      subroutine setspc(isite,ise)
      implicit none
      integer isite,ise
c
      include 'nlsdim.inc'
      include 'eprprm.inc'
      include 'expdat.inc'
      include 'parcom.inc'
      include 'basis.inc'
      include 'errmsg.inc'
      include 'iterat.inc'
      include 'rndoff.inc'
      include 'stdio.inc'
      include 'pidef.inc'
c
c
      integer i
c
c######################################################################
c
      if (ise.le.0 .or. ise.gt.MXSPC
     #   .or. isite.le.0 .or. isite.gt.MXSITE) return
c
c     --------------------------
c     Set spectral parameters
c     --------------------------
c
      fparm(ILB,   isite) = slb(ise)
      fparm(IPHASE,isite) = sphs(ise)
      fparm(IPSI,  isite) = spsi(ise)
      fparm(IB0,   isite) = sb0(ise)
      fparm(IFLDI, isite) = sbi(ise)-shft(ise)-tmpshft(ise)
      fparm(IDFLD, isite) = sdb(ise)
      iparm(INFLD, isite) = npts(ise)
      iparm(IIDERV,isite) = idrv(ise)
c
c     -------------------------------------- 
c     set the SERIES variable if necessary
c     -------------------------------------- 
      if (iser.gt.0) fparm(iser,isite)=serval(ise)
c
c     ------------------------------------------------------
c     Set the matrix dimension if using a pruned basis set
c     ------------------------------------------------------
      if (basno(isite,ise).gt.0) then
         iparm(INDIM,isite) = ltbas( basno(isite,ise) ) 
         do i=0,6
            iparm(ILEMX+i,isite) = mts( i+1, basno(isite,ise) )
         end do
      end if
c
c     ------------------------------------------------------------
c     If B0 is zero, 
c       find B0 from the center field of the spectrum of the data,
c       resetting series value if necessary
c       If data are not available, report an error.
c     ------------------------------------------------------------
      if (sb0(ise).lt.RNDOFF) then
         if (ise.le.nspc) then
            sb0(ise)=sbi(ise)+sdb(ise)*(npts(ise)+1)/2.0D0
            fparm(IB0,isite)=sb0(ise)
            if (iser.eq.IB0) serval(ise)=sb0(ise)
c
            write (luout,1000) ise,ise,sb0(ise)
            if (luttyo.ne.luout) write (luttyo,1000)ise,ise,sb0(ise)
         else
            write (luout,1001) ise
            if (luttyo.ne.luout) write (luttyo,1000) ise
         end if
      end if
c
c     ------------------------------------------------------------
c     Reset negative GIB parameters to absolute values if needed
c     ------------------------------------------------------------
      if (fparm(IGIB0,isite).lt.0.0d0 .or. 
     #    ( (iser.eq.IPSI .or.iparm(INORT,isite).gt.1).and.
     #     fparm(IGIB2,isite).lt.0.0d0 ) ) then
         fparm(IGIB2,isite)=abs(fparm(IGIB0,isite)+fparm(IGIB2,isite))
     #                      -abs(fparm(IGIB0,isite))
         fparm(IGIB0,isite)=abs(fparm(IGIB0,isite))
         xreset=.true.
      end if


      return
c
 1000 format('*** B0(',i1,') outside range of spectrum ',i1,
     #': reset to',f8.1)
 1001 format('*** B0(',i1,') has not been specified ***')
c
      end
