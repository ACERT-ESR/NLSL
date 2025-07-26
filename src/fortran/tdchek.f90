! Version 1.4 10/10/94
!----------------------------------------------------------------------
!                    =========================
!                        subroutine TDCHEK
!                    =========================
!
!  This routine checks whether a tridiagonal matrix is available for
!  a given site and spectrum and allocates one if necessary.
!
!----------------------------------------------------------------------
      subroutine tdchek(isite,ispc,ierr)
!
      use nlsdim
      use eprprm
      use parcom
      use tridag
      use basis
      use mtsdef
      use errmsg
!
      implicit none
      integer isite,ispc,ierr
!
      integer idummy,lthtd,mtstmp(MXMTS)
      character*30 ndummy
!
      integer setmts
      external setmts
!
!------------------------------------------------------------------------
!
      if (ispc.le.0 .or. ispc.gt.MXSPC
     #   .or. isite.le.0 .or. isite.gt.MXSITE) return
!
!   -------------------------------------------------------------------
!   Check whether new tridiagonal matrix is needed for the calculation
!   -------------------------------------------------------------------
!
!     ------------------------------------------------------------
!     SPECIAL CHECK: Avoid repeating a MOMD calculation.
!     If a MOMD calculation is being performed for later spectra
!     in a PSI series, check whether it is possible to use the 
!     tridiagonal matrix calculated from the first spectrum in the 
!     series (i.e. make sure the first matrix is unmodified).
!     ------------------------------------------------------------
      if (iparm(INORT,isite).gt.1 .and. iser.eq.IPSI
     #    .and. ispc.gt.1 .and. modtd(isite,1).eq.0)
     #    then
         ixtd(isite,ispc)=ixtd(isite,1)
         ltd(isite,ispc)=ltd(isite,1)
          modtd(isite,ispc)=0
      end if
!
      if (modtd(isite,ispc).ne.0) then
!
!        ------------------------------------------------ 
!        if nstep has been set, lthtd <- nstep
!        else if pruned basis set in use lthtd <- ltbas
!        else if ndim has been set lthtd <- ndim
!        else determine ndim and set lthtd, nstep <- ndim
!        ------------------------------------------------ 
!
         lthtd=iparm(INSTEP,isite)
         if (lthtd.le.0) then
            if (basno(isite,ispc).eq.0) then
               lthtd=iparm(INDIM,isite)
!
!              set NDIM if necessary
!
               if (lthtd.le.0) then
                  ierr=setmts(fparm(1,isite),iparm(1,isite),mtstmp)
                  ndummy=' '
                  idummy=0
                  call lbasix(ndummy,idummy,mtstmp,lthtd,0,1,ierr)
                  iparm(INDIM,isite)=lthtd
               end if
            else
               lthtd=ltbas(basno(isite,ispc))
            end if
         end if
!
         if (iparm(INORT,isite).gt.1) lthtd=lthtd*iparm(INORT,isite)
!
!       -------------------------------------------------------
!       Check that there is room for the new tridiagonal matrix
!       -------------------------------------------------------
!
         if (nexttd+lthtd.gt.MXTDG) then
            ierr=TDGBIG
            return
         end if
!
!       ----------------------------------
!       Allocate a new tridiagonal matrix
!       ----------------------------------
!
         modtd(isite,ispc)=1
         ltd(isite,ispc)=lthtd
          ixtd(isite,ispc)=nexttd
         nexttd=nexttd+lthtd
         ntd=ntd+1
         tdsite(ntd)=isite
         tdspec(ntd)=ispc
      end if
      return
      end
