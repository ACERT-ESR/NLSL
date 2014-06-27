c Version 1.4 10/10/94
c----------------------------------------------------------------------
c                    =========================
c                        subroutine TDCHEK
c                    =========================
c
c  This routine checks whether a tridiagonal matrix is available for
c  a given site and spectrum and allocates one if necessary.
c
c----------------------------------------------------------------------
      subroutine tdchek(isite,ispc,ierr)
      implicit none
      integer isite,ispc,ierr
c
      include 'nlsdim.inc'
      include 'eprprm.inc'
      include 'parcom.inc'
      include 'tridag.inc'
      include 'basis.inc'
      include 'mtsdef.inc'
      include 'errmsg.inc'
c
      integer idummy,lthtd,mtstmp(MXMTS)
      character*30 ndummy
c
      integer setmts
      external setmts
c
c------------------------------------------------------------------------
c
      if (ispc.le.0 .or. ispc.gt.MXSPC
     #   .or. isite.le.0 .or. isite.gt.MXSITE) return
c
c   -------------------------------------------------------------------
c   Check whether new tridiagonal matrix is needed for the calculation
c   -------------------------------------------------------------------
c
c     ------------------------------------------------------------
c     SPECIAL CHECK: Avoid repeating a MOMD calculation.
c     If a MOMD calculation is being performed for later spectra
c     in a PSI series, check whether it is possible to use the 
c     tridiagonal matrix calculated from the first spectrum in the 
c     series (i.e. make sure the first matrix is unmodified).
c     ------------------------------------------------------------
      if (iparm(INORT,isite).gt.1 .and. iser.eq.IPSI
     #    .and. ispc.gt.1 .and. modtd(isite,1).eq.0)
     #    then
         ixtd(isite,ispc)=ixtd(isite,1)
         ltd(isite,ispc)=ltd(isite,1)
          modtd(isite,ispc)=0
      end if
c
      if (modtd(isite,ispc).ne.0) then
c
c        ------------------------------------------------ 
c        if nstep has been set, lthtd <- nstep
c        else if pruned basis set in use lthtd <- ltbas
c        else if ndim has been set lthtd <- ndim
c        else determine ndim and set lthtd, nstep <- ndim
c        ------------------------------------------------ 
c
         lthtd=iparm(INSTEP,isite)
         if (lthtd.le.0) then
            if (basno(isite,ispc).eq.0) then
               lthtd=iparm(INDIM,isite)
c
c              set NDIM if necessary
c
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
c
         if (iparm(INORT,isite).gt.1) lthtd=lthtd*iparm(INORT,isite)
c
c       -------------------------------------------------------
c       Check that there is room for the new tridiagonal matrix
c       -------------------------------------------------------
c
         if (nexttd+lthtd.gt.MXTDG) then
            ierr=TDGBIG
            return
         end if
c
c       ----------------------------------
c       Allocate a new tridiagonal matrix
c       ----------------------------------
c
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
