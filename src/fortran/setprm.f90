! NLSL Version 1.5 beta 11/23/95
!----------------------------------------------------------------------
!                    =========================
!                      subroutine SETPRM
!                    =========================
!>@brief  This file contains two routines that set a given parameter,
!>  specified by an index into the fparm or iparm array, given the
!>  parameter value and a site/spectrum index. The secondary index
!>  is interpreted according to whether or not the parameter belongs to
!>  a set of parameters that must be the same for all components of a
!>  given spectrum. This set of "spectral parameters" includes
!>
!>           PHASE   (spectral phase angle)
!>           PSI     (director tilt angle)
!>           B0      (spectrometer field)
!>           LB      (spectral line broadening)
!>           RANGE   (spectral range)
!>
!>  and is identified by the logical function spcpar(index) coded 
!>  in this file. These parameters are kept arrays separate from the
!>  fparm arrays in which the site parameters are kept.
!>
!>  The index ixsite is interpreted as referring to a spectrum for
!>  spectral parameters; otherwise, the index refers to a component
!>  within a given spectrum. 
!>  
!>  Whenever a parameter is changed, if recalculation of a spectrum
!>  with the new parameter value would require recalculation of the
!>  tridiagonal matrix, a flag is set for the given sites/spectra
!>  affected by that parameter.
!
!----------------------------------------------------------------------
      subroutine setprm(ixparm,ixsite,fval)
!
      use nlsdim
      use eprprm
      use expdat
      use tridag
      use parcom
      use symdef
      use stdio
!
      implicit none
      integer ixparm,ixsite
      double precision fval
!
      integer i,ixmax,jx,jx1,jx2
!
      logical spcpar,spcprm
      external spcpar
!
      jx1=1
      jx2=1
      spcprm=spcpar(ixparm)
      if (spcprm) then
         ixmax=MXSPC
      else
         ixmax=MXSITE
      endif
!
!     --- Index between 0 and maximum: single site or spectrum
!
      if (ixsite.gt.0 .and. ixsite.le.ixmax) then
         jx1=ixsite
         jx2=ixsite
!
!     --- Zero index: all sites or spectra
! 
      else if (ixsite.eq.0) then
         jx1=1
         jx2=ixmax
!
!     --- Illegal index
!
      else
         write(luout,1000) jx
         if (luout.ne.luttyo) write (luttyo,1000) jx
      end if
!
      do jx=jx1,jx2
!
!       --- Spectral parameters are stored separately:
!
         if (spcprm) then
            if (ixparm.eq.IB0) then
               sb0(jx)=fval
            else if (ixparm.eq.IPHASE) then
               sphs(jx)=fval
            else if (ixparm.eq.IPSI) then
               spsi(jx)=fval
            else if (ixparm.eq.ILB) then
               slb(jx)=fval
            else if (ixparm.eq.IRANGE) then
               srng(jx)=fval
               sbi(jx)=sb0(jx)-0.5d0*fval
               if (npts(jx).gt.1) sdb(jx)=fval/float(npts(jx)-1)
            end if
!
!           Mark tridiagonal matrix for recalculation only if tilt angle
!           has been changed (spectra can be calculated from an existing
!           tridiagonal matrix if the other parameters above are changed)
!
            if (ixparm.eq.IPSI) then
               do i=1,nsite
                  modtd(i,jx)=1
               end do
            end if
!            
!       --- Site parameter: mark tridiagonal matrix for recalculation
!           (except for Gaussian inhomogeneous broadening 
!            and isotropic Lorentzian broadening parameters) 
!   
         else
            fparm(abs(mod(ixparm,100)),jx)=fval
            if (ixparm.ne.IGIB0 .and. ixparm.ne.IGIB2 .and. .not.
     #         (ixparm.eq.IWXX .and.iparm(IIWFLG,ixsite).eq.SPHERICAL) )
     #       then
               do i=1,nser
                  modtd(jx,i)=1
               end do
            end if
!
         end if
!
      end do
!
      return
 1000 format(' *** Illegal index: ',i2,' ***')
      end


!----------------------------------------------------------------------
!                    =========================
!                      subroutine SETIPR
!                    =========================
!
!>@brief Analogous routine to setprm for integer parameters
!> There are only two user-settable integer spectrum parameters: nfield,
!> and ideriv. These are normally determined by the input data file, but
!> are needed when calculations are to performed without data or fitting.
!
!----------------------------------------------------------------------
      subroutine setipr(ixparm,ixsite,ival)
!
      use nlsdim
      use eprprm
      use expdat
      use parcom
      use tridag
      use basis
      use stdio
!
      implicit none
      integer ixparm,ixsite,ival
!
      integer i,ixmax,ixp,jx,jx1,jx2
!
      integer itrim
      logical spcpar,spcprm
      external itrim,spcpar
!
      jx1=1
      jx2=1
      ixp=abs(mod(ixparm,100))
      spcprm=spcpar(ixparm)
!
      if (spcprm) then
         ixmax=MXSPC
      else
         ixmax=MXSITE
      endif
!
!     Single site/spectrum
!
      if (ixsite.gt.0 .and. ixsite.le.MXSITE) then
         jx1=ixsite
         jx2=ixsite
!
!     All sites/spectra
!
      else if (ixsite.eq.0) then
         jx1=1
         jx2=ixmax
!
!     Illegal index
!
      else
         write(luout,1000) jx
         if (luout.ne.luttyo) write (luttyo,1000) jx
      end if
!
      do jx=jx1,jx2
!
!        Spectrum parameters
!
         if (spcprm) then
            if (ixp.eq.INFLD) then
!
!              Cannot change number of points in an existing datafile
!
               if (jx.le.nspc) then
                  write(luout,1001) dataid(jx)(:itrim(dataid(jx)))
                  if (luout.ne.luttyo) write (luttyo,1001)
     #                                 dataid(jx)(:itrim(dataid(jx)))
               else
                  npts(jx)=ival
                  if (ival.gt.1) sdb(jx)=srng(jx)/float(ival-1)
                  iparm(INFLD,jx)=ival
!
                  nft(jx)=1
 9                nft(jx)=nft(jx)*2
                  if (nft(jx).lt.npts(jx)) go to 9
               end if
!
            else if (ixp.eq.IIDERV) then
               idrv(jx)=ival
               iparm(IIDERV,jx)=ival
            end if
!
!        Site parameters
!        Note that changes in site-related integer quantities always
!        require recalculation of the tridiagonal matrix
!
         else
            iparm(ixp,jx)=ival
            do i=1,MXSPC
               modtd(jx,i)=1
            end do
         end if
      end do
!
      return
 1000 format('*** Illegal index: ',i2,' ***')
 1001 format('*** Number of data points for file ',a,
     #       ' cannot be changed ***')
      end



!----------------------------------------------------------------------
!                    =========================
!                         function SPCPAR
!                    =========================
!
!>@brief Returns .true. if the index argument corresponds to a floating
!> point or integer parameter that cannot change for an individual spectrum,
!> regardless of the number of sites for which the calculation is made.
!> These include:
!>
!>       PHASE   (spectral phase angle)
!>       PSI     (director tilt angle)
!>       LB      (spectral line broadening)
!>       B0      (spectrometer field)
!>       FIELDI  (initial field of spectrum)
!>       DFLD    (field step per point in spectrum)
!>       RANGE   (field range of spectrum)
!>
!>       NFIELD  (number of points in spectrum)
!>       IDERIV  (0th/1st derivative flag for spectrum)
!
!
!      Includes:
!       nlsdim.inc
!       eprprm.inc
!       parcom.inc
!
!----------------------------------------------------------------------
      function spcpar( ix )
!
      use nlsdim
      use eprprm
      use parcom
!
      implicit none
      integer ityp,ix,ixp
      logical spcpar
!
      ixp = mod(ix,100)
      ityp = ix/100
      spcpar = (ityp.eq.0 .and.
     #            (ixp.eq.IPHASE .or. ixp.eq.IPSI .or. ixp.eq.ILB
     #            .or. ixp.eq.IB0 .or. ixp.eq.IFLDI .or. ixp.eq.IDFLD
     #            .or. ixp.eq.IRANGE) )
     #     .or.(ityp.eq.1 .and.
     #            (ixp.eq.INFLD .or. ixp.eq.IIDERV) )
!
      return
      end

      function ptype( ix )
      implicit none
      integer ix
      character*7 ptype
!
      logical spcpar
      external spcpar
      if (spcpar(ix)) then
         ptype='spectra'
      else
         ptype='sites'
      end if
      return
      end

      function ixlim( ix )
!
      use nlsdim
      use expdat
      use parcom
!
      implicit none
      integer ixlim, ix
      logical spcpar
      external spcpar
!
      if (spcpar(ix)) then
         ixlim=nser
      else
         ixlim=nsite
      end if
      return
      end

!
!------------------------------------------------------------------------
!                    =========================
!                         function GETPRM
!                    =========================
!>@brief Given a parameter index and a site/spectrum index, this function
!> returns the value of the parameter from the fparm array (for
!> site parameters) or from the spectral parameter arrays for spectrum
!> parameters. 
!------------------------------------------------------------------------
!

      function getprm(ixparm,ixsite)
!
      use nlsdim
      use eprprm
      use expdat
      use tridag
      use parcom
      use symdef
      use stdio
!
      implicit none
      integer ixparm,ixsite
      double precision getprm
!
      integer ixmax
!
      logical spcpar,spcprm
      external spcpar
!   
      getprm=0.0D0
      spcprm=spcpar(ixparm)
      if (spcprm) then
         ixmax=nser
      else
         ixmax=nsite
      endif
!
!     --- Illegal index
!
      if (ixsite.le.0.or.ixsite.gt.ixmax) then
         write(luttyo,1000) ixsite
         return
      end if
!
      if (spcprm) then
         if (ixparm.eq.IB0) then
            getprm=sb0(ixsite)
         else if (ixparm.eq.IPHASE) then
            getprm=sphs(ixsite)
         else if (ixparm.eq.IPSI) then
            getprm=spsi(ixsite)
         else if (ixparm.eq.ILB) then
            getprm=slb(ixsite)
         else if (ixparm.eq.IRANGE) then
            getprm=srng(ixsite)
         end if
!     
      else
         getprm=fparm(abs(mod(ixparm,100)),ixsite)
      end if
!
      return
!
 1000 format(' *** Illegal index: ',i2,' ***')
      end
