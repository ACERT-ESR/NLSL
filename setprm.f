c NLSL Version 1.5 beta 11/23/95
c----------------------------------------------------------------------
c                    =========================
c                      subroutine SETPRM
c                    =========================
c>@brief  This file contains two routines that set a given parameter,
c>  specified by an index into the fparm or iparm array, given the
c>  parameter value and a site/spectrum index. The secondary index
c>  is interpreted according to whether or not the parameter belongs to
c>  a set of parameters that must be the same for all components of a
c>  given spectrum. This set of "spectral parameters" includes
c>
c>           PHASE   (spectral phase angle)
c>           PSI     (director tilt angle)
c>           B0      (spectrometer field)
c>           LB      (spectral line broadening)
c>           RANGE   (spectral range)
c>
c>  and is identified by the logical function spcpar(index) coded 
c>  in this file. These parameters are kept arrays separate from the
c>  fparm arrays in which the site parameters are kept.
c>
c>  The index ixsite is interpreted as referring to a spectrum for
c>  spectral parameters; otherwise, the index refers to a component
c>  within a given spectrum. 
c>  
c>  Whenever a parameter is changed, if recalculation of a spectrum
c>  with the new parameter value would require recalculation of the
c>  tridiagonal matrix, a flag is set for the given sites/spectra
c>  affected by that parameter.
c
c----------------------------------------------------------------------
      subroutine setprm(ixparm,ixsite,fval)
      implicit none
      integer ixparm,ixsite
      double precision fval
c
      integer i,ixmax,jx,jx1,jx2
c
      include 'nlsdim.inc'
      include 'eprprm.inc'
      include 'expdat.inc'
      include 'tridag.inc'
      include 'parcom.inc'
      include 'symdef.inc'
      include 'stdio.inc'
c
      logical spcpar,spcprm
      external spcpar
c   
      spcprm=spcpar(ixparm)
      if (spcprm) then
         ixmax=MXSPC
      else
         ixmax=MXSITE
      endif
c
c     --- Index between 0 and maximum: single site or spectrum
c
      if (ixsite.gt.0 .and. ixsite.le.ixmax) then
         jx1=ixsite
         jx2=ixsite
c
c     --- Zero index: all sites or spectra
c 
      else if (ixsite.eq.0) then
         jx1=1
         jx2=ixmax
c
c     --- Illegal index
c
      else
         write(luout,1000) jx
         if (luout.ne.luttyo) write (luttyo,1000) jx
      end if
c
      do jx=jx1,jx2
c
c       --- Spectral parameters are stored separately:
c
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
c
c           Mark tridiagonal matrix for recalculation only if tilt angle
c           has been changed (spectra can be calculated from an existing
c           tridiagonal matrix if the other parameters above are changed)
c
            if (ixparm.eq.IPSI) then
               do i=1,nsite
                  modtd(i,jx)=1
               end do
            end if
c            
c       --- Site parameter: mark tridiagonal matrix for recalculation
c           (except for Gaussian inhomogeneous broadening 
c            and isotropic Lorentzian broadening parameters) 
c   
         else
            fparm(abs(mod(ixparm,100)),jx)=fval
            if (ixparm.ne.IGIB0 .and. ixparm.ne.IGIB2 .and. .not.
     #         (ixparm.eq.IWXX .and.iparm(IIWFLG,ixsite).eq.SPHERICAL) )
     #       then
               do i=1,nser
                  modtd(jx,i)=1
               end do
            end if
c
         end if
c
      end do
c
      return
 1000 format(' *** Illegal index: ',i2,' ***')
      end


c----------------------------------------------------------------------
c                    =========================
c                      subroutine SETIPR
c                    =========================
c
c>@brief Analogous routine to setprm for integer parameters
c> There are only two user-settable integer spectrum parameters: nfield,
c> and ideriv. These are normally determined by the input data file, but
c> are needed when calculations are to performed without data or fitting.
c
c----------------------------------------------------------------------
      subroutine setipr(ixparm,ixsite,ival)
      implicit none
      integer ixparm,ixsite,ival
c
      integer i,ixmax,ixp,j,jx,jx1,jx2
c
      integer itrim
      logical spcpar,spcprm
      external itrim,spcpar
c
      include 'nlsdim.inc'
      include 'eprprm.inc'
      include 'expdat.inc'
      include 'parcom.inc'
      include 'tridag.inc'
      include 'basis.inc'
      include 'stdio.inc'
c
      ixp=abs(mod(ixparm,100))
      spcprm=spcpar(ixparm)
c
      if (spcprm) then
         ixmax=MXSPC
      else
         ixmax=MXSITE
      endif
c
c     Single site/spectrum
c
      if (ixsite.gt.0 .and. ixsite.le.MXSITE) then
         jx1=ixsite
         jx2=ixsite
c
c     All sites/spectra
c
      else if (ixsite.eq.0) then
         jx1=1
         jx2=ixmax
c
c     Illegal index
c
      else
         write(luout,1000) jx
         if (luout.ne.luttyo) write (luttyo,1000) jx
      end if
c
      do jx=jx1,jx2
c
c        Spectrum parameters
c
         if (spcprm) then
            if (ixp.eq.INFLD) then
c
c              Cannot change number of points in an existing datafile
c
               if (jx.le.nspc) then
                  write(luout,1001) dataid(jx)(:itrim(dataid(jx)))
                  if (luout.ne.luttyo) write (luttyo,1001)
     #                                 dataid(jx)(:itrim(dataid(jx)))
               else
                  npts(jx)=ival
                  if (ival.gt.1) sdb(jx)=srng(jx)/float(ival-1)
c
                  nft(jx)=1
 9                nft(jx)=nft(jx)*2
                  if (nft(jx).lt.npts(jx)) go to 9
               end if
c
            else if (ixparm.eq.IIDERV) then
               idrv(jx)=ival
            end if
c
c        Site parameters
c        Note that changes in site-related integer quantities always
c        require recalculation of the tridiagonal matrix
c
         else
            iparm(ixp,jx)=ival
            do i=1,MXSPC
               modtd(jx,i)=1
            end do
         end if
      end do
c
      return
 1000 format('*** Illegal index: ',i2,' ***')
 1001 format('*** Number of data points for file ',a,
     #       ' cannot be changed ***')
      end



c----------------------------------------------------------------------
c                    =========================
c                         function SPCPAR
c                    =========================
c
c>@brief Returns .true. if the index argument corresponds to a floating
c> point or integer parameter that cannot change for an individual spectrum,
c> regardless of the number of sites for which the calculation is made.
c> These include:
c>
c>       PHASE   (spectral phase angle)
c>       PSI     (director tilt angle)
c>       LB      (spectral line broadening)
c>       B0      (spectrometer field)
c>       FIELDI  (initial field of spectrum)
c>       DFLD    (field step per point in spectrum)
c>       RANGE   (field range of spectrum)
c>
c>       NFIELD  (number of points in spectrum)
c>       IDERIV  (0th/1st derivative flag for spectrum)
c
c
c      Includes:
c       nlsdim.inc
c       eprprm.inc
c       parcom.inc
c
c----------------------------------------------------------------------
      function spcpar( ix )
      implicit none
      integer ityp,ix,ixp
      logical spcpar
c
      include 'nlsdim.inc'
      include 'eprprm.inc'
      include 'parcom.inc'
c
      ixp = mod(ix,100)
      ityp = ix/100
      spcpar = (ityp.eq.0 .and.
     #            (ixp.eq.IPHASE .or. ixp.eq.IPSI .or. ixp.eq.ILB
     #            .or. ixp.eq.IB0 .or. ixp.eq.IFLDI .or. ixp.eq.IDFLD
     #            .or. ixp.eq.IRANGE) )
     #     .or.(ityp.eq.1 .and.
     #            (ixp.eq.INFLD .or. ixp.eq.IIDERV) )
c
      return
      end

      function ptype( ix )
      implicit none
      integer ix
      character*7 ptype
c
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
      implicit none
      integer ixlim, ix
c
      include 'nlsdim.inc'
      include 'expdat.inc'
      include 'parcom.inc'
      logical spcpar
      external spcpar
c
      if (spcpar(ix)) then
         ixlim=nser
      else
         ixlim=nsite
      end if
      return
      end

c
c------------------------------------------------------------------------
c                    =========================
c                         function GETPRM
c                    =========================
c>@brief Given a parameter index and a site/spectrum index, this function
c> returns the value of the parameter from the fparm array (for
c> site parameters) or from the spectral parameter arrays for spectrum
c> parameters. 
c------------------------------------------------------------------------
c

      function getprm(ixparm,ixsite)
      implicit none
      integer ixparm,ixsite
      double precision getprm
c
      include 'nlsdim.inc'
      include 'eprprm.inc'
      include 'expdat.inc'
      include 'tridag.inc'
      include 'parcom.inc'
      include 'symdef.inc'
      include 'stdio.inc'
c
      integer ixmax
c
      logical spcpar,spcprm
      external spcpar
c   
      spcprm=spcpar(ixparm)
      if (spcprm) then
         ixmax=nser
      else
         ixmax=nsite
      endif
c
c     --- Illegal index
c
      if (ixsite.le.0.or.ixsite.gt.ixmax) then
         write(luttyo,1000) ixsite
         return
      end if
c
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
c     
      else
         getprm=fparm(abs(mod(ixparm,100)),ixsite)
      end if
c
      return
c
 1000 format(' *** Illegal index: ',i2,' ***')
      end
