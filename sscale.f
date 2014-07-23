c  VERSION 1.0  (NLSPMC version)   2/5/99
c----------------------------------------------------------------------
c                    =========================
c                      subroutine SSCALE
c                    =========================
c
c     For spectra in series, this subroutine scales the simulated
c  spectrum so that it matches the experimental spectrum.  The scale
c  factor may have vary for spectra with idepnd=0.  
c
c     In case more than two spectra should share the same scale
c  factor, the dependency flag for the first spectrum in that group
c  should be 0 and those for the rest of the spectra in that group
c  should be 1.  This situation happens when one want to fit Sc+ &
c  Sc- spectra from single experiment simutaneously or when one tries
c  to estimate We (T1 relaxation constant) from the intensities of
c  several ELDOR spectra obtained with the identical experimental
c  condition (e.g. amplifier gain, # of averages, same dead time,
c  etc.).  For this strategy to work, one needs to put those dependent
c  spectra in the contiguous space and only the first one should have
c  idepnd=0 in data command.
c  
c  Includes:
c     nlsdim.inc
c     expdat.inc
c 
c---------------------------------------------------------------------- 
c
      subroutine sscale(spct)
c
      implicit none
c
      include 'limits.inc'
      include 'datas.inc'
c
      double precision spct(mxpt)
c
      integer i,isp,ispi,ixs,mpts
      double precision asum,bsum
c
      double precision zero
      parameter (zero=0.0d0)
c
c######################################################################
c
      isp=1
      ispi=1
      ixs=ixsp(1)
      mpts=ndata(1)
c                   * search next independent spectra (idepnd=0)
 10   isp=isp+1
      if (idepnd(isp).eq.1) then	! if 1, share this scale factor
         mpts=mpts+ndata(isp)
         go to 10
      end if
c                   * obtain scale factor from group of data sets
      asum=zero
      bsum=zero
      do 20 i=ixs,ixs+mpts-1
         asum=asum+spct(i)*spct(i)
         bsum=bsum+data(i)*spct(i)
 20   continue
c
      do 30 i=ispi,isp-1
         sfac(i)=bsum/asum
 30   continue
c                              * scale the spectra
      do 40 i=ixs,ixs+mpts-1
         spct(i)=spct(i)*sfac(ispi)
 40   continue
c
      if (isp.gt.nspc) return
c
      ispi=isp
      ixs=ixsp(isp)
      mpts=ndata(isp)
      go to 10
c
      end
