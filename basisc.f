c NLSL Version 1.5 beta 11/23/95
c----------------------------------------------------------------------
c                    =========================
c                       subroutine BASISC
c                    =========================
c
c  Reads in the basis index contained in the given file name
c
c basis <name> { spectrum <specindx> site <siteindx> } 
c
c      name    : name of the basis index file (without .ind extension)
c      specindx: optional index of the spectrum to which the basis set
c                is to be assigned (may be a dataID)
c      siteindx: optional index of the spectrum to which the basis set
c                is to be assigned (may be a dataID)
c
c----------------------------------------------------------------------
      subroutine basisc(line)
      implicit none
      character*80 line
c
      include 'nlsdim.inc'
      include 'eprprm.inc'
      include 'parcom.inc'
      include 'basis.inc'
      include 'tridag.inc'
      include 'stdio.inc'
c
      integer i,ierr,isi,isp,ival,ixsp1,ixsp2,ixsi1,ixsi2,ixss(2),
     #        lth,maxb,new
      character token*30,ixname*30
c
      integer iroot
      logical itoken
      external iroot,itoken
c
c----------------------------------------------------------------------
c  Get the name of the basis set file
c----------------------------------------------------------------------
      call gettkn(line,token,lth)
c
c                           *** No filename specified
      if (lth.eq.0) then
         write(luttyo,1000)
         if (lth.eq.0) return
      end if
c
c                            *** Max # basis sets exceeded
      if (nbas.ge.MXTDM) then
         write(luout,1004) MXTDM
         if (luout.ne.luttyo) write (luttyo,1004) MXTDM
         return
      end if
      nbas=nbas+1
      basisid(nbas)=token
      ixname=token
c
      call getss(line,ixss)
c
c----------------------------------------------------------------------
c     Get the basis set
c----------------------------------------------------------------------
c 
      maxb=MXDIM-nextbs+1
      new=0
      call lbasix(ixname,ibasis(1,nextbs),mts(1,nbas),
     #            ltbas(nbas),maxb,new,ierr)
c
      if (ierr.eq.0) then
         ixbas(nbas)=nextbs
         nextbs=nextbs+ltbas(nbas)
c
         lth = iroot(ixname)
         write (luout,1001) ixname(:lth)//'.ind',(mts(i,nbas),i=1,7)
         if (luout.ne.luttyo) 
     *      write (luttyo,1001) ixname(:lth)//'.ind',(mts(i,nbas),i=1,7)
c
         if (ixss(1).le.MXSPC .and. ixss(2).le.MXSITE) then
c
c          -----------------------------------
c          Set ranges of spectra, site indices
c          -----------------------------------
            if (ixss(2).le.0) then
               ixsi1=1
               ixsi2=MXSITE
            else
               ixsi1=ixss(2)
               ixsi2=ixss(2)
            end if
c
            if (ixss(1).le.0) then
               ixsp1=1
               ixsp2=MXSPC
            else
               ixsp1=ixss(1)
               ixsp2=ixss(1)
            end if
c
c          -----------------------
c           Now assign basis set 
c          ----------------------
            do isi=ixsi1,ixsi2 
               do isp=ixsp1,ixsp2
                  basno(isi,isp)=nbas
                  modtd(isi,isp)=1
               end do
            end do
c
         end if
c
      else
         nbas=nbas-1
         write (luout,1006) ixname(:lth)
         if (luout.ne.luttyo) write (luttyo,1006) ixname(:lth)
      end if
c
      return
c
c ###### format statements ############################################
c
 1000 format('*** No basis set file ID specifed ***')
 1001 format('*** File ',a,' read ***'/' Lemx=',i3,' Lomx=',i3,
     *       ' K{mn,mx}= {',i3,',',i3,'} M{mn,mx}= {',i3,',',i3,
     *       '} ipnmx =',i2)
 1004 format('*** Maximum of',i2,' basis sets exceeded ***')
 1006 format('*** Error reading ''',a,''': no basis set added ***')
      end

c----------------------------------------------------------------------
c                    =========================
c                       subroutine DELETC
c                    =========================
c
c  Deletes the specified basis set from the buffer
c
c  Syntax:
c   delete <name>|index
c
c      name    : I.D. of the basis set
c      index   : Index number of the basis set
c
c----------------------------------------------------------------------
      subroutine deletc(line)
      implicit none
      character*80 line
c
      include 'nlsdim.inc'
      include 'eprprm.inc'
      include 'parcom.inc'
      include 'basis.inc'
      include 'tridag.inc'
      include 'stdio.inc'
c
      integer istrt,isi,isp,ixsm,ival,lth,i,j,k,m
      character token*30,ixname*30
c
      logical itoken
      integer isfind
      external itoken,isfind
c
c----------------------------------------------------------------------
c  Get the name of the basis set file
c----------------------------------------------------------------------
      call gettkn(line,token,lth)
c
c                           *** No basis set specified
      if (lth.eq.0) then
         write(luttyo,1000)
         if (lth.eq.0) return
      end if
c
c     ----------------------------------------------------------
c      Token is not a keyword: check whether it is a symbolic or
c      integer value
c     -----------------------------------------------------------
      ixsm=isfind(token,lth)
      if (ixsm.eq.0) then
         if (.not.itoken(token,lth,ival)) then
c                                           *** Illegal index
            write(luout,1002) token(:lth)
            if (luout.ne.luttyo)write(luttyo,1002) token(:lth)
            return
         end if
      else
         ival=abs(ixsm)
      end if
c                                            *** Undefined index
      if (ival.le.0 .or. ival.gt.nbas) then
         write (luout,1002) token(:lth)
         if (luout.ne.luttyo) write (luttyo,1002) token(:lth)
         return
      end if
c
c     ------------------------------
c     Move basis set parameters down
c     ------------------------------
      istrt=ixbas(ival)
      lth=ltbas(ival)
      nbas=nbas-1
      do i=ival,nbas
         ixbas(i)=ixbas(i+1)
         ltbas(i)=ltbas(i+1)
         bsused(i)=bsused(i+1)
         basisID(i)=basisID(i+1)
         do m=1,9
            mts(m,i)=mts(m,i+1)
         end do
      end do
c
c     ------------------------------
c     Move basis set indices down
c     ------------------------------
      j=istrt
      k=istrt+lth
      do i=1,lth
         do m=1,5
            ibasis(m,j)=ibasis(m,k)
         end do
         j=j+1
         k=k+1
      end do
      nextbs=nextbs-lth
c
c    ------------------------------------
c     Now deassign the deleted basis set
c    ------------------------------------
      do isi=1,MXSITE
         do isp=1,MXSPC
            if (basno(isi,isp).eq.ival) then
               basno(isi,isp)=0
               modtd(isi,isp)=1
            end if
         end do
      end do
c
      return
c
c ###### format statements ############################################
c
 1000 format('*** No basis set specifed ***')
 1002 format('*** Basis set ''',a,''' not defined ***') 
      end

