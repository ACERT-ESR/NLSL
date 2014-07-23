c  VERSION 1.0  (NLSPMC version)   2/5/99
c----------------------------------------------------------------------
c                    =========================
c                       subroutine WRFIT
c                    =========================
c
c Writes the NLS fit 2D spectrum into the file.
c
c----------------------------------------------------------------------
c      subroutine wrfit( spec,npt1,npt2,inif2,res2,inif1,res1,
c     #                  amin,amax,name,lth )
      subroutine wrfit( spec,icnt,amin,amax,name )
      implicit none
c
      include 'limits.inc'
      include 'stdio.inc'
      include 'datas.inc'
c
      integer i,j,npt1,npt2,icnt
      double precision inif2,res2,inif1,res1,amin,amax
cc      double precision spec(snpt1(icnt),snpt2(icnt))
      double precision spec(128,128) 
c
      character*(WORDLG) name
      character*30 xaxisf,yaxisf
      data xaxisf/'f2 (MHz)'/
      data yaxisf/'f1 (MHz)'/
c
c######################################################################
c
c the following is moved here from fitp.f
c
      npt1=snpt1(icnt)
      npt2=snpt2(icnt)
c calc amax if not done already
      if(amax.lt.1.0d-10) then
        do 10 j=1,npt2
          do 10 i=1,npt1
            if(spec(i,j).gt.amax) amax=spec(i,j)
 10     continue
      amin=0.0d0
      end if
c
      inif1=-5.0D2/sstept1(icnt)
      inif2=-5.0D2/sstept2(icnt)
      res1=10.0D2/snpt1(icnt)/sstept1(icnt)
      res2=10.0D2/snpt2(icnt)/sstept2(icnt)
c
      open(unit=ludisk,file=name,status='unknown',
     #     access='sequential',form='formatted')
      write (ludisk,3100) 0
      write (ludisk,3110) npt2,1,npt2,xaxisf
      write (ludisk,3110) npt1,1,npt1,yaxisf
      write (ludisk,3120) inif2,res2,inif1,res1,amin,amax
      write (ludisk,3120) ((spec(i,j),j=1,npt2),i=1,npt1)
      close (unit=ludisk)
c
      if(npt2.eq.1) then
        open(unit=ludisk,file='cwspec',status='unknown',
     #     access='sequential',form='formatted')
        write (ludisk,3130) (i,spec(i,1),
     #		(spec(i+1,1)-spec(i-1,1))/2.0d0,i=2,npt1-1)
        close (unit=ludisk)
      end if
      return
c
c #### format statements ########################################
c
 3100 format(i10)
 3110 format(3i10,5x,a30)
 3120 format(8e10.4)
 3130 format(1i10,2x,e10.4,2x,e10.4,/)
c
      end
