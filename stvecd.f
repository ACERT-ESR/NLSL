c  VERSION 1.0  (NLSPMC version)   2/5/99
c*********************************************************************
c
c                       ===================
c                       SUBROUTINE : STVECD
c                       ===================
c
c     This subroutine calculates the starting vector for the diagonal
c     space, which is used to obtain the weighting factors of diagonal
c     eigenvectors.  It is intended for use with nonlinear lease-squares
c     applications.
c
c     On Entry :
c        evalx : eigenvalues of the off-diagonal space
c        evevx : eigenvectors of the off-diagonal space
c
c        Other parameters are passed through common block /eprprm/.
c
c     On Exit  :
c        stvd  : starting vector for diagonal space returned through
c                the common block /stvcom/
c
c     10-APR-93 Sanghyuk Lee
c
c     Includes :
c        nlsdim.inc
c        eprprm.inc
c        indexf.inc
c        stvcom.inc
c        wkspcm.inc
c
c*********************************************************************
c
      subroutine stvecd(evalx,evecx)
c
      implicit none
c
      include 'limits.inc'
      include 'simparm.inc'
      include 'basis.inc'
      include 'stvcom.inc'
      include 'wkspcm.inc'
c
      integer i,j,k,nt1
      double precision t1,cfact
c
      complex*16 czero
      parameter (czero=(0.0D0,0.0D0))
c
      complex*16 evalx,evecx,cnorm,znormu
      dimension evalx(mxegv),evecx(mxdim,mxegv)
c
      external znormu
c
c#####################################################################
c
      cfact=1.760844783D1
c
      do 10 i=1,nevo
         cwsp1(i)=czero
         do 20 j=1,ndimo
 20         cwsp1(i)=cwsp1(i)+evecx(j,i)*stvo(j)
 10   continue
c
      do 40 nt1=1,5
         t1=(init1+(nt1-1)*stept1)*cfact*1.0D-3
c
         do 50 i=1,nevo
 50         cwsp2(i)=cdexp(-evalx(i)*t1)
c
         do 60 i=1,ndimo
            stvd(i,nt1)=czero
            do 70 j=1,nevo
 70            stvd(i,nt1)=stvd(i,nt1)+evecx(i,j)*cwsp2(j)*cwsp1(j)
 60      continue
c
 40   continue
      do 100 nt1=1,5
         j=ndimd+1
         do 110 i=ndimo,1,-1
         if (pid(i).ne.1.and.pid(i).ne.2) then
           write(*,*)'pid error in stvecd ',i,pid(i)
           stop
         end if
         do 110 k=1,pid(i)
            j=j-1
            if(j.gt.MXDIM)then
		write(*,*)'in stvecd, MXDIM dimension exceeded ',MXDIM,j
                stop
            end if
            stvd(j,nt1)=pp(j)*stvd(i,nt1)
 110     continue
c
         if (nt1.eq.5) then
            do 120 i=1,ndimd
               stvd(i,1)=0.3d0*stvd(i,1)+0.2d0*stvd(i,2)+
     #           0.2d0*stvd(i,3)+0.2d0*stvd(i,4)+0.1d0*stvd(i,5)
 120        continue
         end if
c
 100  continue
c
      cnorm=znormu(stvd(1,1),ndimd)
      do 140 i=1,ndimd
 140     stvd(i,1)=stvd(i,1)/cnorm
c
      return
      end

