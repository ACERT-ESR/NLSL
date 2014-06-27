c NLSL Version 1.3.2 2/22/94
c----------------------------------------------------------------------
c                    =========================
c                        subroutine TDSQZ
c                    =========================
c
c   Loops through blocks stored in common /tridag/ and moves tridiagonal
c   matrices that may be re-used to the bottom of the storage space
c
c----------------------------------------------------------------------
c
      subroutine tdsqz()
      implicit none
c
      include 'nlsdim.inc'
      include 'expdat.inc'
      include 'tridag.inc'
      include 'stdio.inc'
c
      integer i,j,next,isp,isi,ixt,k,n
c
      next=1
      n=0
      do i=1,ntd
         isi=tdsite(i)
         isp=tdspec(i)
c
c        Check whether next tridiagonal matrix may be kept
c
         if (modtd(isi,isp).eq.0) then
c
c           Move the block to lowest available storage location
c
            ixt=ixtd(isi,isp)
            if (ixt.gt.next) then
               ixtd(isi,isp)=next
               do j=1,ltd(isi,isp)
                  k=ixt+j-1
                  alpha(next)=alpha(k)
                  beta(next)=beta(k)
                  next=next+1
               end do
            else if (ixt.eq.next) then
               next=next+ltd(isi,isp)
            else 
               write (luout,1005)
            end if
            n=n+1
            tdsite(n)=isi
            tdspec(n)=isp
         end if
c
      end do
      nexttd=next
      ntd=n
c
 1005 format('*** Error in tridiagonal matrix storage ***')
c 
      return
      end










