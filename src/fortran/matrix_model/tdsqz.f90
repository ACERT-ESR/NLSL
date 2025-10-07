! NLSL Version 1.3.2 2/22/94
!----------------------------------------------------------------------
!                    =========================
!                        subroutine TDSQZ
!                    =========================
!
!   Loops through blocks stored in common /tridag/ and moves tridiagonal
!   matrices that may be re-used to the bottom of the storage space
!
!----------------------------------------------------------------------
!
      subroutine tdsqz()
!
      use nlsdim
      use expdat
      use tridag
      use stdio
!
      implicit none
      integer i,j,next,isp,isi,ixt,k,n
!
      next=1
      n=0
      do i=1,ntd
         isi=tdsite(i)
         isp=tdspec(i)
!
!        Check whether next tridiagonal matrix may be kept
!
         if (modtd(isi,isp).eq.0) then
!
!           Move the block to lowest available storage location
!
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
!
      end do
      nexttd=next
      ntd=n
!
 1005 format('*** Error in tridiagonal matrix storage ***')
! 
      return
      end










