c     Version 1.5  5/2/94
c**********************************************************************
c
c          (C)ontinued (F)raction (V)alue and (D)erivative
c          ===============================================
c
c   This routine calculates a continued fraction of the form
c
c                 1
c   v(z) =   ----------------
c             a(1) - z - b(1)
c                       ----------------
c                        a(2) - z - b(2)
c                                  ----------------
c                                   a(3) - z - b(3)
c                                         .
c                                             .
c                                                .
c                                                  - b(n-1)
c                                                   ----------
c                                                    a(n) - z
c
c
c     Notes:
c       1) continued-fraction value OR its first derivative are returned
c          in val for all nz values of z in a single call. Choice of first 
c          derivative is specified by setting ideriv .ne. 0
c       2) z is assumed to be a REAL number (not the general case)
c       3) A complex-valued diagonal shift w may be added to the a vector
c
c     written 30-SEP-90 DJS
c     modified 5-MAY-92 DEB for use with EPRLL/EPRCGL and EPRLF/EPRCGF
c
c**********************************************************************
c
      subroutine cfvd(a,b,n,z0,dz,nz,w,val,ideriv)
c
      implicit none
c
      integer ideriv,n,nz
      real*8 w,z0,dz
      complex*16 a(n),b(n),val(nz)
c
      integer iz,k
      real*8 z
      complex*16 tv,td,s,x
c
      complex*16 UNITY,CI
      parameter (UNITY=(1.0D0,0.0D0),CI=(0.0D0,1.0D0))
c
c######################################################################
c
c----------------------------------------------------------------------
c    Compute 0th derivative of continued fraction with respect to z
c----------------------------------------------------------------------
      if (ideriv.eq.0) then
        z=z0-dz
        do iz=1,nz
           z=z+dz
           x=dcmplx(w,z)
           s=unity/(a(n)+x)
           tv=s*b(n-1)*b(n-1)
           td=-tv*s
           do k=n-1,2,-1
              s=unity/(a(k)+x-tv)
              tv=s*b(k-1)*b(k-1)
           end do
           s=unity/(a(1)+x-tv)
           tv=s
           val(iz)=tv
        end do
c
c----------------------------------------------------------------------
c    Compute 1st derivative of continued fraction with respect to z
c----------------------------------------------------------------------
      else
         z=z0-dz
         do iz=1,nz
            z=z+dz
            x=dcmplx(w,z)
            s=unity/(a(n)+x)
            tv=s*b(n-1)*b(n-1)
            td=-tv*s
            do k=n-1,2,-1
               s=unity/(a(k)+x-tv)
               tv=s*b(k-1)*b(k-1)
              td=-tv*(unity-td)*s
            end do
            s=unity/(a(1)+x-tv)
            tv=s
            val(iz)=ci*tv*(unity-td)*s
         end do
      end if
      return
      end
