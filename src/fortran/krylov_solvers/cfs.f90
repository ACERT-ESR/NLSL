! Version 1.3.2       5/1/92
!**********************************************************************
!
!          (C)ontinued (F)raction (S)pectral calculation
!          ===============================================
!
!   This routine calculates a continued fraction of the form
!
!                 1
!   v(z) =   ----------------
!             a(1) - z - b(1)
!                       ----------------
!                        a(2) - z - b(2)
!                                  ----------------
!                                   a(3) - z - b(3)
!                                         .
!                                             .
!                                                .
!                                                  - b(n-1)
!                                                   ----------
!                                                    a(n) - z
!
!
!     Notes:
!       1) continued-fraction value OR its first derivative are returned
!          in val for all nz values of z in a single call. Choice of first 
!          derivative is specified by setting ideriv .ne. 0
!       2) z is assumed to be a REAL number (not the general case)
!       3) A complex-valued diagonal shift w may be added to the a vector
!
!     written 30-SEP-90 DJS
!     modified 5-MAY-92 DEB for use with EPRLL/EPRCGL and EPRLF/EPRCGF
!
!     The routine returns a spectrum calculated from both the absorption 
!     and dispersion parts of the spectrum using the given phase angle
!     (in degrees)
!
!     Code modified from CFVD program of D.J. Schneider,
!
!**********************************************************************
!
      subroutine cfs(a,b,n,z0,dz,nz,w,ideriv,phs,val)
!
      implicit none
!
      integer ideriv,n,nz
      double precision w,z,z0,dz,phs,val(nz)
      double complex a(n),b(n),ephs
!
      integer iz,k,ictr
      double complex tv,td,s,x
!
      double complex unity,ci
      parameter (unity=(1.0D0,0.0D0),ci=(0.0D0,1.0D0))
!
      double precision RADIAN
      parameter (RADIAN=0.174532925199433D-01)
!
!######################################################################
!
      ephs=dcmplx( dcos(RADIAN*phs), dsin(RADIAN*phs) )
!
!    -----------------------------------------------------------------
!     Compute 0th derivative of continued fraction with respect to z
!    -----------------------------------------------------------------
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
            val(iz)=dreal(ephs*tv)
         end do
!
!    -----------------------------------------------------------------
!     Compute 1st derivative of continued fraction with respect to z
!    -----------------------------------------------------------------
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
            val(iz)=dreal( ephs*ci*tv*(unity-td)*s )
         end do
      end if
!
!    -----------------------------------------------------------
!     Special check to eliminate spikes near center of spectrum
!    -----------------------------------------------------------
      ictr=int(1.5d0-z0/dz)
      if(ictr .gt. 1 .and. ictr .lt. iz-1) then
        if ( abs(val(ictr)-val(ictr-1))+abs(val(ictr)-val(ictr+1))
     #    .gt. 3.0d0*abs(val(ictr+1)-val(ictr-1)) )
     #    val(ictr)=0.5D0*(val(ictr+1)+val(ictr-1) ) 
      end if
!
      return
      end
