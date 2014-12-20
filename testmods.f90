program testmods

use eprprm
implicit none

! Build instructions for the initial version of this program:
! gfortran -g -ffixed-form -c nlsdim.f90
! gfortran -g -ffixed-form -c parcom.f90
! gfortran -g -ffixed-form -c eprprm.f90
! gfortran -g testmods.f90 nlsdim.o parcom.o eprprm.o -o testmods

integer :: i, j, itest
double precision :: ftest

print *, 'Initializing fparm and iparm'
do j = 1, MXSITE
   do i = 1, NFPRM
      fparm(i,j) = (j-1) * NFPRM + i
   end do
   do i = 1, NIPRM
      iparm(i,j) = (j-1) * NIPRM + i
   end do
end do

print *, 'Starting test of modules'
do j = 1, MXSITE
   call select_site(j)
   itest = ndim
   ftest = pmzz
   print ("('Fetched by pointers for site', i3, ':', i4, f6.1)"), &
         j, itest, ftest
   print ("('The expected values for site', i3, ':', i4, f6.1)"), &
         j, iparm(indim,j), fparm(ipmzz,j)
   call arbsub(mmn,dx,j)
   print ("('From arrays in main for site', i3, ':', i4, f6.1)"), &
         j, iparm(immn,j), fparm(idx,j)
end do

contains

   subroutine arbsub(iarb,farb,jjj)
      integer :: iarb, jjj
      double precision :: farb
      print ("('Via passed pointers for site', i3, ':', i4, f6.1)"), &
            jjj, iarb, farb
   end subroutine arbsub

end program testmods
