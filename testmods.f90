      program testmods

      use eprprm
      use errmsg
      use lpnam
      implicit none

! Build instructions for the initial version of this program:
! gfortran -g -ffixed-form -c nlsdim.f90
! gfortran -g -ffixed-form -c parcom.f90
! gfortran -g -ffixed-form -c eprprm.f90
! gfortran -g testmods.f90 nlsdim.o parcom.o eprprm.o -o testmods
! To build for testing in Windows, using the Intel compiler:
! ifort /debug /fixed /c nlsdim.f90
! ifort /debug /fixed /c parcom.f90
! ifort /debug /fixed /c eprprm.f90
! ifort /debug /o testmods testmods.f90 nlsdim.obj parcom.obj eprprm.obj

      integer :: i, j, itest, iitest, iftest, info
      double precision :: ftest

      print *, 'Initializing fparm and iparm...'
      do j = 1, MXSITE
         do i = 1, NFPRM
            fparm(i,j) = (j-1) * NFPRM + i
         end do
         do i = 1, NIPRM
            iparm(i,j) = (j-1) * NIPRM + i
         end do
      end do

      print *
      print *, 'Starting test of eprprm pointers with changing sites...'
      print *, 'Correct values advance by 24 and 43.0 for each new site'
      do j = 1, MXSITE
         print ("('Spot testing values for site', i3 ,'...')"), j
         call select_site(j)
         iitest = immn
         iftest = ipmzz
         itest = mmn
         ftest = pmzz
         print ("('The expected values for site', i3,                    &
     &            ', indexes', i3, ' and', i3, ':', i5, f7.1)"),         &
     &            j, iitest, iftest, iparm(iitest,j), fparm(iftest,j)
         print ("('Fetched by pointers for site', i3,                    &
     &            ', indexes', i3, ' and', i3, ':', i5, f7.1)"),         &
     &            j, iitest, iftest, itest, ftest
         ! now try some tests from within a subroutine
         call arbsub(j,immn,ipmzz,mmn,pmzz)
         print ("('From array pointers for site', i3,                    &
     &            ', indexes', i3, ' and', i3, ':', i5, f7.1)"),         &
     &            j, iitest, iftest, iepr(iitest), fepr(iftest)
      end do

      print *
      print *, 'Starting spot checks of name arrays...'
      print ("('The name associated with immn =', i3, ' is ', a6)"),     &
     &         immn, iprnam(immn)
      print ("('The name associated with ipmzz=', i3, ' is ', a6)"),     &
     &         ipmzz, parnam(ipmzz)

      print *
      print *, 'Starting spot checks of error message arrays...'
      print ("('The eprerr message associated with DIMBIG=', i3, ': ',   &
     &         a32)"), DIMBIG, eprerr(DIMBIG)(1:30)
      INFO = 11
      print ("('The minerr message associated with info  =', i3, ': ',   &
     &         a32)"), info, minerr(info)(1:30)
      print *

      contains

      subroutine arbsub(jjj,kiarb,kfarb,iarb,farb)
      integer :: jjj, kiarb, kfarb, iarb
      double precision :: farb
      print ("('Via passed pointers for site', i3,                       &
     &         ', indexes', i3, ' and', i3, ':', i5, f7.1)"),            &
     &         jjj, kiarb, kfarb, iarb, farb
      farb = farb * 1.5
      iarb = iarb + 100
      print ("('Manipulated in call for site', i3,                       &
     &         ', indexes', i3, ' and', i3, ':', i5, f7.1)"),            &
     &         jjj, kiarb, kfarb, iarb, farb
      end subroutine arbsub

      end program testmods
