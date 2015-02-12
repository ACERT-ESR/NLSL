! NLSL Version 1.9.0 beta 2/12/15
!
      program testmods

      use eprprm
      use errmsg
      use lpnam
      implicit none

      integer, external :: ipfind, isfind, itrim

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

      integer :: i, j, itest, iitest, iftest, info, isel, ltok, lnam
      double precision :: ftest
      character*6 :: token='      '

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
         iitest = IMMN
         iftest = IPMZZ
         itest = mmn
         ftest = pmzz
         print ("('The expected values for site', i3,                    &
     &            ', indexes', i3, ' and', i3, ':', i5, f7.1)"),         &
     &            j, iitest, iftest, iparm(iitest,j), fparm(iftest,j)
         print ("('Fetched by pointers for site', i3,                    &
     &            ', indexes', i3, ' and', i3, ':', i5, f7.1)"),         &
     &            j, iitest, iftest, itest, ftest
         ! now try some tests from within a subroutine
         call arbsub(j,IMMN,IPMZZ,mmn,pmzz)
         print ("('From array pointers for site', i3,                    &
     &            ', indexes', i3, ' and', i3, ':', i5, f7.1)"),         &
     &            j, iitest, iftest, iepr(iitest), fepr(iftest)
      end do

      print *
      print *, 'Starting spot checks of lpnam name arrays...'
      isel = IMMN
      token = iprnam(isel)
      ltok = itrim(token)
      lnam = ipfind(token,ltok)
      print ("('The iprnam associated with index', i3, ' is ', a4)"),    &
     &         isel, token(:ltok)
      print ("('When given token ', a4, ', ipfind returns ', i3,         &
     &         ', equating to index ', i3)"),                            &
     &         token(:ltok), lnam, lnam-100
      isel = IPMZZ
      token = parnam(isel)
      ltok = itrim(token)
      lnam = ipfind(token,ltok)
      print ("('The parnam associated with index', i3, ' is ', a4)"),    &
     &         isel, token(:ltok)
      print ("('When given token ', a4, ', ipfind returns ', i3,         &
     &         ', equating to index ', i3)"),                            &
     &         token(:ltok), lnam, lnam
      isel = 4
      token = symbol(4)
      ltok = itrim(token)
      lnam = isfind(token,ltok)
      print ("('The symbol associated with index', i3, ' is ', a4)"),    &
     &         isel, token(:ltok)
      print ("('When given token ', a4, ', isfind returns ', i3,         &
     &         ', equating to index ', i3)"),                            &
     &         token(:ltok), lnam, lnam

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
