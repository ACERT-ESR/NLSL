! NLSL Version 1.5 beta 11/24/95
!----------------------------------------------------------------------
!                    =========================
!                      subroutine SITEC
!                    =========================
!
! sites <n>
!
!      n   : number of different "sites" or components in each spectrum
!
!----------------------------------------------------------------------
      subroutine sitec( line )
!
      use nlsdim
      use expdat
      use parcom
      use tridag
      use basis
      use stdio
!
      implicit none
      character*80 line
!
      character token*30
!
      integer i,ixpar,ixsite,j,lth,ival
!
      integer itrim
      logical itoken
      external itrim,itoken
!
!----------------------------------------------------------------------
! Get the number of sites required
!----------------------------------------------------------------------
!
      call gettkn(line,token,lth)
!
      if (lth.gt.0 .and. itoken(token,lth,ival)) then

         if (ival.le.MXSITE.and.ival.ge.1) then
!
!           ------------------------------------------------------------
!           If the number of sites is being reduced, mark all basis
!           sets and tridiagonal matrices belonging to the sites that
!           are being removed 
!           ------------------------------------------------------------
            if (ival.lt.nsite) then
               do i=ival+1,nsite
                  do j=1,MXSPC
                     bsused( basno(i,j) )=0
                     modtd(i,j)=1
                  end do
               end do
            end if
!
!           --------------------------------------------------
!           Remove any variable parameter belonging to a site
!           that is no longer defined
!           --------------------------------------------------
            do i=1,nprm
               if (ixst(i).gt.ival) then
                  ixpar=ixpr(i)
                  ixsite=ixst(i)
                  token=tag(i)
                  call rmvprm(ixpar,ixsite,token)
               end if
            end do
!
            nsite=ival
!
!
!           ------------------------------------------------------------
!           Issue warning message if some of the data in a multiple-site
!           multiple spectrum fit are not normalized
!           ------------------------------------------------------------
            if (nsite.gt.1 .and. nser.gt.1) then
               do i=1,nser
                  if (nrmlz(i).eq.0) 
     #               write (luttyo,1002) i,dataid(i)(:itrim(dataid(i)))
               end do
            end if
!
!                            *** Bad # of sites
         else
            write(luttyo,1000) MXSITE
         end if
!
!                   *** Expecting integer value
      else
         write (luttyo,1001)
      end if
      return
!
! ###### format statements ########################################
!
 1000 format('*** Number of sites must be ',i2,' or less ***')
 1001 format('*** Expecting integer value ***')
 1002 format('*** WARNING: spectrum ',i2,' (',a,
     #     ') not normalized ***')
      end
