c NLSL Version 1.5 beta 11/24/95
c----------------------------------------------------------------------
c                    =========================
c                      subroutine SITEC
c                    =========================
c
c sites <n>
c
c      n   : number of different "sites" or components in each spectrum
c
c----------------------------------------------------------------------
      subroutine sitec( line )
c
      use nlsdim
      use expdat
      use parcom
      use tridag
      use basis
      use stdio
c
      implicit none
      character*80 line
c
      character token*30
c
      integer i,ixpar,ixsite,j,lth,ival
c
      integer itrim
      logical itoken
      external itrim,itoken
c
c----------------------------------------------------------------------
c Get the number of sites required
c----------------------------------------------------------------------
c
      call gettkn(line,token,lth)
c
      if (lth.gt.0 .and. itoken(token,lth,ival)) then

         if (ival.le.MXSITE.and.ival.ge.1) then
c
c           ------------------------------------------------------------
c           If the number of sites is being reduced, mark all basis
c           sets and tridiagonal matrices belonging to the sites that
c           are being removed 
c           ------------------------------------------------------------
            if (ival.lt.nsite) then
               do i=ival+1,nsite
                  do j=1,MXSPC
                     bsused( basno(i,j) )=0
                     modtd(i,j)=1
                  end do
               end do
            end if
c
c           --------------------------------------------------
c           Remove any variable parameter belonging to a site
c           that is no longer defined
c           --------------------------------------------------
            do i=1,nprm
               if (ixst(i).gt.ival) then
                  ixpar=ixpr(i)
                  ixsite=ixst(i)
                  token=tag(i)
                  call rmvprm(ixpar,ixsite,token)
               end if
            end do
c
            nsite=ival
c
c
c           ------------------------------------------------------------
c           Issue warning message if some of the data in a multiple-site
c           multiple spectrum fit are not normalized
c           ------------------------------------------------------------
            if (nsite.gt.1 .and. nser.gt.1) then
               do i=1,nser
                  if (nrmlz(i).eq.0) 
     #               write (luttyo,1002) i,dataid(i)(:itrim(dataid(i)))
               end do
            end if
c
c                            *** Bad # of sites
         else
            write(luttyo,1000) MXSITE
         end if
c
c                   *** Expecting integer value
      else
         write (luttyo,1001)
      end if
      return
c
c ###### format statements ########################################
c
 1000 format('*** Number of sites must be ',i2,' or less ***')
 1001 format('*** Expecting integer value ***')
 1002 format('*** WARNING: spectrum ',i2,' (',a,
     #     ') not normalized ***')
      end
