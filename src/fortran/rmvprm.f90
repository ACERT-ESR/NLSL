c  Version 1.3  11/23/93
c----------------------------------------------------------------------
c                    =========================
c                       subroutine RMVPRM
c                    =========================
c
c Remove a parameter from the list of parameters being varied for nonlinear
c least-squares. Also maintain the lists of 
c    (1) ixst    : secondary parameter index (for multiple sites/spectra)
c    (2) ibnd    : boundary flag for variable parameter
c    (3) prmin   : minimum for variable parameter
c    (4) prmax   : maximum for variable parameter
c    (5) prscl   : desired accuracy for given parameter
c    (6) xfdstp  : Size of forward-differences step
c    (7) tag     : Name of the parameter (character*9)
c    (8) ixx     : index of each element in the fparm array into the x vector
c
c Notes:
c
c ix2= -1 means remove all occurrences of the parmeter in the variable list
c ix2=  0 means remove the parameter only if it applies for all sites
c
c----------------------------------------------------------------------
      subroutine rmvprm(ix,ix2,ident)
c
      use nlsdim
      use parcom
      use eprprm
      use stdio
c
      implicit none
      integer ix,ix2
      logical match
      character ident*30,tmptag*9
c
      integer i,ixa,j,j1,l
c
      integer itrim
      logical tcheck
      external itrim,tcheck
c
c######################################################################
c
      match = .false.
      ixa=abs(mod(ix,100)) 
      if (nprm.eq.0) go to 5
c
      if (ix2.le.0) then
         tmptag=ident
      else
         write (tmptag,1000) ident(:itrim(ident)),ix2
      end if
c
c----------------------------------------------------------------------
c Loop through the list of parameters being varied (it's in sort order)
c----------------------------------------------------------------------
      do i=1,nprm
         if (ixpr(i).gt.ixa) go to 5
         if (ixpr(i).eq.ixa) then
            if (ixst(i).gt.ix2.and. ix2.gt.0) go to 5
            if (ixst(i).eq.ix2 .or. ix2.eq.-1) then
c
c------------------------------------------------------------------------
c  Parameter found: first, make sure it matches the specified symmetry.
c  If it does, delete it and move up any elements below it in the list.
c------------------------------------------------------------------------
               if ( .not.tcheck(ix,ix2,ident,0) ) go to 5
c
               match=.true.
               if (ixst(i).le.0) then
                  do l=1,MXSITE
                     ixx(ixpr(i),l)=0
                  end do
               else
                  ixx(ixpr(i),ixst(i))=0
               end if
c
               do j=i,nprm-1
                  j1=j+1
                  ixpr(j)=ixpr(j1)
                  ixst(j)=ixst(j1)
                  prmin(j)=prmin(j1)
                  prmax(j)=prmax(j1)
                  prscl(j)=prscl(j1)
                  xfdstp(j)=xfdstp(j1)
                  ibnd(j)=ibnd(j1)
                  tag(j)=tag(j1) 
c
                  if (ixst(j).le.0) then
                     do l=1,MXSITE
                        ixx(ixpr(j),l)=j
                     end do
                  else
                     ixx(ixpr(j),ixst(j))=j
                  end if
               end do
c
               nprm=nprm-1
               write (luttyo,1002) tmptag(:itrim(tmptag)),nprm
c
            end if
         end if
c
      end do
c
 5    if (.not.match) write (luttyo,1001) tmptag(:itrim(tmptag))
      return
c
c #### format statements ###########################################
c
 1000 format(a,'(',i1,')')
 1001 format('*** Parameter ''',a,''' is not being varied ***')
 1002 format('*** Parameter ''',a,''' fixed: ',
     #i3,' variable parameters ***')
      end


