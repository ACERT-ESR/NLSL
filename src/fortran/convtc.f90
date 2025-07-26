! Version 1.4 10/10/94
!----------------------------------------------------------------------
!                         ====================
!                          subroutine CONVTC
!                         ====================
!
!     Subroutine to convert the tensor named on the passed command line
!     according to the symmetry specified by iflg. The following are
!     valid symmetries: 
!         1:  Cartesian
!         2:  Spherical
!         3:  Axial
!     See source file tensym.f for an explanation of these symmetries.
!
!     Uses:
!         gettkn, touppr (in strutl.f)
!         ipfind.f
!         itrim.f
!         indtkn.f
!         rmvprm.f
!         tocart, tosphr, toaxil (in tensym.f)
!
!---------------------------------------------------------------------- 
!
      subroutine convtc( line, iflg )
!
      use nlsdim
      use eprprm
      use parcom
      use stdio
      use lpnam
      use symdef
!
      implicit none
!
      character line*80, token*30, ident*9, varstr*9, tnstr*9
      integer i,iflg,ix,ixa,ixf,ixten,ix2,j,jx,jx1,jx2,lth,k
!
      integer ipfind,indtkn,itrim
      external ipfind,indtkn,itrim
!
!######################################################################
!
      if (iflg.lt.CARTESIAN .or. iflg.gt.AXIAL) then
         write (luttyo,1006)
         return
      end if
!
      call gettkn(line,token,lth)
      call touppr(token,lth)
!
!                                          *** No tensor specified
      if (lth.eq.0) then

         write (luout,1000)
         return
      end if
!
      lth=1
      ix=ipfind(token,lth)
      ixa=abs(mod(ix,100))
      if (ix.gt.100) ixa=0
!     
      ix2=indtkn(line)
      if (ix2.le.0) then
         jx1=1
         jx2=MXSITE
      else
         jx1=ix2
         jx2=ix2
      end if
!
      if (ixa.eq.0 .or. ixa-IWXX.ge.NALIAS) then
!
!                                          *** Unknown tensor
         write (luout,1001) token(1:1)
         return
!
!----------------------------------------------------------------------
!     Tensor found: Check existing symmetry
!----------------------------------------------------------------------
      else
         ixf=IIWFLG+(ixa-IWXX)/3
         ixten=IWXX+3*(ixf-IIWFLG)
!
!----------------------------------------------------------------------
!       Check whether any tensor components of another symmetry 
!       are in the list of variable parameters
!----------------------------------------------------------------------
!     
         do 12 jx=jx1,jx2
            if (ix2.le.0) then
               tnstr=token(1:1)
            else
               write(tnstr,1004) token(1:1),jx
            end if
!
            do i=0,2
               j=ixten+i
               if (ixx(j,jx).ne.0) then
                  varstr=tag(ixx(j,jx))
                  write (luttyo,1005) tnstr(:itrim(tnstr)),
     #                                varstr(:itrim(varstr))
                  return
               end if 
            end do
!
!----------------------------------------------------------------------
!           Symmetry not set yet: set it
!----------------------------------------------------------------------
            if(iparm(ixf,jx) .eq. 0) then
               iparm(ixf,jx)=iflg
               if (jx.eq.jx1) 
     #              write(luttyo,1002) tnstr(:itrim(tnstr)),
     #                              symstr(iflg)(:itrim(symstr(iflg)))
               go to 12
!
!----------------------------------------------------------------------
!           Symmetry setting is same as the one specified: skip
!----------------------------------------------------------------------
            else if (iparm(ixf,jx).eq.iflg) then
               go to 12
            end if
!
!----------------------------------------------------------------------
!           Now...convert the tensor symmetry!
!----------------------------------------------------------------------
            if (iflg.eq.SPHERICAL) then
               call tosphr( fparm(ixten,jx),iparm(ixf,jx) )
            else if (iflg.eq.AXIAL) then
               call toaxil( fparm(ixten,jx),iparm(ixf,jx) )
            else
               call tocart( fparm(ixten,jx),iparm(ixf,jx) )
            end if
            iparm(ixf,jx)=iflg
            if (jx.eq.jx1) 
     #         write (luttyo,1003) tnstr(:itrim(tnstr)),
     #                             symstr(iflg)(:itrim(symstr(iflg)))
 12      continue
      end if
!
      return
!     
 1000 format('*** No tensor specified ***')
 1001 format('*** Unknown tensor: ''',a,''' ***')
 1002 format('*** ',a,' tensor set to ',a,' ***')
 1003 format('*** ',a,' tensor converted to ',a,' ***')
 1004 format(a1,'(',i1,') ')
 1005 format('*** ',a,' tensor symmetry unchanged: ',a,
     # ' is being varied ***')
 1006 format('*** CONVERT called with illegal symmetry type ***')
      end
