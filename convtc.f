c Version 1.4 10/10/94
c----------------------------------------------------------------------
c                         ====================
c                          subroutine CONVTC
c                         ====================
c
c     Subroutine to convert the tensor named on the passed command line
c     according to the symmetry specified by iflg. The following are
c     valid symmetries: 
c         1:  Cartesian
c         2:  Spherical
c         3:  Axial
c     See source file tensym.f for an explanation of these symmetries.
c
c     Uses:
c         gettkn, touppr (in strutl.f)
c         ipfind.f
c         itrim.f
c         indtkn.f
c         rmvprm.f
c         tocart, tosphr, toaxil (in tensym.f)
c
c---------------------------------------------------------------------- 
c
      subroutine convtc( line, iflg )
      implicit none
c
      character line*80, token*30, ident*9, varstr*9, tnstr*9
      integer i,iflg,ix,ixa,ixf,ixten,ix2,j,jx,jx1,jx2,lth,k
c
      include 'nlsdim.inc'
      include 'eprprm.inc'
      include 'parcom.inc'
      include 'stdio.inc'
      include 'lpnam.inc'
      include 'symdef.inc'
c
      integer ipfind,indtkn,itrim
      external ipfind,indtkn,itrim
c
c######################################################################
c
      if (iflg.lt.CARTESIAN .or. iflg.gt.AXIAL) then
         write (luttyo,1006)
         return
      end if
c
      call gettkn(line,token,lth)
      call touppr(token,lth)
c
c                                          *** No tensor specified
      if (lth.eq.0) then

         write (luout,1000)
         return
      end if
c
      lth=1
      ix=ipfind(token,lth)
      ixa=abs(mod(ix,100))
      if (ix.gt.100) ixa=0
c     
      ix2=indtkn(line)
      if (ix2.le.0) then
         jx1=1
         jx2=MXSITE
      else
         jx1=ix2
         jx2=ix2
      end if
c
      if (ixa.eq.0 .or. ixa-IWXX.ge.NALIAS) then
c
c                                          *** Unknown tensor
         write (luout,1001) token(1:1)
         return
c
c----------------------------------------------------------------------
c     Tensor found: Check existing symmetry
c----------------------------------------------------------------------
      else
         ixf=IIWFLG+(ixa-IWXX)/3
         ixten=IWXX+3*(ixf-IIWFLG)
c
c----------------------------------------------------------------------
c       Check whether any tensor components of another symmetry 
c       are in the list of variable parameters
c----------------------------------------------------------------------
c     
         do 12 jx=jx1,jx2
            if (ix2.le.0) then
               tnstr=token(1:1)
            else
               write(tnstr,1004) token(1:1),jx
            end if
c
            do i=0,2
               j=ixten+i
               if (ixx(j,jx).ne.0) then
                  varstr=tag(ixx(j,jx))
                  write (luttyo,1005) tnstr(:itrim(tnstr)),
     #                                varstr(:itrim(varstr))
                  return
               end if 
            end do
c
c----------------------------------------------------------------------
c           Symmetry not set yet: set it
c----------------------------------------------------------------------
            if(iparm(ixf,jx) .eq. 0) then
               iparm(ixf,jx)=iflg
               if (jx.eq.jx1) 
     #              write(luttyo,1002) tnstr(:itrim(tnstr)),
     #                              symstr(iflg)(:itrim(symstr(iflg)))
               go to 12
c
c----------------------------------------------------------------------
c           Symmetry setting is same as the one specified: skip
c----------------------------------------------------------------------
            else if (iparm(ixf,jx).eq.iflg) then
               go to 12
            end if
c
c----------------------------------------------------------------------
c           Now...convert the tensor symmetry!
c----------------------------------------------------------------------
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
c
      return
c     
 1000 format('*** No tensor specified ***')
 1001 format('*** Unknown tensor: ''',a,''' ***')
 1002 format('*** ',a,' tensor set to ',a,' ***')
 1003 format('*** ',a,' tensor converted to ',a,' ***')
 1004 format(a1,'(',i1,') ')
 1005 format('*** ',a,' tensor symmetry unchanged: ',a,
     # ' is being varied ***')
 1006 format('*** CONVERT called with illegal symmetry type ***')
      end
