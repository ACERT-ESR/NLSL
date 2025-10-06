!  NLSL Version 1.5 beta 11/23/95
!----------------------------------------------------------------------
!                    =========================
!                       subroutine ADDPRM
!                    =========================
!
!>@brief Add a parameter to list of parameters being varied for nonlinear
!> least-squares. Also maintain the lists of 
!>
!> 1. ixst    : secondary parameter index (for multiple sites/spectra)
!> 2. ibnd    : boundary flag for variable parameter
!> 3. prmin   : minimum for variable parameter
!> 4. prmax   : maximum for variable parameter
!> 5. prscl   : desired accuracy for given parameter
!> 6. xfdstp  : Size of forward-differences step
!> 7. tag     : Name of the parameter (character*9)
!> 8. ixx     : index of each element in the fparm array into the x vector
!>  
!> The values of these parameters for the given variable parameter
!> are passed to ADDPRM:
!>      subroutine parameter |kept in list
!>      ---------------------|------------
!>      ix2                  |ixst
!>      ibd                  |ibnd
!>      prmn                 |prmin
!>      prmx                 |prmax
!>      prsc                 |prscl
!>      step                 |xfdstp
!>      ident                |tag
!> 
!> Notes: 
!> 
!>  - Each addition to the list is included in sort order with respect to
!>     the index ix. This is an artifact of the original implementation
!>     of the program, and no longer required.
!> 
!>  - If ix2 = -1, the parameter specified by ix is to be varied individually
!>     for EACH existing site/spectrum
!>     If ix2 =  0, the parameter is to be varied for ALL sites/spectra at once
! 
!     Includes:  
!        nlsdim.inc
!        eprprm.inc
!        expdat.inc
!        parcom.inc
!        lpnam.inc
!        lmcom.inc
!        stdio.inc
!        rndoff.inc
!        prmeqv.inc
!
!----------------------------------------------------------------------
      subroutine addprm(ix,ix2,ibd,prmn,prmx,prsc,step,ident)
!
      use nlsdim
      use eprprm
      use expdat
      use parcom
      use lpnam
      use lmcom
      use stdio
      use rnddbl
!      use prmeqv
!
      implicit none
      integer ix,ix2,ibd
      double precision prmn,prmx,prsc,step
      character*9 ident
!
      integer i,ixabs,j,jx,jx1,jx2,k,k1,l,lu,nadd,nmov
!
      integer itrim
      logical spcpar,tcheck
      external itrim,spcpar,tcheck
!
!######################################################################
!
!
!  --- Check index bounds
!
      if (    (spcpar(ix) .and. ix2.gt.MXSPC)
     #    .or. (.not.spcpar(ix) .and. ix2.gt.MXSITE) ) then 
        write (luttyo,1002) 
        return
      end if
!
!----------------------------------------------------------------------
! ix2= -1: vary the parameter individually for EACH exiting site/spectrum
! ix2=  0: vary the same parameter for ALL existing sites/spectra
!----------------------------------------------------------------------
      if (ix2.lt.0) then
         jx1=1
         if (spcpar(ix)) then
            jx2=nser
         else
            jx2=nsite
         end if
      else
         jx1=ix2
         jx2=ix2
      end if
!
      nadd=jx2-jx1+1
!
!----------------------------------------------------------------------
! See if there is enough room in the parameter list
!----------------------------------------------------------------------
      if (nprm+nadd.gt.MXVAR) then
        write(luttyo,1001) MXVAR
        return
      end if
!
!----------------------------------------------------------------------
! Loop over all parameters to be added
!----------------------------------------------------------------------
      lu=luttyo
      do 21 jx=jx1,jx2
!
      if (.not.tcheck(ix,jx,ident,lu)) return
      lu=0
!
!----------------------------------------------------------------------
! Search parameter list to see if parameter is already there
!----------------------------------------------------------------------
        ixabs=abs(mod(ix,100))
        do i=1,nprm
!                                        *** parameter already in list
           if (ixpr(i).eq.ixabs .and. ixst(i).eq.jx) then
              write(luttyo,1000) tag(i)
              return
           end if
!                                        *** parameter not in list
           if ( ixpr(i).gt.ixabs .or. 
     #          (ixpr(i).eq.ixabs.and.ixst(i).gt.jx) ) go to 14
        end do
!
!----------------------------------------------------------------------
! Reached end of list: parameter will be added to the end
!----------------------------------------------------------------------
        i=nprm+1
        nprm=nprm+1
        go to 18
!
!----------------------------------------------------------------------
! Found proper location for new parameter: move top of list up
!----------------------------------------------------------------------
   14 nmov=nprm-i+1
        do j=1,nmov
           k=nprm-j+1
           k1=k+1
           x(k1)=x(k)
           ixpr(k1)=ixpr(k)
           ixst(k1)=ixst(k)
           prmin(k1)=prmin(k)
           prmax(k1)=prmax(k)
           prscl(k1)=prscl(k)
           xfdstp(k1)=xfdstp(k)
           ibnd(k1)=ibnd(k)
           tag(k1)=tag(k)
           if (ixst(k1).le.0) then
              do l=1,MXSITE
                 ixx(ixpr(k1),l)=k1
              end do
           else
              ixx(ixpr(k1),ixst(k1))=k1
           end if
!
       end do
       nprm=nprm+1
!
!----------------------------------------------------------------------
!  Add new parameter(s) to proper location in list
!----------------------------------------------------------------------
 18     if (jx.eq.0) then
           x(i)=fparm(ixabs,1)
           do l=1,MXSITE
              ixx(ixabs,l)=i
           end do
        else
           x(i)=fparm(ixabs,jx)
           ixx(ixabs,jx)=i
        end if
        ixpr(i)=ixabs
        ixst(i)=jx
        prmin(i)=prmn
        prmax(i)=prmx
        prscl(i)=prsc
        xfdstp(i)=step*x(i)
        if(dabs(xfdstp(i)).lt.RNDOFF) xfdstp(i)=step
        ibnd(i)=ibd
        tag(i)=ident
!
        if (jx.ne.0) then
           k=itrim(tag(i))
           write(tag(i)(k+1:),2000) jx
        end if
!
 21   continue
!
      write(luttyo,1003) nprm
      return
!
! #### format statements ###########################################
!
 1000 format('*** Parameter ''',a,''' is already being varied ***')
 1001 format('*** Limit of',i2,' variable parameters exceeded ***')
 1002 format('*** Index out of range ***')
 1003 format(i3,' parameters now being varied')
 2000 format('(',i1,')') 
      end
