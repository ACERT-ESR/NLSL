c  VERSION 1.0  (NLSPMC version)   2/5/99
c----------------------------------------------------------------------
c                    =========================
c                       subroutine RMVPRM
c                    =========================
c
c Remove a parameter from the list of parameters being varied for nonlinear
c least-squares. Also maintain the lists of 
c    (1) ibnd    : boundary flag for variable parameter
c    (2) prmin   : minimum for variable parameter
c    (3) prmax   : maximum for variable parameter
c    (4) prscl   : desired accuracy for given parameter
c    (5) xfdstp  : Size of forward-differences step
c    (6) tag     : Name of the parameter (character*6)
c    (7) ixx     : index of each element in the fparm array into the x vector
c
c ixr identifies variable, with axis identification,
c ident is the the location in xxx of the variable to be removed.
c Note the xxx array itself is set elsewhere on calling LM routine
c
c Includes:
c
c      limits.inc
c      parms.inc
c      simparm.inc
c      stdio.inc
c      miscel.inc
c
c Spectrum and site information is checked against current variable 
c parameter specification.
c
c----------------------------------------------------------------------
      subroutine rmvprm(ident)
      implicit none
      integer varid,ident,icomp,ispec
c
      include 'limits.inc'
      include 'parms.inc'
      include 'simparm.inc'
      include 'stdio.inc'
      include 'miscel.inc'
c
      integer j,j1
      character*6 idparm
c
c######################################################################
c
      if (nprm.eq.0) then
        write (luttyo,1000) 
 1000   format('*** No parameters are being varied ***')
	return	! if no variable parameters, return
      end if
c
      idparm=tag(ident)		! get ascii parameter name
      write (luttyo,1001) idparm,nprm-1
 1001 format('*** Parameter ''',a,''' fixed: ',
     #i3,' variable parameters remaining ***')
c
c  Delete it and move up elements below it in the list.
c  (Note that move loop will not execute if we are removing last element;
c  i.e. if ident=nprm)
c------------------------------------------------------------------------
c
      varid=ixpr(ident)		! identity of the variable
      do 1 j=ident,nprm-1
        j1=j+1
        ixpr(j)=ixpr(j1)
        prmin(j)=prmin(j1)
        prmax(j)=prmax(j1)
        prscl(j)=prscl(j1)
        xfdstp(j)=xfdstp(j1)
        ibnd(j)=ibnd(j1)
        tag(j)=tag(j1) 
 1      continue
        do 2 icomp=1,ncomps
          do 2 ispec=1,nspectra
c zero ixx pointer for parameter being fixed
            if (ixx(varid,ispec,icomp).eq.ident) then
              ixx(varid,ispec,icomp)=0 
            end if
 2      continue
c decrement ixx pointer for moved parms
        do 3 j=1,nprm-1
          do 3 icomp=1,ncomps
            do 3 ispec=1,nspectra
              if (ixx(ixpr(j),ispec,icomp).gt.ident) 
     #         ixx(ixpr(j),ispec,icomp)=ixx(ixpr(j),ispec,icomp)-1
 3      continue
c
      nprm=nprm-1
      return
c
c #### format statements ###########################################
c
      end
