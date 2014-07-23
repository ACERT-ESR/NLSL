c  VERSION 1.0  (NLSPMC version)   2/5/99
c*********************************************************************
c
c                       ===============
c                       SUBROUTINE LUMP
c                       ===============
c
c     Subroutine LUMP combines the eigenvalues of T-matrix using the
c     relative tolerance reltol.  Lumping is necessary because it is
c     impossible to accurately precict the accuracy of the CMTQLI
c     routine.  LUMP combines T-eigenvalues that have slipped by the
c     tolerance that was used in the T-multiplicity tests.  In parti-
c     cular if for some j,
c
c          |eval(j)-eval(j-1)| < max(reltol*|eval(j)|,scale2*multol)
c
c     then these T-eigenvalues are combined.  multol is the tolerance
c     that was used in the T-multiplicity test in COMPEV.
c
c     If in a set of T-eigenvalues to be combined there is an eigenvalue
c     with lindex=1, then the value of the combined T-eigenvalues is set
c     equal to the value of that eigenvalue.  Note that if a spurious
c     T-eigenvalue is to be 'combined' with a good eigenvalue, then this
c     is done only by increasing the index, lindex, for that eigenvalue
c     numerical values of spurious T-eigenvalues are never combined with
c     those of good T-eigenvalues.
c
c     arguments :
c
c         vc(j) = jth distinct T-eigenvalue
c         va(j) = |vc(j)|, in order of increasing magnitude
c         lindex(j) = T-multiplicity of jth distinct T-eigenvalue
c         loop = number of distinct T-eigenvalues
c         value of reltol is 1.d-8.
c
c     Includes :
c               nlsdim.inc
c               eprprm.inc
c
c*********************************************************************
c
      subroutine lump(vc,v1,va,w,reltol,sputol,scale2,lindex,
     #                tflag,loop)
c
      implicit none
c
      include 'limits.inc'
      include 'simparm.inc'
c
      complex*16 vc(mxstep),v1(mxstep),w(mxstep),csum,wsum,czero
      complex*16 wt(100)
c
      double precision  va(mxstep),reltol,sputol,scale2
      double precision  thold,th1,th2,dgap,one
      parameter (one=1.0D0,czero=(0.0D0,0.0D0))
c
      integer  i,icount,idif,in,indsum,ispur,j,jn,k,loop,nloop
      integer  lindex(mxstep),tflag(mxstep)
c
c-----------------------------------------------------------------------
c
      th2=scale2*sputol
c
      do 10 k=1,loop
 10   tflag(k)=0

      nloop=0
      j=0
c
 20   j=j+1
      if (j.gt.loop) go to 130
      if (tflag(j).eq.1) go to 20
      nloop=nloop+1
      tflag(j)=1
      v1(1)=vc(j)
      wt(1)=w(j)
      icount=1
      jn=lindex(j)
      th1=reltol*va(j)
      thold=dmax1(th1,th2)
c     thold=reltol*dmax1(one,va(j))
c
      if (jn.ne.0) then
        indsum=jn
        ispur=0
        csum=dfloat(jn)*vc(j)
        wsum=dfloat(jn)*w(j)
      else
        indsum=1
        ispur=1
        csum=czero
        wsum=czero
      end if
c
      if (j.eq.loop) go to 70
      i=j
 50   i=i+1
      if (i.gt.loop) go to 70
      if (tflag(i).eq.1) go to 50
      dgap=va(i)-va(j)
      if (dgap.ge.thold) go to 70
      dgap=cdabs(vc(i)-vc(j))
      if (dgap.ge.thold) go to 50
c
c     lump vc(i) with vc(j)
c
      icount=icount+1
      tflag(i)=1
      v1(icount)=vc(i)
      wt(icount)=w(i)      
      in=lindex(i)
c
      if (in.eq.0) then
        ispur=ispur+1
        indsum=indsum+1
      else
        indsum=indsum+in
        csum=csum+dfloat(in)*vc(i)
        wsum=wsum+dfloat(in)*w(i)
      end if
      go to 50
c
c     compute the 'combined' T-eigenvalue and the resulting
c     T-multiplicity
c
 70   continue
c
      if (icount.eq.1) indsum=jn
c
      idif=indsum-ispur
c
      if (icount.eq.1) then
        vc(nloop)=vc(j)
        va(nloop)=va(j)
        w(nloop)=w(j)
        lindex(nloop)=indsum
      else if (idif.ne.0) then
        csum=csum/dfloat(idif)
        wsum=wsum/dfloat(idif)
        vc(nloop)=csum
        va(nloop)=cdabs(csum)
        w(nloop)=wsum
        lindex(nloop)=indsum
      else
        do 90 k=1,icount
          vc(nloop+k-1)=v1(k)
          va(nloop+k-1)=cdabs(v1(k))
          w(nloop+k-1)=wt(k)
 90     continue
        nloop=nloop+icount-1
      end if
c
      go to 20
c
c     index j is finished
c
c     on return vc contains the distinct T-eigenvalues  va=|vc|
c     lindex contains the corresponding T-multiplicities
c
 130  continue
      loop=nloop
c
      return
c
      end
