c----------------------------------------------------------------------
c                    =========================
c                       common block INDEXF
c                    =========================
c
c  Common block containing a list of the basis set indices
c  for off-diagonal and diagonal subspace.
c
c pidptr - points to location within pi1 etc of basis set information
c dpidptr - points into dpi1 etc.
c basisptr - indicates which basis is currently mapped into jqe1 etc 
c
c----------------------------------------------------------------------
c
      integer jqe1,pi1,qi1,l1,jk1,k1,jm1,m1
      integer mjqe1,mpi1,mqi1,ml1,mjk1,mk1,mjm1,mm1,
     #		pidptr,dpidptr,basisptr
c allow for multiple sites, but not multiple spec basis sets.
      common /indexo/ jqe1(MXDIM),pi1(MXDIM),qi1(MXDIM),
     #        l1(MXDIM),jk1(MXDIM),k1(MXDIM),jm1(MXDIM),
     #        m1(MXDIM)
      common /mindexo/ mjqe1(MMXDIM),mpi1(MMXDIM),mqi1(MMXDIM),
     #        ml1(MMXDIM),mjk1(MMXDIM),mk1(MMXDIM),mjm1(MMXDIM),
     #        mm1(MMXDIM),pidptr(MXSPEC*MXSITE),basisptr
c
      integer djqe1,dpi1,dqi1,dl1,djk1,dk1,djm1,dm1
      integer mdjqe1,mdpi1,mdqi1,mdl1,mdjk1,mdk1,mdjm1,mdm1
      common /indexd/ djqe1(2*MXDIM),dpi1(2*MXDIM),dqi1(2*MXDIM),
     #        dl1(2*MXDIM),djk1(2*MXDIM),dk1(2*MXDIM),djm1(2*MXDIM),
     #        dm1(2*MXDIM)
      common/mindexd/mdjqe1(2*MMXDIM),mdpi1(2*MMXDIM),mdqi1(2*MMXDIM),
     #        mdl1(2*MMXDIM),mdjk1(2*MMXDIM),mdk1(2*MMXDIM),
     #        mdjm1(2*MMXDIM),mdm1(2*MMXDIM),dpidptr(MXSPEC*MXSITE)
c
c----------------------------------------------------------------------
c                    =========================
c                       common block PULSEP
c                    =========================
c
c  Common block containining the pulse propagator
c
c  pp(i)  : the actual pulse propagator element for the i-th basis
c           element in the diagonal space.
c  pid(j) : the number of the corresponding basis elements in the
c           diagonal space for the j-th basis element in off-diagonal
c           basis.  It can be either 1 or 2.
c
c----------------------------------------------------------------------
c
      integer pid,mpid
      double precision pp,mpp
      common /pulsep/ pid(MXDIM),pp(2*MXDIM),mpid(MMXDIM),mpp(2*MMXDIM)
