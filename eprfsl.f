c Version 1.6  8/12/94
c--------------------------------------------------------------------------
c                         ====================
c                          subroutine EPRFSL
c                         ====================
c
c   Field-swept conjugate gradients calculation used to determine a
c   minimum truncation set for an EPRLL slow-motional spectrum.
c   The main part of this routine consists of code from the EPRBL
c   program in versions 1.0 through 1.3 of the EPRLL program family.
c   It was separated into a subroutine because multiple calls to this
c   routine are necessitated by the MOMD capability of version 1.4.
c
c   The routine uses the conjugate gradients method to solve the matrix
c   equation x = A*b, where the diagonal of A depends on the field for
c   which the spectrum is being calculated. The problem is solved for
c   the specified field values across the spectrum, and at each point
c   the overlap of each basis vector with the solution x is calculated.
c   The maximum overlap of each vector with the solution is returned
c   and the field value at which the maximum occurred are returned
c   in the basis array.
c  
c   Upon entry to subroutine EPRFSL, it is assumed that:
c      (1) the basis set indices in /lindex/ have been initialized
c      (2) the stochastic Liouville matrix has been calculated in zmat
c      (3) the starting vector has been calculated and is passed in the
c          b vector
c      (4) if the "report" flag is true, the log file has already been 
c          opened on unit lulog
c
c     Writted by DEB using program EPRBL from version 1.3
c
c--------------------------------------------------------------------------
c
      subroutine eprfsl( b,ixfld,bsswt,spect,bmax,terrmx,nstpmx,cpufs,
     #                   ntotal,report )
      implicit none
      include 'stddim.inc'
c
      integer ixfld(MXDIM),nstpmx,ntotal
      double precision b(2,MXDIM),bsswt(MXDIM),spect(2,MXPT),bmax,
     #                 terrmx
      logical report
c
      include 'stdio.inc'
      include 'eprdat.inc'
      include 'eprmat.inc'
      include 'vectrs.inc'
      include 'pidef.inc'
      include 'fnames.inc'
      include 'indexl.inc'
      include 'timer.inc'
c
      double precision ZERO,ONE
      parameter (ZERO=0.0D0,ONE=1.0D0)
c
      logical pcflag
c
      real cpufs
      double precision invrs(MXDIM)
c
      double precision test,scale,truerr,error,field,
     #                 dfield,cputmp
c
      integer i,ifld,j,ncalc,ndone
c
      complex*16 temp,zdotu
      external zdotu
c
c#####################################################################
c
      cputmp = dtime( tarray )
c
      nstpmx=0
      ntotal=0
      bmax=ZERO
      terrmx=ZERO
      cpufs=0.0
c
c    -----------------------------
c     Zero basis projection array 
c    -----------------------------
c
      do i=1,ndim
         ixfld(i)=0
         bsswt(i)=ZERO
      end do
c
c----------------------------------------------------------------------
c     Calculate the inverse of the preconditioning matrix
c----------------------------------------------------------------------
c
      do j=1,ndim
         invrs(j)=ONE/sqrt(zdiag(1,j)+shiftr)
      end do
c

      if (report) then
         write (luttyo,3000)
         write (lulog,3000)
      end if
c
c======================================================================
c     loop over field positions
c======================================================================
c
      dfield=(fieldf-fieldi)/dble(nfield-1)
      do 200 ifld=1,nfield
c
c----------------------------------------------------------------------
c     do the preconditioned CG iteration
c----------------------------------------------------------------------
c
         do j=1,ndim
	    x(1,j)=ZERO
	    x(2,j)=ZERO
         end do
c
         field=fieldi+(ifld-1)*dfield
         temp=dcmplx(shiftr,shifti-field)
c
         pcflag=.false.
         call cspccg( b,ndim,nstep,cgtol,temp,invrs,x,ndone,
     #                error,pcflag )
c
         cputmp = dtime( tarray )
         cpufs = cpufs+cputmp
c
c---------------------------------------------------------------------
c     find the true error in the solution
c----------------------------------------------------------------------
c
         call scmvm(x,y,ndim)
         call zaxpy(x,y,temp,ndim)
c
	 truerr=ZERO
         do j=1,ndim
            y(1,j)=b(1,j)-y(1,j)
            truerr=truerr+y(1,j)*y(1,j)  
            y(2,j)=b(2,j)-y(2,j)
	    truerr=truerr+y(2,j)*y(2,j)
         end do
c
         truerr=sqrt(truerr)
c
	 if (truerr.gt.terrmx) terrmx=truerr
         ntotal=ntotal+abs(ndone)+1
c
c----------------------------------------------------------------------
c     write out messages if required
c----------------------------------------------------------------------
c
         if ((ndone.gt.0).and.report) then
            write(lulog,3010) ifld,ndone,error,truerr,cputmp
            write(luttyo,3010) ifld,ndone,error,truerr,cputmp
         endif
c
         if (ndone.le.0) then
            write(luttyo,3020) ifld,nstep,-ndone,cgtol,error,truerr,
     #                         cputmp
            if (report) write(lulog,3020) ifld,nstep,-ndone,cgtol,error,
     #                                    truerr,cputmp
         end if
         if (iabs(ndone).gt.nstpmx) nstpmx=iabs(ndone)
c
c----------------------------------------------------------------------
c     store spectral information
c----------------------------------------------------------------------
c
         temp=zdotu(b,x,ndim)
         spect(1,ifld)=spect(1,ifld)+dreal(temp)/PI
         spect(2,ifld)=spect(2,ifld)+dimag(temp)/PI
c
c----------------------------------------------------------------------
c     sift through the solution vector
c----------------------------------------------------------------------
c
         scale=1.0D0/abs(temp)
         do j=1,ndim
            test=sqrt(x(1,j)*x(1,j)+x(2,j)*x(2,j))*scale
            if (test.gt.bsswt(j)) then
               bsswt(j)=test
               ixfld(j)=ifld
            end if
            if (bsswt(j).gt.bmax) bmax=bsswt(j)
          end do
c
c======================================================================
c     end of loop over field positions
c======================================================================
c
 200    continue
c
        if (report) then 
           write(lulog,3030)
           write(lulog,3100) ntotal,terrmx
           write (lulog,3110) cpufs
c
           write(luttyo,3030)
           write(luttyo,3100) ntotal,terrmx
           write(luttyo,3110) cpufs
        end if
c
        return
c
c  ##### Format statements ############################################
c
 3000 format(2x,'+--------+-----------+----------------+',
     #            '----------------+',
     #       /,2x,'| Step   | Number of |   Modulus of   |',
     #            '  L2 norm of    |',
     #       /,2x,'| number | CG steps  |   rectanorm    |',
     #            '  true residual |',
     #       /,2x,'|--------+-----------+----------------+',
     #            '----------------|')
 3010 format(2x,'| ',i6,' |   ',i6,'  | ',g14.7,' | ',g14.7,' | ',
     #       f10.2)
 3020 format(10x,'*** CG calculation ',i3,' did not converge ***',
     #     /,10x,'*** Maximum number of CG steps  : ',i4,' ***',
     #     /,10x,'*** Number of CG steps executed : ',i4,' ***',
     #     /,10x,'*** Maximum CG error allowed    : ',g14.7,' ***',
     #     /,10x,'*** Final CG error              : ',g14.7,' ***',
     #     /,10x,'*** True error                  : ',g14.7,' ***',
     #     /,10x,'*** CPU time (seconds)          : ',g14.2,' ***')
 3030 format(2x,'+--------+-----------+----------------+',
     #          '----------------+',/)
 3100 format(2x,'Total number of maxtrix-vector multiplies : ',i5,
     #       /,2x,'Maximum true error : ',g14.7)
 3110 format(2x,'Total CPU time for swept-field calc (seconds) : ',
     #     f10.2)
        end

