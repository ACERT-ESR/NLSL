c       VERSION 1.0     2/5/99
c**********************************************************************
c
c                    ===================
c                       function W3J
c                    ===================
c
c       This double precision function will calculate the values of
c       the Wigner 3-J symbols used in the Stochastic Liouville matrix
c       formulation of the slow-motional EPR calculation, i.e.
c
c                         ( J1  J2  J3 )
c                         ( M1  M2  M3 )
c
c       For J2 <= 2 and arbitrary J1 and J3, this routine explicitly
c       evaluates formulae given in Table 2 of "Angular Momentum in 
c       Quantum Mechanics" by A. R. Edmonds, revised 2nd printing, 
c       Princeton Univ. Press, 1968.
c
c       For J2 > 2 and J1,J3 < Lmax (mxlval defined in 'maxl.inc'),
c       this functions calls a modified version of code from
c       Prof. C.C.J. Roothan, function wig3j, which appears in this file. 
c
c       NOTE: Before any calls to function w3j are made, it is 
c       necessary to set up a list of binomial coefficients in
c       common block /bincom/ by calling subroutine bincof. This
c       is now handled automatically by function wig3j (see below).
c
c       Written by D.J. Schneider, 5/7/90
c       Updated by D.E.Budil 1/10/93 to fix bug for J2=1
c
c  Coded in this file:
c     w3j(n1,n2,n3,n4,n5,n6)
c     bincof()
c     wig3j(n1,n2,n3,n4,n5,n6)
c
c  Includes:
c    maxl.inc
c
c**********************************************************************
c
      function w3j(n1,n2,n3,n4,n5,n6)
      double precision w3j
      integer n1,n2,n3,n4,n5,n6
c
      include 'maxl.inc'
c
      integer j1,j2,j3,m1,m2,m3,jdelta,k
      double precision phase,parity,x,y,z,ztemp1,ztemp2
c
      double precision wig3j
c
c######################################################################
c
c----------------------------------------------------------------------
c       check triangle conditions, etc.
c----------------------------------------------------------------------
c
      if ((n1.lt.0).or.(n2.lt.0).or.(n3.lt.0).or.
     #     ((n2.gt.2).and.(n1+n2+n3+1.gt.2*(mxlval+8)+1))) then
        write (*,*) '*** quantum nubers too large in w3j ***'
        stop
      end if
c
      if ((abs(n4).gt.n1).or.(abs(n5).gt.n2).or.(abs(n6).gt.n3).or.
     #     ((n4+n5+n6).ne.0).or.
     #     ((n1+n2).lt.n3).or.((n1+n3).lt.n2).or.((n2+n3).lt.n1)) then
        w3j=0.0d0
        return
      end if
c
      j1=n1
      j2=n2
      j3=n3
      m1=n4
      m2=n5
      m3=n6
c
c----------------------------------------------------------------------
c       use wig3j to calculate if j2 > 2
c----------------------------------------------------------------------
c
      if (n2.gt.2) then
        w3j=wig3j(j1,j2,j3,m1,m2,m3)
        return
      end if
c
c----------------------------------------------------------------------
c       permute variables if necessary to get m2 => 0 and j1 <= j3,
c       keep track of phases, and store variables
c----------------------------------------------------------------------
c
      j1=n1
      j2=n2
      j3=n3
      m1=n4
      m2=n5
      m3=n6
c
      phase=1.0d0
c
      if (mod(j1+j2+j3,2).eq.0) then
        parity=1.d0
      else
        parity=-1.d0
      end if
c
      if (m2.lt.0) then
        m1=-m1
        m2=-m2
        m3=-m3
        phase=parity
      else
        phase=1.d0
      end if
c
      if(j1.gt.j3) then
        k=j1
        j1=j3
        j3=k
        k=m1
        m1=m3
        m3=k
        phase=phase*parity
      end if
c
      if (mod(j1-m3,2).ne.0) then
        phase=-phase
      end if
c
c----------------------------------------------------------------------
c        calculate wigner 3-j symbols
c----------------------------------------------------------------------
c
      jdelta=j3-j1
      x=dble(j1+j3-1)
      y=dble(j1-m3)
      z=dble(j1+m3)
c
      if (j2.eq.0) then
        w3j=phase/sqrt(dble(2*j1+1))
      else if (j2.eq.2) then
        ztemp2=x*(x+1.0d0)*(x+2.0d0)*(x+3.0d0)*(x+4.0d0)
        if (m2.eq.0) then
          if (jdelta.eq.0) then
            ztemp1=2.0d0*dble(3*m3*m3-j1*(j1+1))
            w3j=phase*ztemp1/sqrt(ztemp2)
          else if (jdelta.eq.1) then
            ztemp1=6.0d0*(z+1.0d0)*(y+1.0d0)
            w3j=-phase*2.0d0*dble(m3)*sqrt(ztemp1/ztemp2)
          else
            ztemp1=6.0d0*(z+2.0d0)*(z+1.0d0)
     #           *(y+2.0d0)*(y+1.0d0)
            w3j=phase*sqrt(ztemp1/ztemp2)
          end if
        else if (m2.eq.1) then
          if (jdelta.eq.0) then
            ztemp1=6.0d0*(z+1.0d0)*y
            w3j=phase*dble(2*m3+1)
     #           *sqrt(ztemp1/ztemp2)
          else if (jdelta.eq.1) then
            ztemp1=y*(y+1.0d0)
            w3j=-phase*dble(2*j1+4*m3+4)*
     #           sqrt(ztemp1/ztemp2)
          else
            ztemp1=(z+2.0d0)*y*(y+1.0d0)*(y+2.0d0)
            w3j=phase*2.0d0*sqrt(ztemp1/ztemp2)
          end if
        else
          if (jdelta.eq.0) then
            ztemp1=6.0d0*(y-1.0d0)*y
     #           *(z+1.0d0)*(z+2.0d0)
            w3j=phase*sqrt(ztemp1/ztemp2)
          else if (jdelta.eq.1) then
            ztemp1=(y-1.0d0)*y*(y+1.0d0)*(z+2.0d0)
            w3j=-phase*2.d0*sqrt(ztemp1/ztemp2)
          else
            ztemp1=(y-1.0d0)*y*(y+1.0d0)*(y+2.0d0)
            w3j=phase*sqrt(ztemp1/ztemp2)
          end if
        end if
c
      else
         ztemp2=(x+1.0d0)*(x+2.0d0)*(x+3.0d0)
         if (m2.eq.0) then
            if (jdelta.eq.0) then
               w3j=-phase*2.0*dble(m3)/sqrt(ztemp2)
            else
               ztemp1=2.0d0*(y+1.0d0)*(z+1.0d0)
               w3j=-phase*sqrt(ztemp1/ztemp2)
            end if
         else
            if (jdelta.eq.0) then
               ztemp1=2.0d0*y*(z+1.0d0)
               w3j=-phase*sqrt(ztemp1/ztemp2)
            else
               ztemp1=y*(y+1.0d0)
               w3j=-phase*sqrt(ztemp1/ztemp2)
            end if
         end if
      end if
c
      return
      end



c----------------------------------------------------------------------
c
c                      =========================
c                                BINCOF
c                      =========================
c
c   Generates an array of binomial coefficients for use with function
c   This must be called before wig3j can calculate a 3-J symbol,
c   and is automatically called by wig3j if the array has not been
c   initialized.
c
c                         Designed and coded by
c                        Clemens C. J. Roothaan
c                        Department of Chemistry
c                         University of Chicago
c                        5735 South Ellis Avenue
c                        Chicago, Illinois 60637
c                           December 9, 1988
c
c                   modified by DJS 5/7/90, DEB 5/20/92, DEB 1/10/93
c
c     MODIFIED 5/20/92 by DEB: added scaling to prevent overflows
c       at high values of mxlval
c       The scale factor should be the maximum of (machep,cubmax)
c       where "machep" is the machine precision and  "cubmax" is
c       the inverse cube root of the largest number that can be
c       represented on a given machine. In absence of other information,
c       set scale to unity. 
c
c     Basic vector coupling coefficient package for calculating 3n-j
c     symbols.  All the entries are fortran functions; formally:
c
c     bincof()                        generate binomial coefficients in bincom
c     wig3j(a,b,c,alpha,beta,gamma)   calculate a 3-j symbol
c----------------------------------------------------------------------

      subroutine bincof()
c
      include 'maxl.inc'
c
c --- common BINCOM ------------------------
      integer nb,nbncf
      parameter (nb=2*(mxlval+8)+2)
      parameter (nbncf=nb*(nb+1)+1)
c
      integer bncfx
      double precision bncf,scale,scal3,scal15
      common /bincom/bncf(nbncf),bncfx(nb),scale,scal3,scal15
c ------------------------------------------
c
      integer i,j,ij
      double precision bncf0,temp
c
      scale=1.0D-25
      scal3=(1.0D0/scale)**3
      scal15=dsqrt(scal3)
c
      ij=1
      do 20 i=0,nb-1
        bncfx(i+1)=ij
        if (i.ne.0) then
          bncf0=0.0d0
          do 10 j=1,i
            temp=bncf(ij-i)
            bncf(ij)=bncf0+temp
            bncf0=temp
            ij=ij+1
 10       continue
        end if
        bncf(ij)=scale
        ij=ij+1
 20   continue
c
      return
      end

c----------------------------------------------------------------------
c                    =========================
c                         function WIG3J
c                    =========================
c
c     Calculates arbitrary 3-J symbol using binomial coefficients 
c     for J values up to an limit determined by the machine precision
c     (approx. 125 for IEEE format real*8 floating-point numbers -DEB)
c     
c     This replaces the function of the same name written by G. Moro 
c     for use with EPR spectral simulations.
c
c                         Designed and coded by
c                        Clemens C. J. Roothaan
c                        Department of Chemistry
c                         University of Chicago
c                        5735 South Ellis Avenue
c                        Chicago, Illinois 60637
c                           December 9, 1988
c
c                   modified by DJS 5/7/90, DEB 1/10/93
c----------------------------------------------------------------------

      function wig3j(j1,j2,j3,m1,m2,m3)
      double precision wig3j
      integer j1,j2,j3,m1,m2,m3
      logical notset
      data notset / .true. /
c
      include 'maxl.inc'
c
c --- common BINCOM ------------------------
      integer nb,nbncf
      parameter (nb=2*(mxlval+8)+2)
      parameter (nbncf=nb*(nb+1)+1)
c
      integer bncfx
      double precision bncf,scale,scal3,scal15
      common /bincom/bncf(nbncf),bncfx(nb),scale,scal3,scal15
c ------------------------------------------
c
c......................................................................
c
c
      integer i,j,k,l,m,n,p,q,z,zmin,zmax,bp,bnj,bmk
      double precision sum
c
      if (notset) then
         call bincof
         notset = .false.
      end if
c
      i=j1+m1
      j=j1-m1
      k=j2+m2
      l=j2-m2
      m=j2+j3-j1
      n=j3+j1-j2
      p=j1+j2-j3
      q=j1+j2+j3 + 1
c
      if(q+1.gt.nb) then
        write (*,1000) j1,j2,j3,m1,m2,m3
        return
      endif
c
      bp=bncfx(p+1)
      bnj=bncfx(n+1)+n-j
      bmk=bncfx(m+1)+m-k
      zmin=max(0,j-n,k-m)
      zmax=min(p,j,k)
c
      sum=0.0d0
        do 30 z=zmin,zmax
        sum=-sum+bncf(bp+z)*bncf(bnj+z)*bncf(bmk+z)*scal3
 30     continue
c
      if (sum.ne.0.0d0) then
        if (mod(i+l+zmax,2).ne.0) sum=-sum
        sum=sum*sqrt(bncf(bncfx(q)+m)/
     &               bncf(bncfx(q)+i))/
     &          sqrt(bncf(bncfx(q-i)+k))*
     &          sqrt(bncf(bncfx(q-m)+n)/
     &               bncf(bncfx(q)+j))/
     &          sqrt(bncf(bncfx(q-j)+l))/
     &          sqrt(bncf(bncfx(q+1)+1))  
      endif
c
      wig3j=sum/scal15
      return
c
c######################################################################
c
 1000   format(' wig3j called with (',6(i3,','),'):'/
     #         '   sum of L values exceeds limit')
      end
