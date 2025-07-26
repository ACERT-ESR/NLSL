!       VERSION 1.1     1/12/93
!**********************************************************************
!
!                    ===================
!                       function W3J
!                    ===================
!
!       This double precision function will calculate the values of
!       the Wigner 3-J symbols used in the Stochastic Liouville matrix
!       formulation of the slow-motional EPR calculation, i.e.
!
!                         ( J1  J2  J3 )
!                         ( M1  M2  M3 )
!
!       For J2 <= 2 and arbitrary J1 and J3, this routine explicitly
!       evaluates formulae given in Table 2 of "Angular Momentum in 
!       Quantum Mechanics" by A. R. Edmonds, revised 2nd printing, 
!       Princeton Univ. Press, 1968.
!
!       For J2 > 2 and J1,J3 < Lmax (mxlval defined in 'maxl.inc'),
!       this functions calls a modified version of code from
!       Prof. C.C.J. Roothan, function wig3j, which appears in this file. 
!
!       NOTE: Before any calls to function w3j are made, it is 
!       necessary to set up a list of binomial coefficients in
!       common block /bincom/ by calling subroutine bincof. This
!       is now handled automatically by function wig3j (see below).
!
!       Written by D.J. Schneider, 5/7/90
!       Updated by D.E.Budil 1/10/93 to fix bug for J2=1
!
!  Coded in this file:
!     w3j(n1,n2,n3,n4,n5,n6)
!     bincof()
!     wig3j(n1,n2,n3,n4,n5,n6)
!
!  Includes:
!    maxl.inc
!
!**********************************************************************
!
      function w3j(n1,n2,n3,n4,n5,n6)
!
      use maxl
!
      double precision w3j
      integer n1,n2,n3,n4,n5,n6
!
      integer j1,j2,j3,m1,m2,m3,jdelta,k
      double precision phase,parity,x,y,z,ztemp1,ztemp2
!
      double precision wig3j
!
!######################################################################
!
!----------------------------------------------------------------------
!       check triangle conditions, etc.
!----------------------------------------------------------------------
!
      if ((n1.lt.0).or.(n2.lt.0).or.(n3.lt.0).or.
     #     ((n2.gt.2).and.(n1+n2+n3+1.gt.2*(mxlval+8)+1))) then
        write (*,*) '*** quantum nubers too large in w3j ***'
        stop
      end if
!
      if ((abs(n4).gt.n1).or.(abs(n5).gt.n2).or.(abs(n6).gt.n3).or.
     #     ((n4+n5+n6).ne.0).or.
     #     ((n1+n2).lt.n3).or.((n1+n3).lt.n2).or.((n2+n3).lt.n1)) then
        w3j=0.0d0
        return
      end if
!
      j1=n1
      j2=n2
      j3=n3
      m1=n4
      m2=n5
      m3=n6
!
!----------------------------------------------------------------------
!       use wig3j to calculate if j2 > 2
!----------------------------------------------------------------------
!
      if (n2.gt.2) then
        w3j=wig3j(j1,j2,j3,m1,m2,m3)
        return
      end if
!
!----------------------------------------------------------------------
!       permute variables if necessary to get m2 => 0 and j1 <= j3,
!       keep track of phases, and store variables
!----------------------------------------------------------------------
!
      j1=n1
      j2=n2
      j3=n3
      m1=n4
      m2=n5
      m3=n6
!
      phase=1.0d0
!
      if (mod(j1+j2+j3,2).eq.0) then
        parity=1.d0
      else
        parity=-1.d0
      end if
!
      if (m2.lt.0) then
        m1=-m1
        m2=-m2
        m3=-m3
        phase=parity
      else
        phase=1.d0
      end if
!
      if(j1.gt.j3) then
        k=j1
        j1=j3
        j3=k
        k=m1
        m1=m3
        m3=k
        phase=phase*parity
      end if
!
      if (mod(j1-m3,2).ne.0) then
        phase=-phase
      end if
!
!----------------------------------------------------------------------
!        calculate wigner 3-j symbols
!----------------------------------------------------------------------
!
      jdelta=j3-j1
      x=dble(j1+j3-1)
      y=dble(j1-m3)
      z=dble(j1+m3)
!
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
!
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
!
      return
      end



!----------------------------------------------------------------------
!
!                      =========================
!                                BINCOF
!                      =========================
!
!   Generates an array of binomial coefficients for use with function
!   This must be called before wig3j can calculate a 3-J symbol,
!   and is automatically called by wig3j if the array has not been
!   initialized.
!
!                         Designed and coded by
!                        Clemens C. J. Roothaan
!                        Department of Chemistry
!                         University of Chicago
!                        5735 South Ellis Avenue
!                        Chicago, Illinois 60637
!                           December 9, 1988
!
!                   modified by DJS 5/7/90, DEB 1/10/93
!----------------------------------------------------------------------

      subroutine bincof()
!
      use maxl
      use bincom
!
      integer i,j,ij
      double precision bncf0,temp
!
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
        bncf(ij)=1.0d0
        ij=ij+1
 20   continue
!
      return
      end

!----------------------------------------------------------------------
!                    =========================
!                         function WIG3J
!                    =========================
!
!     Calculates arbitrary 3-J symbol using binomial coefficients 
!     for J values up to an limit determined by the machine precision
!     (approx. 125 for IEEE format real*8 floating-point numbers -DEB)
!     
!     This replaces the function of the same name written by G. Moro 
!     for use with EPR spectral simulations.
!
!                         Designed and coded by
!                        Clemens C. J. Roothaan
!                        Department of Chemistry
!                         University of Chicago
!                        5735 South Ellis Avenue
!                        Chicago, Illinois 60637
!                           December 9, 1988
!
!                   modified by DJS 5/7/90, DEB 1/10/93
!----------------------------------------------------------------------

      function wig3j(j1,j2,j3,m1,m2,m3)
!
      use maxl
      use bincom
!
      double precision wig3j
      integer j1,j2,j3,m1,m2,m3
      integer i,j,k,l,m,n,p,q,z,zmin,zmax,bp,bnj,bmk
      double precision sum
      logical notset
      data notset / .true. /
!
!......................................................................
!
!
!
      if (notset) then
         call bincof
         notset = .false.
      end if
!
      i=j1+m1
      j=j1-m1
      k=j2+m2
      l=j2-m2
      m=j2+j3-j1
      n=j3+j1-j2
      p=j1+j2-j3
      q=j1+j2+j3 + 1
!
      if(q+1.gt.nb) then
        write (*,1000) j1,j2,j3,m1,m2,m3
        wig3j=0.0d0
        return
      endif
!
      bp=bncfx(p+1)
      bnj=bncfx(n+1)+n-j
      bmk=bncfx(m+1)+m-k
      zmin=max(0,j-n,k-m)
      zmax=min(p,j,k)
!
      sum=0.0d0
        do 30 z=zmin,zmax
        sum=-sum+bncf(bp+z)*bncf(bnj+z)*bncf(bmk+z)
 30     continue
!
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
!
      wig3j=sum
      return
!
!######################################################################
!
 1000   format(' wig3j called with (',6(i3,','),'):'/
     #         '   sum of L values exceeds limit')
      end
