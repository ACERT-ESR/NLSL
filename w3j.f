c	Version 1.5   5/2/94
c***********************************************************************
c
c                       ===================
c                          function W3J
c                       ===================
c
c       This double precision function will calculate the values of
c       the Wigner 3-J symbols used in the Stochastic Liouville matrix
c       calculation, i.e.
c
c                         ( J1  J2  J3 )
c                         ( M1  M2  M3 )
c
c       For J2 <= 2 and arbitrary J1 and J3, this routine explicitly
c       evaluates formulae given in Table 2 of "Angular Momentum in 
c       Quantum Mechanics" by A. R. Edmonds, revised 2nd printing, 
c       Princeton Univ. Press, 1968.
c
c       For J2 > 2 and J1,J3 < 48 this functions calls a modified version 
c       function J3J originally by  G. Moro.
c
c       Written by D.E. Budil 1/8/92 based on programs by D.J. Schneider 
c       and G. Moro 
c
c**********************************************************************
c
        double precision function w3j(n1,n2,n3,n4,n5,n6)
c
        implicit integer (i-n)
        implicit double precision (a-h,o-z)
c
c#######################################################################
c
c-----------------------------------------------------------------------
c       check triangle conditions, etc.
c-----------------------------------------------------------------------
c
        if ((n1.lt.0).or.(n2.lt.0).or.(n3.lt.0).or.
     #     ((n2.gt.2).and.((n1.gt.48).or.(n3.gt.48))))
     #       stop 'W3J called with illegal parameters'
c
        if ((abs(n4).gt.n1).or.(abs(n5).gt.n2).or.(abs(n6).gt.n3).or.
     #     ((n4+n5+n6).ne.0).or.
     #     ((n1+n2).lt.n3).or.((n1+n3).lt.n2).or.((n2+n3).lt.n1)) then
           w3j=0.0d0
           return
        end if
c
c-----------------------------------------------------------------------
c       use wig3j to calculate if j2 > 2
c-----------------------------------------------------------------------
c
        if (n2.gt.2) then
           w3j=wig3j(n1,n2,n3,n4,n5,n6)
           return
        end if
c
c-----------------------------------------------------------------------
c       permute variables if necessary to get m2 => 0 and j1 <= j3,
c       keep track of phases, and store variables
c-----------------------------------------------------------------------
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
        if (mod(j1-m3,2).ne.0) phase=-phase
c
c-----------------------------------------------------------------------
c        calculate wigner 3-j symbols
c-----------------------------------------------------------------------
c
        jdelta=j3-j1
        x=dble(2*j1-1+jdelta)
        y=j1-m3
        z=j1+m3
c
        if (j2.eq.0) then
           w3j=phase/sqrt(2.d0*j1+1.d0)
        else if (j2.eq.2) then
           ztemp2=x*(x+1.d0)*(x+2.d0)*(x+3.d0)*(x+4.d0)
           if (m2.eq.0) then
              if (jdelta.eq.0) then
                 ztemp1=2.d0*(3.d0*m3*m3-j1*(j1+1.d0))
                 w3j=phase*ztemp1/sqrt(ztemp2)
              else if (jdelta.eq.1) then
                 ztemp1=6.d0*(z+1.d0)*(y+1.d0)
                 w3j=-phase*2.d0*m3*sqrt(ztemp1/ztemp2)
              else
                 ztemp1=6.d0*(z+2.d0)*(z+1.d0)
     #                *(y+2.d0)*(y+1.d0)
                 w3j=phase*sqrt(ztemp1/ztemp2)
              end if
           else if (m2.eq.1) then
              if (jdelta.eq.0) then
                 ztemp1=6.d0*(z+1.d0)*y
                 w3j=phase*(2.d0*m3+1.d0)*sqrt(ztemp1/ztemp2)
              else if (jdelta.eq.1) then
                 ztemp1=y*(y+1.d0)
                 w3j=-phase*(2.d0*j1+4.d0*m3+4.d0)*sqrt(ztemp1/ztemp2)
              else
                 ztemp1=(z+2.d0)*y*(y+1.d0)*(y+2.d0)
                 w3j=phase*2.d0*sqrt(ztemp1/ztemp2)
              end if
           else
              if (jdelta.eq.0) then
                 ztemp1=6.d0*(y-1.d0)*y*(z+1.d0)*(z+2.d0)
                 w3j=phase*sqrt(ztemp1/ztemp2)
              else if (jdelta.eq.1) then
                 ztemp1=(y-1.d0)*y*(y+1.d0)*(z+2.d0)
                 w3j=-phase*2.d0*sqrt(ztemp1/ztemp2)
              else
                 ztemp1=(y-1.d0)*y*(y+1.d0)*(y+2.d0)
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
c
c
c******************************************************************
c
c  This function calculates the Wigner 3-J symbols for the program
c  eprll for the special cases handled above.  This function is an
c  extensive revision of G. Moro's routine of the same name.  The 
c  modifications mainly have to do with streamlining the function 
c  for use with matrll, which assumes certain relationships between
c  the L and M values of the 3-J symbol.
c
c  Calling this program in any other applications may lead to serious
c  errors.
c
c*******************************************************************
c
c
      double precision function wig3j(ni,j,k,l,m,n)
      implicit real*8(a-h,o-z)
      implicit integer*4(i-n)
      dimension h(101),jj(101)
      logical notset
c
      data notset /.true./
c
c -- function definitions
c
      intptf(q)=q+q+dsign(.1d0,q)
      iparf(ij)=1-mod(ij,4)
c
c......................................................................
c
c -- Check whether array of factorials has already been computed
c 
      if (notset) then
         h(1)=1.0d0
         jj(1)=0
         x=0.d0
c
c     --- Build array of factorials and keep track of magnitude
c
         do 10 i=2,101
            x=x+1.0d0
            h(i)=h(i-1)*x
            jj(i)=jj(i-1)
c
 9          if(h(i).gt.10.0d0) then
               h(i)=0.01d0*h(i)
               jj(i)=jj(i)+2
               go to 9
            end if
 10      continue
         notset=.false.
      end if
c
      a=dble(ni)
      b=dble(j)
      c=dble(k)
      xx=dble(l)
      yy=dble(m)
      zz=dble(n)
      zz=-zz
c
      k1=intptf(a)
      k2=intptf(b)
      k3=intptf(c)
      k4=intptf(xx)
      k5=intptf(yy)
      k6=intptf(zz)
c
      if((k4+k5-k6).eq.0) then
         m1=k1+k2-k3
         m2=k2+k3-k1
         m3=k3+k1-k2
         m4=k1+k4
         m5=k1-k4
         m6=k2+k5
         m7=k2-k5
         m8=k3+k6
         m9=k3-k6
         m10=k1+k2+k3+2
c
         if ( m1.ge.0 .and. m2.ge.0 .and. m3.ge.0 .and. m4.ge.0 .and. 
     #        m5.ge.0 .and. m6.ge.0 .and. m7.ge.0 .and. m8.ge.0 .and.
     #        m9.ge.0 .and. mod(m4,2).eq.0 .and. mod(m10,2).eq.0 ) then
c
            y=k3+1
            m1=m1/2+1
            m2=m2/2+1
            m3=m3/2+1
            m4=m4/2+1
            m5=m5/2+1
            m6=m6/2+1
            m7=m7/2+1
            m8=m8/2+1
            m9=m9/2+1
            m10=m10/2+1
c
            y=dsqrt( y * h(m1) * h(m2) * h(m3) * h(m4) * h(m5) *
     #               h(m6) * h(m7) * h(m8) * h(m9) / h(m10) )
            iy=( jj(m1) + jj(m2) + jj(m3) + jj(m4) + jj(m5) + 
     #           jj(m6) + jj(m7) + jj(m8) + jj(m9) - jj(m10) ) / 2
c
            n4=m1
            if(n4.gt.m5) n4=m5
            if(n4.gt.m6)n4=m6
            n4=n4-1
c
            m2=k2-k3-k4
            m3=k1+k5-k3
c
            n5=0
            if(n5.lt.m2) n5=m2
            if(n5.lt.m3) n5=m3
            n5par=iparf(n5)
            n5=n5/2
            z=0.0d0
c
c --  Calculate w3j symbol as a sum of products of factorials
c
 11         if (n5.le.n4) then
               mm1=m1-n5
               mm2=m5-n5
               mm3=m6-n5
               mm4=n5-(m2/2)+1
               mm5=n5-(m3/2)+1
c
               x=1.d0/(h(mm1)*h(mm2)*h(mm3)*h(mm4)*h(mm5)*h(n5+1))
               ix=-jj(mm1)-jj(mm2)-jj(mm3)-jj(mm4)-jj(mm5)-jj(n5+1)
c
c         ---  Adjust magnitude of result
c
 12            if(ix+iy.lt.0) then
                  x=0.1d0*x
                  ix=ix+1
                  go to 12
               end if
c
 13            if (ix+iy.gt.0) then
                  x=10.0d0*x
                  ix=ix-1
                  go to 13
               end if 
c
               if (n5par.lt.0) x=-x
               z=z+x
               n5par=-n5par
               n5=n5+1
               go to 11
            end if
            clebsh=z*y
c
         else
            clebsh=0.0d0
         end if
c
      else
         clebsh=0.0d0
      end if
c
      js=abs(k1-k2+k6)
      jspar=iparf(js)
      wig3j=jspar*clebsh/dsqrt(k3+1.0d0)
      return
      end
