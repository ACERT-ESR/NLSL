! NLSL Version 1.4 10/10/94 
!----------------------------------------------------------------------
!                    ======================
!                      subroutine TWOFFT
!                    ======================
!     From Numerical Recipes by Press et al.
!----------------------------------------------------------------------
      subroutine twofft(data1,data2,fft1,fft2,n)
      implicit none
      integer j,n,n2
      double precision data1(n),data2(n)
      double complex fft1(n),fft2(n),h1,h2,c1,c2
      double precision rfft1(2*n)
!
      intrinsic dcmplx,dreal,dimag,dconjg
!
      double precision ZERO
      parameter(ZERO=0.0d0)
!
!......................................................................
!
      c1=dcmplx(0.5d0,ZERO)
      c2=dcmplx(ZERO,-0.5d0)
      do 11 j=1,n
        fft1(j)=dcmplx(data1(j),data2(j))
        rfft1(2*j-1)=data1(j)
        rfft1(2*j)=data2(j)
11    continue
      call four1(rfft1,n,1)
      do 21 j=1,n
        fft1(j)=dcmplx(rfft1(2*j-1),rfft1(2*j))
21    continue
      fft2(1)=dcmplx(dimag(fft1(1)),ZERO)
      fft1(1)=dcmplx(dreal(fft1(1)),ZERO)
      n2=n+2
      do 12 j=2,n/2+1
        h1=c1*(fft1(j)+dconjg(fft1(n2-j)))
        h2=c2*(fft1(j)-dconjg(fft1(n2-j)))
        fft1(j)=h1
        fft1(n2-j)=dconjg(h1)
        fft2(j)=h2
        fft2(n2-j)=dconjg(h2)
12    continue
      return
      end

!----------------------------------------------------------------------
!                    ======================
!                      subroutine REALFT
!                    ======================
!     From Numerical Recipes by Press et al.
!----------------------------------------------------------------------
      subroutine realft(data,n,isign)
      implicit none
      integer i,i1,i2,i3,i4,isign,n,n2p3
      double precision c1,c2,h1r,h1i,h2r,h2i,wrs,wis,wr,wi,wpr,wpi,
     #                 wtemp,theta
      double precision data(2*n)
!
      theta=6.28318530717959d0/2.0d0/dble(n)
      c1=0.5d0
      if (isign.eq.1) then
        c2=-0.5d0
        call four1(data,n,+1)
      else
        c2=0.5d0
        theta=-theta
      endif
      wpr=-2.0d0*dsin(0.5d0*theta)**2
      wpi=dsin(theta)
      wr=1.d0+wpr
      wi=wpi
      n2p3=2*n+3
      do 11 i=2,n/2+1
        i1=2*i-1
        i2=i1+1
        i3=n2p3-i2
        i4=i3+1
        wrs=sngl(wr)
        wis=sngl(wi)
        h1r=c1*(data(i1)+data(i3))
        h1i=c1*(data(i2)-data(i4))
        h2r=-c2*(data(i2)+data(i4))
        h2i=c2*(data(i1)-data(i3))
        data(i1)=h1r+wrs*h2r-wis*h2i
        data(i2)=h1i+wrs*h2i+wis*h2r
        data(i3)=h1r-wrs*h2r+wis*h2i
        data(i4)=-h1i+wrs*h2i+wis*h2r
        wtemp=wr
        wr=wr*wpr-wi*wpi+wr
        wi=wi*wpr+wtemp*wpi+wi
11    continue
      if (isign.eq.1) then
        h1r=data(1)
        data(1)=h1r+data(2)
        data(2)=h1r-data(2)
      else
        h1r=data(1)
        data(1)=c1*(h1r+data(2))
        data(2)=c1*(h1r-data(2))
        call four1(data,n,-1)
      endif
      return
      end

!----------------------------------------------------------------------
!                    ======================
!                      subroutine FOUR1
!                    ======================
!     From Numerical Recipes by Press et al.
!----------------------------------------------------------------------
      subroutine four1(data,nn,isign)
      implicit none
      integer i,istep,j,m,mmax,n,nn,isign
      double precision wr,wi,wpr,wpi,wtemp,tempi,tempr,theta
      double precision data(2*nn)
!
      n=2*nn
      j=1
      do 11 i=1,n,2
        if(j.gt.i)then
          tempr=data(j)
          tempi=data(j+1)
          data(j)=data(i)
          data(j+1)=data(i+1)
          data(i)=tempr
          data(i+1)=tempi
        endif
        m=n/2
1       if ((m.ge.2).and.(j.gt.m)) then
          j=j-m
          m=m/2
        go to 1
        endif
        j=j+m
11    continue
      mmax=2
2     if (n.gt.mmax) then
        istep=2*mmax
        theta=6.28318530717959d0/(isign*mmax)
        wpr=-2.d0*dsin(0.5d0*theta)**2
        wpi=dsin(theta)
        wr=1.d0
        wi=0.d0
        do 13 m=1,mmax,2
          do 12 i=m,n,istep
            j=i+mmax
            tempr=sngl(wr)*data(j)-sngl(wi)*data(j+1)
            tempi=sngl(wr)*data(j+1)+sngl(wi)*data(j)
            data(j)=data(i)-tempr
            data(j+1)=data(i+1)-tempi
            data(i)=data(i)+tempr
            data(i+1)=data(i+1)+tempi
12        continue
          wtemp=wr
          wr=wr*wpr-wi*wpi+wr
          wi=wi*wpr+wtemp*wpi+wi
13      continue
        mmax=istep
      go to 2
      endif
      return
      end
