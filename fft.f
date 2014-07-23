c NLSPMC Version 1.0 2/5/99
      subroutine fft(data,n)
c*********************************************************************
c     This subroutine performs complex fast Fourier transform with
c     double precision.  It replaces dffcf routine in imsl library.
c
c            data  :  array of length nn
c            n     :  integer which is power of 2
c
c*********************************************************************
c
      implicit none
c
      integer n,m,nd2,nm1,i,j,k,l,le,le1,ip
      double precision pi,arg
      complex*16 data(n),u,w,temp
c
c#####################################################################
c
      pi=4.0D0*datan(1.0D0)
      m=int(dlog(dfloat(n))/dlog(2.0D0)+1.0D-1)
      nd2=n/2
      nm1=n-1
      j=1
      do 10 i=1,nm1
        if (i.lt.j) then
          temp=data(i)
          data(i)=data(j)
          data(j)=temp
        end if
        k=nd2
c
  15    if (k.lt.j) then
          j=j-k
          k=k/2
          go to 15
        else
          j=j+k
        end if
c
  10  continue
c
      do 20 l=1,m
        le=2**l
        le1=le/2
        u=(1.0D0,0.0D0)
        arg=pi/le1
        w=dcmplx(dcos(arg),-dsin(arg))
        do 30 j=1,le1
          do 40 i=j,n,le
            ip=i+le1
            temp=data(ip)*u
            data(ip)=data(i)-temp
            data(i)=data(i)+temp
  40      continue
          u=u*w
  30    continue
  20  continue
c
      return
      end
