c  subroutine ftest1 (subroutine reporting modified trace of a vector)
c  vec = vector to be tested 
c  siz = dimension of vector
c  txt = identification to be written
c  ltg = length in characters of txt array.
c
      subroutine ftest1(vec,siz,txt,ltg)
      integer siz,i,ltg
      character*(ltg) txt
      double precision vec(siz),ftest
      ftest=0.0d0
      do 10 i=1,siz
        ftest=ftest+vec(i)+dsqrt(vec(i)*vec(i))
10    continue
      write(*,*)txt,ftest
      return
      end
c
c  subroutine ftest2 (subroutine returning modified norm of an array)
c  arr = arr to be tested 
c  sizr,sizc = dimension of array
c  txt = idnetification to be written
c
      subroutine ftest2(arr,sizr,sizc,txt,lgt)
      integer sizr,sizc,i,j,lgt
      character*(lgt) txt
      double precision arr(sizr,sizc),ftest
      ftest=0.0d0
      do 10 i=1,sizc
        do 10 j=1,sizr
        ftest=ftest+arr(j,i)+dsqrt(arr(j,i)*arr(j,i))
10    continue
      write(*,*)txt,ftest
      return
      end
c
c  subroutine ftest3 (returns modified trace of a complex vector) 
c  vec = complex vector to be tested 
c  siz = dimension of vector
c  txt = idnetification to be written
c
      subroutine ftest3(vec,siz,txt,lgt)
      integer siz,i,lgt
      character*(lgt) txt
      complex*16 vec(siz),ftest
      ftest=(0.0d0,0.0d0)
      do 10 i=1,siz
        ftest=ftest+vec(i)+cdsqrt(vec(i)*vec(i))
10    continue
      write(*,*)txt,cdabs(ftest)
      return
      end
c
c  subroutine ftest4 (subroutine returning modified norm of an array)
c  arr = complex arr to be tested 
c  sizr,sizc = dimension of array
c  txt = idnetification to be written
c
      subroutine ftest4(arr,sizr,sizc,txt,lgt)
      integer sizr,sizc,i,j,lgt
      character*(lgt) txt
      complex*16 arr(sizr,sizc),ftest
      ftest=(0.0d0,0.0d0)
      do 10 i=1,sizc
        do 10 j=1,sizr
        ftest=ftest+arr(j,i)+cdsqrt(arr(j,i)*arr(j,i))
10    continue
      write(*,*)txt,cdabs(ftest)
      return
      end
c
