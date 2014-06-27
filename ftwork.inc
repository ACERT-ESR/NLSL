c NLSL Version 1.3
c----------------------------------------------------------------------
c                    =========================
c                      common block FTWORK
c                    =========================
c
c Defines common block for working arrays used for correlation and
c convolution function calculations by Fourier transform methods
c
c    tmpclc     Temporary storage array for calculated spectrum 
c               used by subroutine sshift (and subroutine correl)
c               and by eprls (and subroutine convlv)
c
c    tmpdat     Temporary storage array for data
c               used by subroutine sshift (and subroutine correl)
c               and by eprls (and subroutine convlv)
c
c    work,work2 Temporary work array used by foutines correl and convlv
c
c----------------------------------------------------------------------

      double precision work,work2,tmpclc,tmpdat
      common /ftwork/ work(MXPT),work2(MXPT),tmpclc(MXPT),tmpdat(MXPT)


