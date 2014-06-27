c Version 1.4 1/2/95
c----------------------------------------------------------------------
c                    =========================
c                      common block ITERAT
c                    =========================
c
c   Iteration counter, residual norm, and related statistical quantities
c   for nonlinear least-squares procedure
c----------------------------------------------------------------------
c  
      integer iter,iwflag
      double precision fnorm,chisqr,rdchsq,ch2bnd,chnbnd,confid,delchi1,
     #                qfit,f2bnd,fnbnd,tbound,fnmin
      logical newitr,covarOK,xreset
      common /iterat/ fnorm,chisqr,rdchsq,ch2bnd,chnbnd,confid,delchi1,
     #                qfit,f2bnd,fnbnd,iter,iwflag,tbound,fnmin,
     #                newitr,covarOK,xreset







