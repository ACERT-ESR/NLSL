      subroutine mapxxx(xxx,nvar,idir)
c
c subroutine to map infinite range variables handled by levenburg 
c routine into finite range physical parameters.
c This routine is only called from lmnls.f.
c
c    xxx: 	variables passed by LM routine
c    bvar: 	Number of parameters in xxx
c    idir:	direction to transform, 1=to physical values, else = LM vars.
c 
      integer i,nvar,idir
      double precision xxx(*)
c
      include 'limits.inc'
      include 'parms.inc'
c
      do 10 i=1,nvar
c check if this variable has min and max set:
        if(abs(prmax(i))+abs(prmin(i)).gt.1D-7) then
          if(idir.eq.1) then
c map to physical values for spectral calculation
            xxx(i)=prmin(i)+(prmax(i)-prmin(i))*
     #  	(1.D0-1.D0/(1.D0+exp(4.*(xxx(i)-
     #		(prmin(i)+prmax(i))/2.D0)/(prmax(i)-prmin(i)))))
          else
c map to infinite range for LM routine
            if(xxx(i).ge.prmax(i).or.xxx(i).le.prmin(i)) then
c values should be within limits, force them.
              write(*,*)'error parameter limits exceeded',
     #			i,xxx(i),prmax(i),prmin(i)
              if(xxx(i).ge.prmax(i)) xxx(i)=prmax(i)* 0.99999D0
              if(xxx(i).le.prmax(i)) xxx(i)=prmin(i)* 1.00001D0
            end if
c ok to do the mapping:
              xxx(i)=(prmin(i)+prmax(i))/2.+(prmax(i)-prmin(i))/4.*
     #          log(1./(1.-(xxx(i)-prmin(i))/(prmax(i)-prmin(i)))-1.)
          end if
        end if
 10   continue
      return
      end
