c this function converts a mean anomaly into an eccentric anomaly

      real*8 function eccanom(capm,e)

      include '../swift.inc'
      real*8 k,E0,f0,f1,f2,f3,del1,del2,del3,nextone,diff,percent
      real*8 Eanom,cosang,acoscorr,Enoiseco,ang
      real*8 e,capm
      integer sign

      k=0.85
      sign=sin(capm)/abs(sin(capm))
      E0=capm+sign*k*e

      Eanom=E0
      diff=Eanom
      percent=diff/Eanom

 10   if(percent.lt.0.001) then
         goto 20
      endif

      f0=Eanom-e*sin(Eanom)-capm
      f1=1-e*cos(Eanom)
      f2=e*sin(Eanom)
      f3=e*cos(Eanom)

      del1=-f0/f1
      del2=-f0/(f1+.5*del1*f2)
      del3=-f0/(f1+.5*del2*f2+1./6.*del2**2.*f3)

      nextone=Eanom+del3
      diff=abs(Eanom-nextone)
      percent=diff/Eanom

      Eanom=nextone

      goto 10

 20   eccanom=Eanom
      return
      end
      end
