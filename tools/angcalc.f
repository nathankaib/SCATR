c this function converts an eccentric anomaly into an angular position
c in an orbit
      real*8 function angcalc(Eanom,e)


      include '../swift.inc'
      real*8 Eanom,cosang,e,ang,acoscorr,Enoiseco

      Eanom=abs(Eanom)
      cosang=(cos(Eanom)-e)/(1-e*cos(Eanom))
      ang=acos(cosang)

      acoscorr=aint(Eanom/PI)*2.*PI
      Enoiseco=aint(Eanom/(2.*PI))*4.*PI
      acoscorr=acoscorr-Enoiseco

      ang=acoscorr-ang
      ang=(ang*ang)**.5

      angcalc=ang

      return
      end
