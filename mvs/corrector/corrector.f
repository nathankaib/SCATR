c*************************************************************************
c                        CORRECTOR.F
c*************************************************************************
c This subroutine calculates the 3rd order symplectic corrector given in
c Wisdom (2006)
c
c             Input:
c                  nbod        ==>  number of massive bodies (int scalor)
c                  ntp         ==>  number of tp bodies (int scalor)
c                  mass        ==>  mass of bodies (real array)
c                  j2rp2,j4rp4 ==>  J2*radii_pl^2 and  J4*radii_pl^4
c                                     (real scalars)
c                  xh,yh,zh    ==>  massive part position 
c                                     (real arrays)
c                  xht,yht,zht ==>  test part position 
c                                     (real arrays)
c                  istat       ==>  status of the test paricles
c                                      (integer array)
c                                      istat(i) = 0 ==> active:  = 1 not
c                                    NOTE: it is really a 2d array but 
c                                          we only use the 1st row
c             Output:
c               xh,yh,zh    ==>  massive part position
c                                     (real arrays)
c               xht,yht,zht ==>  test part position 
c
c Author:  Nathan Kaib 
c Date:    6/2/09
c Last revision: 6/2/09

      subroutine corrector(nbod,npl,ntp,time,mass,j2rp2,j4rp4,xh,yh,zh,
     &     vxh,vyh,vzh,xht,yht,zht,vxht,vyht,vzht,istat,dt,tinc)

      include '../../swift.inc'

c...  Inputs: 
      integer nbod,ntp,istat(NTPMAX,NSTAT),bar,npl
      real*8 mass(NPLMAX),xh(NPLMAX),yh(NPLMAX),zh(NPLMAX)
      real*8 vxh(NPLMAX),vyh(NPLMAX),vzh(NPLMAX)

      real*8 xht(NTPMAX),yht(NTPMAX),zht(NTPMAX),j2rp2,j4rp4
      real*8 vxht(NTPMAX),vyht(NTPMAX),vzht(NTPMAX)

      real*8 time
c...  Outputs:
      
c...  Internals:
      integer i
      real*8 gamma,a1,b1,a2,b2,dt,adt,bdt,tinc
      real*8 xh2(NPLMAX),yh2(NPLMAX),zh2(NPLMAX)
      real*8 vxh2(NPLMAX),vyh2(NPLMAX),vzh2(NPLMAX)

c----
c...  Executable code
c     don't want to double correct pls
      do i=1,nbod
         xh2(i) = xh(i)
         yh2(i) = yh(i)
         zh2(i) = zh(i)
         vxh2(i) = vxh(i)
         vyh2(i) = vyh(i)
         vzh2(i) = vzh(i)
      enddo
c     shut off barycentric tps
      do i=1,ntp
         if (istat(i,NSTAT-2).lt.1.and.istat(i,1).eq.0) then
            istat(i,1) = -1
         endif
      enddo
      
c     correct heliocentric tps and pls
      gamma=sqrt(10.0)
      a1 = 3.0*gamma/10.0
      b1 = gamma/72.0
      adt = dt*a1/tinc
      bdt = dt*b1/tinc
      call correctterm_hel(nbod,npl,ntp,time,mass,j2rp2,j4rp4,xh,yh,zh,
     &     vxh,vyh,vzh,xht,yht,zht,vxht,vyht,vzht,istat,adt,bdt)
      
      a2 = gamma/5.0
      b2 = 0.0-gamma/24.0
      adt = dt*a2/tinc
      bdt = dt*b2/tinc
      call correctterm_hel(nbod,npl,ntp,time,mass,j2rp2,j4rp4,xh,yh,zh,
     &     vxh,vyh,vzh,xht,yht,zht,vxht,vyht,vzht,istat,adt,bdt)

c     turn on barycentric tps, shut off heliocentric tps
      do i=1,ntp
         if (istat(i,1).eq.-1) then
            istat(i,1) = 0
         endif
         if (istat(i,NSTAT-2).ge.1.and.istat(i,1).eq.0) then
            istat(i,1) = -1
         endif
      enddo

c     correct barycentric tps
      adt = dt*a1
      bdt = dt*b1
      call correctterm_bar(nbod,npl,ntp,time,mass,j2rp2,j4rp4,xh2,yh2,
     &     zh2,vxh2,vyh2,vzh2,xht,yht,zht,vxht,vyht,vzht,istat,adt,bdt)

      adt = dt*a2
      bdt = dt*b2
      call correctterm_bar(nbod,npl,ntp,time,mass,j2rp2,j4rp4,xh2,yh2,
     &     zh2,vxh2,vyh2,vzh2,xht,yht,zht,vxht,vyht,vzht,istat,adt,bdt)
      
c     turn all active tps back on
      do i=1,ntp
         if (istat(i,1).eq.-1) then
            istat(i,1) = 0
         endif
      enddo

      return
      end      ! getacch_tp

c---------------------------------------------------------------------




