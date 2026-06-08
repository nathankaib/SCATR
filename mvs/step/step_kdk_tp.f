c*************************************************************************
c                            STEP_KDK_TP.F
c*************************************************************************
c This subroutine takes a step in helio coord.  
c Does a KICK than a DRIFT than a KICK.
c ONLY DOES TEST PARTICLES
c
c             Input:
c                 i1st           ==>  = 0 if first step; = 1 not (int scalar)
c                 nbod           ==>  number of massive bodies (int scalar)
c                 ntp            ==>  number of test bodies (int scalar)
c                 mass           ==>  mass of bodies (real array)
c                 j2rp2,j4rp4    ==>  J2*radii_pl^2 and  J4*radii_pl^4
c                                     (real scalars)
c                 xbeg,ybeg,zbeg ==>  massive part position at beginning of dt
c                                       (real arrays)
c                 xend,yend,zend ==>  massive part position at end of dt
c                                       (real arrays)
c                 xht,yht,zht    ==>  initial part position in helio coord 
c                                      (real arrays)
c                 vxht,vyht,vzht ==>  initial velocity in helio coord 
c                                        (real arrays)
c                 istat           ==>  status of the test paricles
c                                      (2d integer array)
c                                      istat(i,1) = 0 ==> active:  = 1 not
c                                      istat(i,2) = -1 ==> Danby did not work
c                 dt             ==>  time step
c                 longasc,argper,
c                 inc            ==> orbital elements of suns orbit in 
c                                    molecular cloud (real scalars)
c             Output:
c                 xht,yht,zht    ==>  final position in helio coord 
c                                       (real arrays)
c                 vxht,vyht,vzht ==>  final position in helio coord 
c                                       (real arrays)
c
c Remarks: Adopted from martin's nbwhnew.f program
c Authors:  Hal Levison 
c Date:    2/12/93
c Last revision: 2/24/94

      subroutine step_kdk_tp(i1st,time,nbod,npl,ntp,mass,j2rp2,j4rp4,
     &              xbeg,ybeg,zbeg,xend,yend,zend,xht,yht,zht,vxht,
     &              vyht,vzht,istat,dt,bar)	


      include '../../swift.inc'

c...  Inputs Only: 
      integer nbod,npl,ntp,i1st,bar,i
      real*8 mass(NPLMAX),dt,j2rp2,j4rp4,mcent
      real*8 xbeg(NPLMAX),ybeg(NPLMAX),zbeg(NPLMAX)
      real*8 xend(NPLMAX),yend(NPLMAX),zend(NPLMAX)
      real*8 crap1(NPLMAX),crap2(NPLMAX),crap3(NPLMAX)

c...  Inputs and Outputs:
      integer istat(NTPMAX,NSTAT)
      real*8 xht(NTPMAX),yht(NTPMAX),zht(NTPMAX)
      real*8 vxht(NTPMAX),vyht(NTPMAX),vzht(NTPMAX)

c...  Internals:
      real*8 dth 
      real*8 axht(NTPMAX),ayht(NTPMAX),azht(NTPMAX)
      real*8 time

      save axht,ayht,azht     ! Note this !!

c----
c...  Executable code 

      dth = 0.5d0*dt

      mcent = mass(1)
      if (bar.eq.0) then
         do i=2,npl
            mcent = mcent+mass(i)
         enddo
      endif

      if(i1st.eq.0) then 
c...      Get the accelerations in helio frame.
          call getacch_tp(nbod,npl,ntp,time,mass,j2rp2,j4rp4,xbeg,ybeg,
     &        zbeg,xht,yht,zht,istat,axht,ayht,azht,dt,mcent,bar)
          i1st = 1    ! turn this off
      endif

c...  Apply a heliocentric kick for a half dt 
      call kickvh_tp(ntp,vxht,vyht,vzht,axht,ayht,azht,istat,dth) 

c...  Take a drift forward full step
      call drift_tp(ntp,mcent,xht,yht,zht,vxht,vyht,vzht,dt,istat)	

c...  Get the accelerations in helio frame.
      call getacch_tp(nbod,npl,ntp,time,mass,j2rp2,j4rp4,xend,yend,zend,
     &     xht,yht,zht,istat,axht,ayht,azht,dt,mcent,bar)

c...  Apply a heliocentric kick for a half dt 
      call kickvh_tp(ntp,vxht,vyht,vzht,axht,ayht,azht,istat,dth)

      return
      end   ! step_kdk_tp
c---------------------------------------------------------------------

