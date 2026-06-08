c*************************************************************************
c                            STEP_SUN_PLUM.F
c*************************************************************************
c This subroutine integrates the sun in a plummer potential.  
c Does a KICK than a DRIFT than a KICK.
c ONLY DOES MASSIVE PARTICLES
c
c             Input:
c                 i1st          ==>  = 0 if first step; = 1 not (int scalar)
c                 nbod          ==>  number of massive bodies (int scalar)
c                 mass          ==>  mass of bodies (real array)
c                 j2rp2,j4rp4   ==>  J2*radii_pl^2 and  J4*radii_pl^4
c                                     (real scalars)
c                 xsun,ysun,zsun      ==>  initial position
c                                    (real arrays)
c                 vxsun,vysun,vzsun   ==>  initial velocity 
c                                    (real arrays)
c                 dt            ==>  time step
c             Output:
c                 xsun,ysun,zsun      ==>  final position 
c                                       (real arrays)
c                 vxsun,vysun,vzsun   ==>  final velocity 
c                                       (real arrays)
c
c Remarks: 
c Authors:  Nathan Kaib 
c Date:    9/27/07
c Last revision: 9/27/07

      subroutine step_sun_plum(clusm,pluma,xsun,ysun,zsun,vxsun,vysun,
     &     vzsun,dt)

c     input/output
      real*8 clusm,pluma,xsun,ysun,zsun,vxsun,vysun,vzsun,dt

c     internals
      real*8 rsun,acc,ax,ay,az

c----
c     Executable code
      
c     calculating acceleration
      rsun = sqrt(xsun*xsun+ysun*ysun+zsun*zsun)
      acc = -clusm*rsun/(rsun**2.+pluma**2.)**(3./2.)
      ax = acc*xsun/rsun
      ay = acc*ysun/rsun
      az = acc*zsun/rsun

c     applying kick for half-step
      vxsun = vxsun + 0.5*dt*ax
      vysun = vysun + 0.5*dt*ay
      vzsun = vzsun + 0.5*dt*az
      
c     drifting for full step
      xsun = xsun + dt*vxsun
      ysun = ysun + dt*vysun
      zsun = zsun + dt*vzsun

c     recalculating acceleration
      rsun = sqrt(xsun*xsun+ysun*ysun+zsun*zsun)
      acc = -clusm*rsun/(rsun**2.+pluma**2.)**(3./2.)
      ax = acc*xsun/rsun
      ay = acc*ysun/rsun
      az = acc*zsun/rsun

c     applying kick for half-step
      vxsun = vxsun + 0.5*dt*ax
      vysun = vysun + 0.5*dt*ay
      vzsun = vzsun + 0.5*dt*az

      return

      end
