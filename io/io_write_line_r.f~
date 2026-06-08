c*************************************************************************
c                            IO_WRITE_LINE_R
c*************************************************************************
c write out one line to real*4 binary file.
c
c      Input:
c            iu       ==> unit number to write to
C	     a        ==> semi-major axis or pericentric distance if a parabola
c                          (real scalar)
c            e        ==> eccentricity (real scalar)
C            inc      ==> inclination  (real scalar)
C            capom    ==> longitude of ascending node (real scalar)
C	     omega    ==> argument of perihelion (real scalar)
C	     capm     ==> mean anomoly(real scalar)
c
c Remarks: 
c Authors:  Hal Levison 
c Date:    2/22/94
c Last revision: 


c i changed this routine so it now writes the cartesian positions and 
c velocities to the binary file (nak 8/17/04)...and masses (nak 10/26/04)
      subroutine io_write_line_r(iu,id,nvis,mass,a,e,inc,capom,omega,
     &     capm)

      include '../swift.inc'
      include 'io.inc'

c...  Inputs: 
      integer iu,id,nvis
      real*8 mass
      real*8 a,e,inc,capom,omega,capm,x,y,z,vx,vy,vz

c...  Internals
      integer*2 id2
      real*4 a4,m4,e4,inc4,capom4,omega4,capm4,x4,y4,z4,vx4,vy4,vz4

c----
c...  Executable code 

      id2 = id

      a4 = a
      e4 = e
      inc4 = inc
      capom4 = capom
      capm4 = capm
      omega4 = omega
      x4 = x
      y4 = y
      z4 = z
      vx4 = vx
      vy4 = vy
      vz4 = vz
      m4 = mass

      write(iu) id2,nvis,m4,a4,e4,inc4,capom4,omega4,capm4

      return
      end      ! io_write_line_r
c--------------------------------------------------------------------------
