c*************************************************************************
c                            IO_WRITE_FLUX
c*************************************************************************
c writes out info on test particles when they exit or enter a 1 AU sphere
c around the sun
c
c             Input:
c                 i1stflux      ==>  integer indicating 1st flux file write
c                 time          ==>  current time (real scalar)
c                 nbod          ==>  number of massive bodies (int scalar)
c                 x,y,z    ==> current part position in helio coord 
c                                      (real arrays)
c                 vx,vy,vz ==> current velocity in helio coord 
c                                        (real arrays)
c                 iu              ==> unit number to write to
c                 fopenstat       ==>  The status flag for the open 
c                                      statements of the output files.  
c                                          (character*80)
c
c
c Remarks: Based on io_write_frame
c Authors:  Nathan Kaib 
c Date:    8/8/05
c Last revision: 

      subroutine io_write_flux(time,npl,mass,xh,yh,zh,vxh,vyh,vzh,
     &     i,x,y,z,vx,vy,vz,fopenstat)

      include '../swift.inc'
      include 'io.inc'

c...  Inputs:
      integer i,status,npl,j
      real*8 x,y,z
      real*8 vx,vy,vz
      real*8 time
      real*8 xbar,ybar,zbar,vxbar,vybar,vzbar,mtot
      real*8 xh(NPLMAX),yh(NPLMAX),zh(NPLMAX),mass(NPLMAX)
      real*8 vxh(NPLMAX),vyh(NPLMAX),vzh(NPLMAX)
      character*80 fopenstat

c...  Internals
      integer*2 id2,status2
      real*4 x4,y4,z4,t4,m4
      real*4 vx4,vy4,vz4
      integer ierr,iu
      integer flux1st    ! =0 first time through; =1 after
      data flux1st/0/
      save flux1st
c----
c...  Executable code 

c...  if first time through open file
      iu=70
      if(flux1st.eq.0) then
         call io_open(iu,'flux.dat',fopenstat,'UNFORMATTED',ierr)
         if(ierr.ne.0) then
           write(*,*) ' SWIFT ERROR: in io_write_flux: '
           write(*,*) '     Could not open binary output file:'
           call util_exit(1)
         endif
         flux1st = 1
      else
        call io_open(iu,'flux.dat','append','UNFORMATTED',ierr)
      endif
      
c      call io_write_hdr_r(iu,time,npl,1)
c
c      do j=2,npl
c         x4 = xh(j)
c         y4 = yh(j)
c         z4 = zh(j)
c         vx4 = vxh(j)
c         vy4 = vyh(j)
c         vz4 = vzh(j)
c         m4 = mass(j)
c
c         id2 = -1 - j
c         write(iu) id2,m4,x4,y4,z4,vx4,vy4,vz4
c      enddo

      mtot = mass(1)
      xbar = 0.0d0
      ybar = 0.0d0
      zbar = 0.0d0
      vxbar = 0.0d0
      vybar = 0.0d0
      vzbar = 0.0d0
      do j=2,npl
         xbar = xbar + xh(j) * mass(j)
         ybar = ybar + yh(j) * mass(j)
         zbar = zbar + zh(j) * mass(j)
         vxbar = vxbar + vxh(j) * mass(j)
         vybar = vybar + vyh(j) * mass(j)
         vzbar = vzbar + vzh(j) * mass(j)
         mtot = mtot + mass(j)
      enddo
      xbar = xbar / mtot
      ybar = ybar / mtot
      zbar = zbar / mtot
      vxbar = vxbar / mtot
      vybar = vybar / mtot
      vzbar = vzbar / mtot

      x4 = x - xbar
      y4 = y - ybar
      z4 = z - zbar
      vx4 = vx - vxbar
      vy4 = vy - vybar
      vz4 = vz - vzbar
      id2 = i
      m4 = 0.0

c     writing time
      call io_write_hdr_r(iu,time,0,1)

c     writing test particle data
      write(iu) id2,m4,x4,y4,z4,vx4,vy4,vz4

      close(unit = iu)
      return
      end
