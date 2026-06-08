c************************************************************************
c                              IO_PASS.F
c************************************************************************
c IO_PASS reads in the data for the next stellar passage. The last "planet"
c will normally be a passing star. 
c
c             Input:
c                 passfile      ==> File name to read from (character*80)
c                 lclose        ==> .true. --> discard particle if it gets 
c                                    too close to a planet. Read in that 
c                                    distance in io_init_pl
c                                      (logical*2 scalar)
c                 nbod          ==> Number of gravitating bodies
c                 tpass       ==> time to start next passage
c
c             Output:
c                 nbod          ==>  number of massive bodies (int scalar)
c                 mass          ==>  mass of bodies (real array)
c                 xh,yh,zh      ==>  initial position in Helio coord 
c                                    (real arrays)
c                 vxh,vyh,vzh   ==>  initial velocity in Helio coord 
c                                    (real arrays)
c                 rplsq         ==>  min distance^2 that a tp can get from pl
c                                    (real array)
c
c Remarks: 
c Authors:  Nathan Kaib
c Date:    10/25/04
c Last revision: 10/25/04  NAK

	subroutine io_pass(passfile,lclose,mass,xh,yh,zh,vxh,vyh,vzh,
     &     nbod,rplsq,t,tstop,tpass,index)

	include '../swift.inc'
	include 'io.inc'

c...    Input
	character*(*) passfile
	logical*2 lclose

c...    Output
	real*8 mass(NPLMAX),rplsq(NPLMAX)
	real*8 xh(NPLMAX),yh(NPLMAX),zh(NPLMAX)
	real*8 vxh(NPLMAX),vyh(NPLMAX),vzh(NPLMAX)
	real*8 xc,yc,zc,vxc,vyc,vzc,mtot
	real*8 t,tstop,tpass,tstar
	integer nbod

c...   Internal
	integer j,ierr,index,i
        real*8 rpl,masstest

c-----
c...  Executable code      

	write(*,*) 'Stellar passage file is ',passfile
	call io_open(7,passfile,'old','formatted',ierr)

	tstar=0
	
c       reading through all the passes that have already taken place
c       and then choosing the next one
	if(lclose) then
	   do while(tpass.ge.tstar)
 	      read(7,*,iostat=ierr) tstar,mass(index),rpl
	      read(7,*,iostat=ierr) xh(index),yh(index),zh(index)
	      read(7,*,iostat=ierr) vxh(index),vyh(index),vzh(index)
	      if(ierr.ne.0) goto 2
	   enddo
	else
	   do while(tpass.ge.tstar)
	      read(7,*,iostat=ierr) tstar,mass(index)
	      read(7,*,iostat=ierr) xh(index),yh(index),zh(index)
	      read(7,*,iostat=ierr) vxh(index),vyh(index),vzh(index)
	      if(ierr.ne.0) goto 2
	   enddo
	endif

c       converting star's barycentric info to heliocentric
	xc = 0.
	yc = 0.
	zc = 0.
	vxc = 0.
	vyc = 0.
	vzc = 0.
	mtot = 0.

c       calculating location and velocity of center-of-mass
	do i=1,5
	   xc = xc + xh(i)*mass(i)
	   yc = yc + yh(i)*mass(i)
	   zc = zc + zh(i)*mass(i)
	   vxc = vxc + vxh(i)*mass(i)
	   vyc = vyc + vyh(i)*mass(i)
	   vzc = vzc + vzh(i)*mass(i)
	   mtot = mtot + mass(i)
	enddo
	xc = xc/mtot
	yc = yc/mtot
	zc = zc/mtot
	vxc = vxc/mtot
	vyc = vyc/mtot
	vzc = vzc/mtot

c       putting star's pos and vel in heliocentric coords
	xh(index) = xh(index) + xc
	yh(index) = yh(index) + yc
	zh(index) = zh(index) + zc
	vxh(index) = vxh(index) + vxc
	vyh(index) = vyh(index) + vyc
	vzh(index) = vzh(index) + vzc

 2	if(tstar.lt.tpass) then
	   tpass=tstop
	else
	   tpass=tstar
	endif

c       displaying the data for the next stellar passage
	if(lclose) then
	   write(*,*) mass(index),rpl,index
	   rplsq(index) = rpl*rpl
	else
	   write(*,*) mass(index),index
	endif

	write(*,*) xh(index),yh(index),zh(index)
	write(*,*) vxh(index),vyh(index),vzh(index)

	close(unit = 7)
	
	return
	end     ! io_pass.f
c--------------------------------------------------------------------------
