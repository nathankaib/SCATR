c converts binary file to ascii file
c this modified version of the program creates an ascii
c output file for each planet and test particle in the 
c simulation.  planet files are named pl[particle #].out
c and test particles files are named tp[particle #].out

	include 'swift.inc'

	real*8 xht(NTPMAX),yht(NTPMAX),zht(NTPMAX)
	real*8 vxht(NTPMAX),vyht(NTPMAX),vzht(NTPMAX)

	real*8 mass(NPLMAX),dr,peri
	real*8 xh(NPLMAX),yh(NPLMAX),zh(NPLMAX)
	real*8 vxh(NPLMAX),vyh(NPLMAX),vzh(NPLMAX)

	integer istat(NTPMAX,NSTAT)
        real*8 rstat(NTPMAX,NSTATR)
	integer nbod,ntp,ierr,ifol,istep
	integer iflgchk,iu,nleft,i,id,j,labelint
        integer io_read_hdr,io_read_line
        integer io_read_hdr_r,io_read_line_r

	real*8 t0,tstop,dt,dtout,dtdump,t
	real*8 time

	real*8 rmin,rmax,rmaxu,qmin,rplsq(NPLMAX)
        logical*2 lclose
        real*8 a,e,inc,capom,omega,capm,j2rp2,j4rp4,m
        real*8 elh,elk,elp,elq,apo

	character*80 outfile,inparfile,inplfile,intpfile,fopenstat

c       these extra varibles are to generate the different filenames
	character filenamepl1*17,filenamepl2*22,label*5,label2*5
	character filenametp1*17,filenametp2*22
	integer zeros,n

	real*8 x,y,z,vx,vy,vz
	
c Get data for the run and the test particles
	write(*,*) 'Enter name of parameter data file : '
	read(*,999) inparfile
	call io_init_param(inparfile,t0,tstop,dt,dtout,dtdump,
     &         iflgchk,rmin,rmax,rmaxu,qmin,lclose,outfile,fopenstat)

c Prompt and read name of planet data file
	write(*,*) ' '
	write(*,*) 'Enter name of planet data file : '
	read(*,999) inplfile
999 	format(a)
	call io_init_pl(inplfile,lclose,iflgchk,nbod,mass,xh,yh,zh,
     &       vxh,vyh,vzh,rplsq,j2rp2,j4rp4)

c Get data for the run and the test particles
	write(*,*) 'Enter name of test particle data file : '
	read(*,999) intpfile
	call io_init_tp(intpfile,ntp,xht,yht,zht,vxht,vyht,
     &               vzht,istat,rstat)

        iu = 20

        dr = 180.0/PI

        if(btest(iflgchk,0)) then
           write(*,*) ' Reading an integer*2 binary file '
        else if(btest(iflgchk,1)) then
           write(*,*) ' Reading an real*4 binary file '
        else
           write(*,*) ' ERROR: no binary file format specified '
           write(*,*) '        in param file '
           stop
        endif

        do j=-nbod,-2
	   ifol=j

	   write(*,*) ' Following particle ',ifol

c       creating label for output file
	   labelint=-j
	   write (label,300)labelint
 300	   format(I5)
	   zeros=len(label)-int(aint(log10(dble(-j))))-1
	   do n=1,zeros
	      label(n:n)='0'
	   enddo

	   filenamepl1='pl'//label//'.clo'

	   open(unit=iu, file='clo.dat', status='old',form='unformatted')  
	   open(unit=7,file=filenamepl1)
	   
	   write(*,*) '1 2 3 4 5   6     7     8    9    10  11 12 13 14 15 16'
	   write(*,*) 't,m,a,e,inc,capom,omega,capm,peri,apo, x, y, z,vx,vy,vz'
	   
 1	   continue
	   if(btest(iflgchk,0))  then ! bit 0 is set
	      ierr = io_read_hdr(iu,t,nbod,nleft) 
	   else
	      ierr = io_read_hdr_r(iu,t,nbod,nleft) 
	   endif

	   if(ierr.ne.0) goto 2

	   istep = 0
	   do while (ierr.eq.0)
	      if(btest(iflgchk,0))  then ! bit 0 is set
		 ierr = io_read_line(iu,id,a,e,inc,capom,omega,capm) 
	      else
		 ierr = io_read_line_r(iu,id,t,a,e,inc,capom,omega,capm,x,
     &		      y,z,vx,vy,vz)
	      endif
	      if(ierr.ne.0) goto 2
	      if(id.eq.ifol) then
		 istep = 1
		 inc = inc*dr
		 capom = capom*dr
		 omega = omega*dr
		 capm = capm*dr
		 peri = a*(1.0d0-e)
		 apo = a*(1.0d0+e)
c putting upper limit on apo
		 if(apo.gt.1000000.) apo = 999999.
		 write(7,1000) t,m,a,e,inc,capom,omega,capm,peri,apo,x,y,z,
     &		      vx,vy,vz
 1000		 format(16(1x,e15.7))
		 write(8,1001) t,elh,elk,elp,elq
 1001		 format(1x,e13.5,4(1x,f7.5))

	      endif
	   enddo
 2	   continue

c       file access is sequential, so i'm reseting the read position
	   rewind(iu)
	enddo


c---------------------------------------------------------------------
c  getting info for test particles
	do j=1,ntp
	   ifol=j

	   write(*,*) ' Following particle ',ifol

c       creating label for output file
	   labelint=j
	   write (label,300)labelint
	   zeros=len(label)-int(aint(log10(dble(j))))-1
	   do n=1,zeros
	      label(n:n)='0'
	   enddo
	   
	   filenametp1='tp'//label//'.clo'
	   filenametp2='tp_hkpq'//label//'.clo'

	   open(unit=iu, file='clo.dat', status='old',form='unformatted')
	   open(unit=7,file=filenametp1)
	   open(unit=8,file=filenametp2)

	   write(*,*) '1 2 3 4 5   6     7     8    9    10  11 12 13 14 15 16'
	   write(*,*) 't,m,a,e,inc,capom,omega,capm,peri,apo, x ,y, z,vx,vy,vz'
	   
 3	   continue
	   ierr=0
	   istep=0

	   do while (ierr.eq.0)
	      ierr = io_read_line_r(iu,id,t,a,e,inc,capom,omega,capm,x,
     &		      y,z,vx,vy,vz)
	      if(ierr.ne.0) goto 4
	      if(id.eq.ifol) then
		 istep = 1
		 elh = e*cos(omega+capom)
		 elk = e*sin(omega+capom)
		 elp = sin(inc/2.0)*cos(capom)
		 elq = sin(inc/2.0)*sin(capom)
		 inc = inc*dr
		 capom = capom*dr
		 omega = omega*dr
		 capm = capm*dr
		 peri = a*(1.0d0-e)
c putting upper limit on apo
 2000		 format(16(1x,e15.7))
		 apo = a*(1.0d0+e)
		 if(apo.gt.1000000.) apo = 999999.
		 write(7,2000) t,m,a,e,inc,capom,omega,capm,peri,apo,x,y,z,
     &		      vx,vy,vz
		 write(8,2001) t,elh,elk,elp,elq
 2001		 format(1x,e13.5,4(1x,f7.5))
	      endif
	   enddo
	   if(istep.eq.0) goto 4 ! did not find particle this time step
	   
	   goto 3

 4	   continue

	   write(*,*) ' Tmax = ',time
c       file access is sequential so i'm resetting the read position
	   rewind(iu)
	enddo

        stop
        end
