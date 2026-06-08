c**********************************************************************
c		      SWIFT_RMVS3.F
c**********************************************************************
c
c                 INCLUDES CLOSE ENCOUNTERS
c                 To run, need 3 input files. The code prompts for
c                 the file names, but examples are :
c
c                   parameter file like       param.in
c		    planet file like          pl.in
c                   test particle file like   tp.in
c                   cluster environment file  cluster.in
c                   stellar passage file like pass.in
c
c Authors:  Hal Levison \& Martin Duncan, (Nathan Kaib)
c Date:    8/25/94
c Last revision: 8/19/08 NAK

     
	include 'swift.inc'

	real*8 xht(NTPMAX),yht(NTPMAX),zht(NTPMAX),idt(NTPMAX)
	real*8 vxht(NTPMAX),vyht(NTPMAX),vzht(NTPMAX)
	real*8 xht2(NTPMAX),yht2(NTPMAX),zht2(NTPMAX)
	real*8 vxht2(NTPMAX),vyht2(NTPMAX),vzht2(NTPMAX)
	real*8 tpfrac(NTPMAX)

	real*8 mass(NPLMAX),j2rp2,j4rp4
	real*8 xh(NPLMAX),yh(NPLMAX),zh(NPLMAX)
	real*8 vxh(NPLMAX),vyh(NPLMAX),vzh(NPLMAX)
	real*8 xh2(NPLMAX),yh2(NPLMAX),zh2(NPLMAX)
	real*8 vxh2(NPLMAX),vyh2(NPLMAX),vzh2(NPLMAX)

	integer istat(NTPMAX,NSTAT),i1st
	integer nbod,ntp,nleft,i,j,npl
	integer iflgchk,iub,iuj,iud,iue
	integer clflg,nclose,iclose(NTPMAX)
        real*8 rstat(NTPMAX,NSTATR)

	real*8 t0,tstop,dt,dtout,dtdump
	real*8 tout,tdump,tfrac,eoff
	real*8 t

	real*8 rmin,rmax,rmaxu,qmin,rplsq(NPLMAX)
        logical*2 lclose 

	real*8 rcrit,tinc
	real*8 qrec,fluxdist
	integer xi,yi,ftest,ftestall

	character*80 outfile,inparfile,inplfile,intpfile,stfile
	character*80 fopenstat,frame

	integer nbod0,icenc
	character*2 du

	real*8 r2hill(NPLMAX)
	common/c_rhill/r2hill

	real*8 rad
	du='Stellar encouter'
c-----
c...    Executable code 

c...    print version number
        call util_version

c Get data for the run and the test particles
	write(*,*) 'Enter name of parameter data file : '
	read(*,999) inparfile
	call io_init_param(inparfile,t0,tstop,dt,tinc,dtout,
     &      dtdump,iflgchk,rmin,rmax,rmaxu,qmin,lclose,
     &      rcrit,qrec,fluxdist,outfile,frame,fopenstat,ftest)
	ftestall=ftest

c Prompt and read name of planet data file
	write(*,*) ' '
	write(*,*) 'Enter name of planet data file : '
	read(*,999) inplfile
999 	format(a)
	call io_init_pl(inplfile,lclose,iflgchk,nbod,npl,mass,xh,
     &       yh,zh,vxh,vyh,vzh,rplsq,j2rp2,j4rp4,ftest)
	if (ftest.ne.0) then 
	   ftestall=ftest
	endif
	write(*,*) ftest
	
c Get data for the run and the test particles
	write(*,*) 'Enter name of test particle data file : '
	read(*,999) intpfile
	call io_init_tp(intpfile,ntp,idt,tpfrac,xht,yht,zht,vxht,
     &               vyht,vzht,istat,rstat,ftest,frame)
	if (ftest.ne.0) then 
	   ftestall=ftest
	endif
	
	if (ftestall.ne.0) then
	   ftestall = 0
	   call io_init_param('backupparam.dat',t0,tstop,dt,tinc,dtout,
     &      dtdump,iflgchk,rmin,rmax,rmaxu,qmin,lclose,
     &      rcrit,qrec,fluxdist,outfile,frame,fopenstat,ftest)
	   ftestall=ftest
	   call io_init_pl('backuppl.dat',lclose,iflgchk,nbod,npl,mass,
     &       xh,yh,zh,vxh,vyh,vzh,rplsq,j2rp2,j4rp4,ftest)
	   if (ftest.ne.0) then 
	      ftestall=ftest
	   endif
	call io_init_tp('backuptp.dat',ntp,idt,tpfrac,xht,yht,zht,vxht,
     &               vyht,vzht,istat,rstat,ftest,frame)
	   if (ftest.ne.0) then 
	      ftestall=ftest
	   endif
	   if (ftestall.ne.0) then
	      write(*,*) 'Backup files messed up!'
	      stop
	   endif
	endif
	
c Initialize initial time and times for first output and first dump
	t = t0
	tout = t0 + dtout
	tdump = t0 + dtdump

        iub = 20
        iuj = 30
        iud = 40
        iue = 60

c...    Do the initial io write
        if(btest(iflgchk,0))  then ! bit 0 is set
           call io_write_frame(t0,nbod,ntp,mass,xh,yh,zh,vxh,vyh,
     &         vzh,xht,yht,zht,vxht,vyht,vzht,istat,outfile,iub,
     &         fopenstat)
        endif
        if(btest(iflgchk,1))  then ! bit 1 is set
           call io_write_frame_r(t0,nbod,npl,ntp,mass,xh,yh,zh,vxh,vyh,
     &         vzh,idt,tpfrac,xht,yht,zht,vxht,vyht,vzht,istat,
     &         outfile,iub,frame,qrec,fopenstat)
        endif

        if(btest(iflgchk,2))  then    ! bit 2 is set
           eoff = 0.0d0
	   nclose = ntp
           call anal_energy_write(t0,nbod,mass,j2rp2,j4rp4,xh,yh,zh,vxh,
     &          vyh,vzh,iue,fopenstat,eoff,nclose)
        endif
        if(btest(iflgchk,3))  then    ! bit 3 is set
           call anal_jacobi_write(t0,nbod,ntp,mass,xh,yh,zh,vxh,
     &        vyh,vzh,xht,yht,zht,vxht,vyht,vzht,istat,2,iuj,fopenstat)
        endif

c...    must initize discard io routine
        if(btest(iflgchk,4))  then ! bit 4 is set
           call io_discard_write(0,t,nbod,ntp,xh,yh,zh,vxh,vyh,
     &          vzh,idt,tpfrac,xht,yht,zht,vxht,vyht,vzht,istat,rstat,
     &          iud,'discard.out',fopenstat,nleft)
        endif

        nleft = ntp
	i1st = 0
	icenc = 0
	nbod0=npl

c       apply inverse corrector map to numerical hamiltonian
	call invcorrector(nbod,npl,ntp,t,mass,j2rp2,j4rp4,xh,yh,zh,
     &     vxh,vyh,vzh,xht,yht,zht,vxht,vyht,vzht,istat,dt,tinc)
	
c***************here's the big loop *************************************
        write(*,*) ' ************** MAIN LOOP ****************** '
	  do while ( (t .le. tstop) .and.(nleft.gt.0))
c            checking to see if tp's enter heliocentric zone
	     call rmvs3_hel(nbod,npl,ntp,t,mass,j2rp2,j4rp4,xh,yh,zh,
     &         vxh,vyh,vzh,xht,yht,zht,vxht,vyht,vzht,istat,dt,tinc,
     &         iclose,nclose,rcrit)

 	     call rmvs3_step(i1st,t,nbod,npl,ntp,
     &         mass,j2rp2,j4rp4,xh,yh,zh,vxh,vyh,vzh,idt,xht,yht,zht,
     &         vxht,vyht,vzht,istat,rstat,dt,tinc,fopenstat,nclose,
     &         iclose,lclose,rmin,rmax,rmaxu,qmin,rplsq,rcrit,fluxdist)

	     t = t + dt
	     
	     if(btest(iflgchk,4))  then	! bit 4 is set
                call discard(t,dt,nbod,npl,ntp,mass,xh,yh,zh,vxh,vyh,
     &               vzh,idt,xht,yht,zht,vxht,vyht,vzht,rmin,rmax,rmaxu,
     &               qmin,lclose,rplsq,istat,rstat)

                call io_discard_write(1,t,nbod,ntp,xh,yh,zh,vxh,vyh,
     &               vzh,idt,tpfrac,xht,yht,zht,vxht,vyht,vzht,istat,
     &               rstat,iud,'discard.out',fopenstat,nleft)
             else
                nleft = ntp
             endif
	     
c if it is time, output orb. elements, 
	  if(t .ge. tout) then
c            apply corrector
	     do i=1,nbod
		xh2(i)=xh(i)
		yh2(i)=yh(i)
		zh2(i)=zh(i)
		vxh2(i)=vxh(i)
		vyh2(i)=vyh(i)
		vzh2(i)=vzh(i)
	     enddo
	     do i=1,ntp
		xht2(i)=xht(i)
		yht2(i)=yht(i)
		zht2(i)=zht(i)
		vxht2(i)=vxht(i)
		vyht2(i)=vyht(i)
		vzht2(i)=vzht(i)
	     enddo
	     call corrector(nbod,npl,ntp,t,mass,j2rp2,j4rp4,xh2,yh2,zh2,
     &            vxh2,vyh2,vzh2,xht2,yht2,zht2,vxht2,vyht2,vzht2,istat,
     &            dt,tinc)
             if(btest(iflgchk,1))  then    ! bit 1 is set
                call  io_write_frame_r(t,nbod,npl,ntp,mass,xh2,yh2,zh2,
     &               vxh2,vyh2,vzh2,idt,tpfrac,xht2,yht2,zht2,vxht2,
     &               vyht2,vzht2,istat,outfile,iub,frame,qrec,fopenstat)
	     endif
	    tout = tout + dtout
	  endif

c If it is time, do a dump
          if(t.ge.tdump) then
c            apply corrector
	     do i=1,nbod
		xh2(i)=xh(i)
		yh2(i)=yh(i)
		zh2(i)=zh(i)
		vxh2(i)=vxh(i)
		vyh2(i)=vyh(i)
		vzh2(i)=vzh(i)
	     enddo
	     do i=1,ntp
		xht2(i)=xht(i)
		yht2(i)=yht(i)
		zht2(i)=zht(i)
		vxht2(i)=vxht(i)
		vyht2(i)=vyht(i)
		vzht2(i)=vzht(i)
	     enddo
	     call corrector(nbod,npl,ntp,t,mass,j2rp2,j4rp4,xh2,yh2,zh2,
     &            vxh2,vyh2,vzh2,xht2,yht2,zht2,vxht2,vyht2,vzht2,istat,
     &            dt,tinc)

             tfrac = (t-t0)/(tstop-t0)
             write(*,998) t,tfrac,nleft
 998         format(' Time = ',1p1e12.5,': fraction done = ',0pf5.3,
     &            ': Number of active tp =',i4)

	     call io_dump_pl('backuppl.dat',nbod,npl,mass,xh2,yh2,zh2,
     &                 vxh2,vyh2,vzh2,lclose,iflgchk,rplsq,j2rp2,j4rp4)
	     call io_dump_tp('dump_tp.dat',ntp,idt,tpfrac,xht2,yht2,
     &                      zht2,vxht2,vyht2,vzht2,istat,rstat)
	     call io_dump_param('backupparam.dat',t,tstop,dt,tinc,
     &           dtout,dtdump,iflgchk,rmin,rmax,rmaxu,qmin,
     &           lclose,rcrit,qrec,fluxdist,outfile,frame)
	     call io_dump_pl('dump_pl.dat',nbod,npl,mass,xh2,yh2,zh2,
     &                 vxh2,vyh2,vzh2,lclose,iflgchk,rplsq,j2rp2,j4rp4)
	     call io_dump_tp('dump_tp.dat',ntp,idt,tpfrac,xht2,yht2,
     &                      zht2,vxht2,vyht2,vzht2,istat,rstat)
	     call io_dump_param('dump_param.dat',t,tstop,dt,tinc,
     &           dtout,dtdump,iflgchk,rmin,rmax,rmaxu,qmin,
     &           lclose,rcrit,qrec,fluxdist,outfile,frame)
	     tdump = tdump + dtdump

             if(btest(iflgchk,2))  then    ! bit 2 is set
                call anal_energy_write(t,nbod,mass,j2rp2,j4rp4,xh2,
     &               yh2,zh2,vxh2,vyh2,vzh2,iue,fopenstat,eoff,nclose)
             endif
             if(btest(iflgchk,3))  then    ! bit 3 is set
                call anal_jacobi_write(t,nbod,ntp,mass,xh2,yh2,zh2,vxh2,
     &               vyh2,vzh2,xht2,yht2,zht2,vxht2,vyht2,vzht2,istat,2,
     &               iuj,fopenstat)
             endif

	  endif
	enddo
c********** end of the big loop from time 't0' to time 'tstop'

c Do a final dump for possible resumption later 
c       apply corrector
	do i=1,nbod
	   xh2(i)=xh(i)
	   yh2(i)=yh(i)
	   zh2(i)=zh(i)
	   vxh2(i)=vxh(i)
	   vyh2(i)=vyh(i)
	   vzh2(i)=vzh(i)
	enddo
	do i=1,ntp
	   xht2(i)=xht(i)
	   yht2(i)=yht(i)
	   zht2(i)=zht(i)
	   vxht2(i)=vxht(i)
	   vyht2(i)=vyht(i)
	   vzht2(i)=vzht(i)
	enddo
	call corrector(nbod,npl,ntp,t,mass,j2rp2,j4rp4,xh2,yh2,zh2,
     &     vxh2,vyh2,vzh2,xht2,yht2,zht2,vxht2,vyht2,vzht2,istat,
     &     dt,tinc)
	call io_dump_pl('backuppl.dat',nbod,npl,mass,xh2,yh2,zh2,
     &            vxh2,vyh2,vzh2,lclose,iflgchk,rplsq,j2rp2,j4rp4)
	call io_dump_tp('backuptp.dat',ntp,idt,tpfrac,xht2,yht2,zht2,
     &              vxht2,vyht2,vzht2,istat,rstat)
	call io_dump_param('backupparam.dat',t,tstop,dt,tinc,
     &         dtout,dtdump,iflgchk,rmin,rmax,rmaxu,qmin,
     &         lclose,rcrit,qrec,fluxdist,outfile,frame)
	call io_dump_pl('dump_pl.dat',nbod,npl,mass,xh2,yh2,zh2,
     &            vxh2,vyh2,vzh2,lclose,iflgchk,rplsq,j2rp2,j4rp4)
	call io_dump_tp('dump_tp.dat',ntp,idt,tpfrac,xht2,yht2,zht2,
     &              vxht2,vyht2,vzht2,istat,rstat)
	call io_dump_param('dump_param.dat',t,tstop,dt,tinc,
     &         dtout,dtdump,iflgchk,rmin,rmax,rmaxu,qmin,
     &         lclose,rcrit,qrec,fluxdist,outfile,frame)

        call util_exit(0)
        end    ! swift_rmvs3.f
c---------------------------------------------------------------------
