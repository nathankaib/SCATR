c*************************************************************************
c                            IO_WRITE_CLO
c*************************************************************************
c writes out info on test particles and planets either when they are entering
c a close encounter or when they are leaving a close encounter
c
c             Input:
c                 i1stenc       ==>  integer indicating 1st enounter
c                 time          ==>  current time (real scalar)
c                 nbod          ==>  number of massive bodies (int scalar)
c                 ntp            ==>  number of massive bodies (int scalar)
c                 mass          ==>  mass of bodies (real array)
c                 xh,yh,zh      ==>  current position in helio coord 
c                                    (real arrays)
c                 vxh,vyh,vzh   ==>  current velocity in helio coord 
c                                    (real arrays)
c                 xht,yht,zht    ==>  current part position in helio coord 
c                                      (real arrays)
c                 vxht,vyht,vzht ==>  current velocity in helio coord 
c                                        (real arrays)
c                 istat           ==>  status of the test paricles
c                                      (2d integer array)
c                                      istat(i,1) = 0 ==> active:  = 1 not
c                                      istat(i,2) = -1 ==> Danby did not work
c                 oname           ==> output file name (character string) 
c                 iu              ==> unit number to write to
c                 fopenstat       ==>  The status flag for the open 
c                                      statements of the output files.  
c                                          (character*80)
c
c
c Remarks: Based on io_write_frame
c Authors:  Nathan Kaib 
c Date:    6/23/05
c Last revision: 

      subroutine io_write_clo(i1stenc,time,ienc,xht,yht,zht,vxht,vyht,
     &     vzht,xh,yh,zh,vxh,vyh,vzh,fopenstat,ntp,nbod,mass)

      include '../swift.inc'
      include 'io.inc'

c...  Inputs:
      integer i1stenc,i,ntp,ienc(NTPMAX),nbod
      real*8 xht(NTPMAX),yht(NTPMAX),zht(NTPMAX),mass(NPLMAX)
      real*8 vxht(NTPMAX),vyht(NTPMAX),vzht(NTPMAX)
      real*8 xh(NPLMAX),yh(NPLMAX),zh(NPLMAX)
      real*8 vxh(NPLMAX),vyh(NPLMAX),vzh(NPLMAX)
      real*8 time
      character*80 fopenstat

c...  Internals
      integer ierr,iu,ialpha,pl
      real*8 a,e,inc,capom,omega,capm
      real*8 gm

c----
c...  Executable code 

c...  if first time through open file
      iu=70
      gm=mass(1)
      if(i1stenc.eq.0) then
         call io_open(iu,'clo.dat',fopenstat,'UNFORMATTED',ierr)
         if(ierr.ne.0) then
           write(*,*) ' SWIFT ERROR: in io_write_clo: '
           write(*,*) '     Could not open binary output file:'
           call util_exit(1)
         endif
         i1stenc = 1
      else
        call io_open(iu,'clo.dat','append','UNFORMATTED',ierr)
      endif
      
      do i=1,ntp
         if (ienc(i).gt.0) then
            call orbel_xv2el(xht(i),yht(i),zht(i),vxht(i),vyht(i),
     &           vzht(i),gm,ialpha,a,e,inc,capom,omega,capm)
c     writing test particle data
            call io_write_line_r(iu,i,time,a,e,inc,capom,omega,capm,
     &           xht(i),yht(i),zht(i),vxht(i),vyht(i),vzht(i))
c     writing planetary data
            pl=ienc(i)
            call io_write_line_r(iu,pl,time,a,e,inc,capom,omega,capm,
     &           xh(pl),yh(pl),zh(pl),vxh(pl),vyh(pl),vzh(pl))
         endif
      enddo

      close(unit = iu)
      return
      end
