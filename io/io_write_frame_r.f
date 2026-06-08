c*************************************************************************
c                            IO_WRITE_FRAME_R
c*************************************************************************
c write out a whole frame to an real*4 binary file.
c both massive and test particles
c
c             Input:
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
c Authors:  Hal Levison 
c Date:    2/22/94
c Last revision: 

      subroutine io_write_frame_r(time,nbod,npl,ntp,mass,xh,yh,zh,vxh,
     &           vyh,vzh,idt,tpfrac,xht,yht,zht,vxht,vyht,vzht,istat,
     &           oname,iu,frame,qrec,fopenstat)

      include '../swift.inc'
      include 'io.inc'

c...  Inputs: 
      integer nbod,ntp,iu,nleft,npl
      real*8 mass(NPLMAX),dummy,mtot
      real*8 time
      integer istat(NTPMAX,NSTAT),idt(NTPMAX)
      real*8 xh(NPLMAX),yh(NPLMAX),zh(NPLMAX)
      real*8 vxh(NPLMAX),vyh(NPLMAX),vzh(NPLMAX)
      real*8 xht(NTPMAX),yht(NTPMAX),zht(NTPMAX)
      real*8 vxht(NTPMAX),vyht(NTPMAX),vzht(NTPMAX),tpfrac(NTPMAX)
      character*80 oname,fopenstat,frame
      real*8 qrec,q
      real*8 xc,yc,zc,vxc,vyc,vzc,rad

c...  Internals
      integer i,id
      integer ialpha,ierr
      real*8 a,e,inc,capom,omega,capm
      real*8 gm
      integer i1st    ! =0 first time through; =1 after
      data i1st/0/
      save i1st

c----
c...  Executable code

c...  if first time through open file
      if(i1st.eq.0) then
         call io_open(iu,oname,fopenstat,'UNFORMATTED',ierr)
         if(ierr.ne.0) then
           write(*,*) ' SWIFT ERROR: in io_write_frame: '
           write(*,*) '     Could not open binary output file:'
           call util_exit(1)
         endif
         i1st = 1
      else
        call io_open(iu,oname,'append','UNFORMATTED',ierr)
      endif

c...  only outputting particles that have q < 40 AU
      mtot=0.0d0
      do i=1,npl
         mtot = mtot+mass(i)
      enddo
      gm=mtot
      xc=0.0d0
      yc=0.0d0
      zc=0.0d0
      vxc=0.0d0
      vyc=0.0d0
      vzc=0.0d0
      do i=2,npl
         xc=xc+xh(i)*mass(i)
         yc=yc+yh(i)*mass(i)
         zc=zc+zh(i)*mass(i)
         vxc=vxc+vxh(i)*mass(i)
         vyc=vyc+vyh(i)*mass(i)
         vzc=vzc+vzh(i)*mass(i)
      enddo
      xc=xc/gm
      yc=yc/gm
      zc=zc/gm
      vxc=vxc/gm
      vyc=vyc/gm
      vzc=vzc/gm

      nleft = 0
      gm = mtot
      do i=1,ntp
         if(istat(i,1).eq.0) then
            if (istat(i,NSTAT-2).eq.1) then 
               call orbel_xv2el(xht(i)-xc,yht(i)-yc,zht(i)-zc,
     &              vxht(i)-vxc,vyht(i)-vyc,vzht(i)-vzc,gm,ialpha,a,
     &              e,inc,capom,omega,capm)
            else
               call orbel_xv2el(xht(i),yht(i),zht(i),vxht(i),vyht(i),
     &              vzht(i),gm,ialpha,a,e,inc,capom,omega,capm)
            endif
            q=a*(1.0d0 - e)
            if (q.lt.qrec) then
               nleft=nleft+1
            endif
         endif
      enddo

      call io_write_hdr_r(iu,time,nbod,nleft)

c...  write out planets
      do i=2,nbod
         gm = mass(1)
         id = -1*i
 	 call orbel_xv2el(xh(i),yh(i),zh(i),vxh(i),vyh(i),vzh(i),gm,
     &          ialpha,a,e,inc,capom,omega,capm)

         call io_write_line_r(iu,id,0,mass(i),a,e,inc,capom,omega,capm)
      enddo

c...  write out test particles
      if (frame(1:3).eq.'bar'.or.frame(1:3).eq.'BAR') then
         gm = mtot
         do i=1,ntp
            if(istat(i,1).eq.0) then
               if (istat(i,NSTAT-2).eq.1) then 
                  call orbel_xv2el(xht(i)-xc,yht(i)-yc,zht(i)-zc,
     &                 vxht(i)-vxc,vyht(i)-vyc,vzht(i)-vzc,gm,ialpha,a,
     &                 e,inc,capom,omega,capm)
               else
                  call orbel_xv2el(xht(i),yht(i),zht(i),vxht(i),vyht(i),
     &                 vzht(i),gm,ialpha,a,e,inc,capom,omega,capm)
               endif
               q=a*(1.0d0 - e)
               if (q.lt.qrec) then
                  call io_write_line_r(iu,idt(i),istat(i,NSTAT-3),
     &                 tpfrac(i),a,e,inc,capom,omega,capm)

c                  call io_write_line_r(iu,idt(i),tpfrac(i),xht(i),
c     &                 yht(i),zht(i),vxht(i),vyht(i),vzht(i))
               endif
            endif
         enddo
      else
         gm = mass(1)
         do i=1,ntp
            if(istat(i,1).eq.0) then
               if (istat(i,NSTAT-2).eq.0) then 
                  call orbel_xv2el(xht(i)+xc,yht(i)+yc,zht(i)+zc,
     &                 vxht(i)+vxc,vyht(i)+vyc,vzht(i)+vzc,gm,ialpha,a,
     &                 e,inc,capom,omega,capm)
               else
                  call orbel_xv2el(xht(i),yht(i),zht(i),vxht(i),vyht(i),
     &                 vzht(i),gm,ialpha,a,e,inc,capom,omega,capm)
               endif
               q=a*(1.0d0 - e)
               if (q.lt.qrec) then
c                  call io_write_line_r(iu,idt(i),
c     &                 tpfrac(i),xht(i),yht(i),zht(i),vxht(i),vyht(i),
c     &                 vzht(i))
                  call io_write_line_r(iu,idt(i),istat(i,NSTAT-3),
     &                 tpfrac(i),a,e,inc,capom,omega,capm)
               endif
            endif
         enddo
      endif

      close(iu)
      return
      end      ! io_write_frame_r
c----------------------------------------------------------------------
