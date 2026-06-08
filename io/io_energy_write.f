c*************************************************************************
c                          IO_ENERGY_WRITE.F
c*************************************************************************
c Does the write for anal_jacobi_write
c
c      Input:
c            i1st           ==>  =0 if first write, =1 if not (int scalar)
c            t              ==>  current time (real scalar)
c            energy         ==>  Total energy
c            eltot          ==>  components of total angular momentum
c                               (real array)
c            iu             ==>  unit to write to
c            fopenstat      ==>  The status flag for the open 
c                                statements of the output files.  
c                                          (character*80)
c
c Remarks: 
c Authors:  Hal Levison 
c Date:    2/21/94
c Last revision: 3/4/94

c changing this so flux of comets with 1 AU of sun is tacked onto energy.out
c (nathan kaib 8/20/04)

      subroutine io_energy_write(i1st,t,energy,eltot,iu,fopenstat,flux)

      include '../swift.inc'
      include 'io.inc'

c...  Inputs: 
      integer iu,i1st
      real*8 energy,eltot(3)
      real*8 t
      character*(*) fopenstat

c...  Internals
      integer ierr

      integer flux

c----
c...  Executable code 

      if(i1st.eq.0) then

         call io_open(iu,'energy.out',fopenstat,'FORMATTED',ierr)
         if(ierr.ne.0) then
            write(*,*) ' SWIFT ERROR: in anal_energy_write '
            write(*,*) '     Could not open energy.out '
            call util_exit(1)
         endif
         
      else
         
         call io_open(iu,'energy.out','append','FORMATTED',ierr)

      endif

c     tacking on flux of comets
      write(iu,2) t,energy,eltot,flux
 2    format(1x,1p1e12.5,4(2x,1p1e23.16),2x,i6)

      close(iu)

      return
      end                       ! io_energy_write
c-------------------------------------------------------------------------
