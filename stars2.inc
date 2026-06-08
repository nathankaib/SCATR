c*************************************************************************
c                        STARS2.INC
c*************************************************************************
C Include file for the stellar encounter stuff.
c
c Author:  Hal Levison  
c Date:    10/16/98
c Last revision: 

c...   Maximum number of stellar encounters
       integer  NSTENC
       parameter  (NSTENC = 2000000) 

c...   Number of lines in the star model
       integer  NMODEL
       parameter  (NMODEL = 20) 

c...   Radius that encounters will happen
       real*8 R0_STAR
       parameter (R0_STAR=2.10d5)   ! one pc

c...   Ratio of timestep to `orbital period' to resolve peri
       real*8 RAT_T
       parameter (RAT_T=600.0d0)  

c...   Radius of a star in AU
       real*8 RSTAR
       parameter (RSTAR=0.005d0)
