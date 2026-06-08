# SCATR

This is the README file for the SCATR integrator, which is described in Kaib, Quinn, & Brasser (2011): https://iopscience.iop.org/article/10.1088/0004-6256/141/1/3/pdf

The code is nearly identical to the SWIFT RMVS3 integrator, and users should consult the documentation for this package to attain an understanding of most of SCATR's functionality: https://www2.boulder.swri.edu/~hal/swift.html

There are, however, a few key differences:

1) The major one is that particles are drifted about a planetary system's barycenter when they exceed some critical distance. In addition, in this barycentric routine, a larger timestep can be used. The large barycentric timestep, in days, as the 3rd number of the first line of the param.in file (see 'example' simulation directory). The number of smaller heliocentric timesteps taken during a single barycentric one is specified as the 4th number in this same line in param.in.

2) Only particles with perihelia below a critical pericenter are recorded when orbital elements are written to bin.dat. This critical pericenter, in au, is given in the 5th line of param.in (see 'example' simulation directory).

3) The orbits of particles making pericenter passages inside a given pericenter distance are recorded in a binary file named flux.dat. This critical pericenter passage distance is given in line 6th line of param.in (see 'example' simulation directory).

4) The distance, in au, at which the integration drift switches from heliocentric to barycentric is given in the 7th line of param.in (see 'example' simulation directory).

5) The frame (heliocentric or barycentric) in which particles' inputs, dumps and binary outputs are written and read is specified in the 8th line of param.in as 'hel' or 'bar' (see 'example' simulation directory).

6) The total number of massive bodies bound by gravity is given as the 2nd number in the first line of pl.in (see 'example' simulation directory). This is done in case the user wishes to introduce short-lived unbound particles to the simulation, such as passing stars.

7) Each of the test particles is assigned an ID number in the tp.in file. This is the first number in the line above a particle's Cartesian positions in the tp.in file (see 'example' simulation directory).

8) Each of the test particles is also given a 'weight' in the tp.in file. This is the second number in the line above a particle's Cartesian positions in the tp.in file (see 'example' simulation directory). This is done in case a use wishes to implement particle cloning within the running of a simulation.

That's all the differences I can think of!
