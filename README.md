# DLASimulatorMarkVII
Simple Particle Tracking Code for DLA

CURRENT ISSUES:
Particle survival suspiciously low.
plotting function input arguments take phaseSpace.dist instead of phaseSpace
particle survival goes UP when amplitude and phase noise are added. Why?


Current best design: v6_latticeDesign_withShift.gds
from an email sent to Ben Cowan explaining the design
"Here's the first cut of APF lattices. There's two, one where I deliberately ignore the slow phase drift of the acceleration relative to the synchronous particle, and one where I don't. I'm making a lot of assumptions with this lattice - I'm doing all my simulations of structure efficiency in 2D, for one thing. Actually, if you have time, it might be helpful to run my design flow past you via Zoom, to see if you spot any obvious mistakes I've made.

Anyway, I'm currently simulating a plane wave field of 540 MV/m from both sides, which should lead to ~200 MeV/m peak gradient, and ~100 MeV/m average gradient in the structure. My injection energy is 90 keV and wavelength is 2 um. I'm mainly interested in whether any particles at all make it through the lattice. For comparison, my simulation predicts quite low, but nonzero, transmission through the entire structure."


TO USE:
first run "makeTaperedLatticeMatFile.m" in Matlab to generate the base for the lattice (how many APF periods, the tapering of the APF period, etc)
next, run "run_v6.m", which will export a lattice design .csv file
next, run "loadMatlabLatticeOutput_WriteGDS.opynb" in jupyter notebook, which will write the GDS

note that the structure constant fit function in run_v6 and the pillar shift relative to period centroid 
are hardcoded in for a particular pillar design that I simulated in Lumerical. if you want to change
 that you have to redo the lumerical sims

look in C:\Users\dylan\Documents\Lumerical Simulations\2D_Optimized_Pillars for the last time I did this
the analysis scripts are in C:\Users\dylan\Documents\MATLAB\Analyze Lumerical Output