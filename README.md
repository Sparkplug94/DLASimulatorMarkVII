# DLASimulatorMarkVII
Simple Particle Tracking Code for DLA

CURRENT ISSUES:
Outlier particles and complex numbers - currently, some particles appear to get imaginary kicks. i.e. a momentum kick from a DLA that is not a real number. I have no idea why, but i suspect this is the source of the "outlier" particles, with energies waaaay outside the norm. 

Current manifestation of this problem: DLAUpdate function is fine when stepped through line by line. produces no complex vales in phaseSpace object. for eps = 100e6, rn = 0, theta_r = 0, it works up to 13 times iteratively, and on the 14th, it very reliably produces a complex value. I have no idea why. 
