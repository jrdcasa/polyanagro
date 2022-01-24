#/bin/bash

module load PreMCParallel_V09
time PreMCparallel_v09.x psf ../namd.psf trjout.dcd trappeUA 0 1 0 
