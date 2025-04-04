#/bin/bash

module load PreMCParallel_V09
time PreMCparallel_v09.x psf ../00-GENERATE/namd.psf trjout.dcd trappeUA 0 1 0 
