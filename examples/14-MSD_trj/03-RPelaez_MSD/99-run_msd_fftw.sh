#!/bin/bash

TRJ=01-PVL4_C0_LOPLS_trj.spunto
MSD=/home/jramos/Programacion/GITHUB_OTHERS/MSD_CPU_GPU_RPelaez/MeanSquareDisplacement/cmake-build-debug/bin/msd
NATOMS=12060
NFRAMES=501
DIM=3

cat ${TRJ} | ${MSD} -N ${NATOMS} -Nsteps ${NFRAMES} -dimensions ${DIM}  >out_msd.dat

