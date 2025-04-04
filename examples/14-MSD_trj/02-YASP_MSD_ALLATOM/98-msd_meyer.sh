#!/bin/bash

LAMMPS2YASP=/home/jramos/Programacion/hmeyer/yasp_cJ/bin/lammpstxt2yasp
MSD=/home/jramos/Programacion/hmeyer/yasp_cJ/bin/msdmol.x

LAMMPSTRJ=01-PVL4_C0_LOPLS_trj.lammpstrj
YASPTRJ=01-PVL4_C0_LOPLS_trj.trj

NATCH=603
NMOL=20
DT=20  #ps

# Convert trajectory
${LAMMPS2YASP} -o ${YASPTRJ} -i ${LAMMPSTRJ} -t ${DT}

# MSD
$MSD ${NMOL} ${YASPTRJ} <template_all.dat >out_msd.dat
#$MSD ${NBINS} ${RMAX} -inter ${NATCH} ${YASPTRJ} <template_all.dat >out_inter_all.dat

