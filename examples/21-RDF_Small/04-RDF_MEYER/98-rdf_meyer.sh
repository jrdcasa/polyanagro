#!/bin/bash

LAMMPS2YASP=/home/jramos/Programacion/hmeyer/yasp_cJ/bin/lammpstxt2yasp
RDF=/home/jramos/Programacion/hmeyer/yasp_cJ/bin/rdf_all.x

LAMMPSTRJ=01-PVMon.lammpstrj
YASPTRJ=01-PVMon.trj

NATCH=17
NBINS=200
RMAX=20

# Convert trajectory
${LAMMPS2YASP} -o ${YASPTRJ} -i ${LAMMPSTRJ}

# Intramolecular 
echo "$RDF ${NBINS} ${RMAX}   ${YASPTRJ} <template_all.dat >out_intra_all.dat"
$RDF ${NBINS} ${RMAX}   ${YASPTRJ} <template_all.dat >out_intra_all.dat
echo "$RDF ${NBINS} ${RMAX} -inter ${NATCH} ${YASPTRJ} <template_all.dat >out_inter_all.dat"
$RDF ${NBINS} ${RMAX} -inter ${NATCH} ${YASPTRJ} <template_all.dat >out_inter_all.dat

