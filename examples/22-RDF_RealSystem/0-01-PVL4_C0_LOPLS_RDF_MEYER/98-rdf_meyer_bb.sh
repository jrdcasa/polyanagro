#!/bin/bash

LAMMPS2YASP=/home/jramos/Programacion/hmeyer/yasp_cJ/bin/lammpstxt2yasp
RDF=/home/jramos/Programacion/hmeyer/yasp_cJ/bin/rdf_all.x

LAMMPSTRJ=01-PVL4_C0_LOPLS_trj.lammpstrj
YASPTRJ=01-PVL4_C0_LOPLS_trj.trj

NATCH=603
NBINS=100
RMAX=20

# Convert trajectory
${LAMMPS2YASP} -o ${YASPTRJ} -i ${LAMMPSTRJ}

# Intramolecular 
$RDF ${NBINS} ${RMAX}   ${YASPTRJ} <template_bb.dat >out_intra_bb.dat
$RDF ${NBINS} ${RMAX} -inter ${NATCH} ${YASPTRJ} <template_bb.dat >out_inter_bb.dat

