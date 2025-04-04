#!/bin/bash
WD=`pwd`
#============================   DESKTOP  ===============================
TOP="../01-GENERATE_INITIAL/topol.top" 
GRO="../02-GROMACS_MIN-ALL/confout.gro"
MDP="./NPT_Berendsen.mdp"


#============================   ACCELERA (Biophym1) ===============================
source /opt/gromacs-2018.3/bin/GMXRC.bash

EXE1="/opt/gromacs-2018.3/bin/gmx_mpi grompp"
EXE2="/opt/gromacs-2018.3/bin/gmx_mpi mdrun"

module load gromacs-2018.3
${EXE1} -f $MDP -p $TOP -c $GRO -maxwarn 10 >& out_grompp.dat
${EXE2} -ntomp 4 -pin on -v >& out_md.dat


