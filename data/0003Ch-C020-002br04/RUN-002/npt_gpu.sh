#!/bin/bash
WD=`pwd`
#============================   DESKTOP  ===============================
TOP="../01-GENERATE_INITIAL/topol.top" 
GRO="../02-GROMACS_MIN-ALL/confout.gro"
MDP="./NPT_Berendsen.mdp"
TPR="../RUN-001/topol.tpr"
CPT="../RUN-001/state.cpt"

#============================   ACCELERA (Biophym1) ===============================
source /opt/gromacs-2018.3/bin/GMXRC.bash

EXE1="/opt/gromacs-2018.3/bin/gmx_mpi grompp"
EXE2="/opt/gromacs-2018.3/bin/gmx_mpi mdrun"
EXE3="/opt/gromacs-2018.3/bin/gmx_mpi convert-tpr"

module load gromacs-2018.3
${EXE3} -s ${TPR}  -o topol.tpr -until 20 >& out_grompp.dat
${EXE2} -ntomp 4 -pin on -s topol.tpr -cpi ${CPT} -noappend -v >& out_md.dat



