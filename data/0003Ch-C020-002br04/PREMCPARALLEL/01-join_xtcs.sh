#/bin/bash
module load gromacs_2020_2


time gmx_mpi trjcat -f ../RUN-001/traj_comp.xtc ../RUN-002/traj_comp.part0002.xtc  ../RUN-003/traj_comp.part0003.xtc -o trjout.xtc

rm -f \#*
