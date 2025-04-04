#/bin/bash
module load gromacs_2020_2


time gmx_mpi trjcat -f ../01-RESTART-0000-1000ns/traj_comp.xtc ../02-RESTART-1000-2000ns/traj_comp.part0002.xtc  ../03-RESTART-2000-3000ns/traj_comp.part0003.xtc -o trjout.xtc

rm -f \#*
