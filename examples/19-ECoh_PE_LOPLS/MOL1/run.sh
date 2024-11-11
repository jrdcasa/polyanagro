#!/bin/bash

source /opt/gromacs/gromacs-2024_x/bin/GMXRC.bash

# Unwrap trj
echo 0 | gmx trjconv -f SC_traj_ch_0000.xtc -s SC_traj_ch_0000.tpr -pbc whole -o SC_traj_ch_0000_whole.xtc
echo 0 | gmx trjconv -f SC_traj_ch_0000.gro -s SC_traj_ch_0000.tpr -pbc whole -o SC_traj_ch_0000_whole.gro

# Rerun wrap
gmx grompp -f rerun.mdp -c SC_traj_ch_0000_whole.gro -p polymer_SC.top -o SC_traj_ch_0000_whole.tpr
gmx mdrun -rerun SC_traj_ch_0000_whole.xtc -s SC_traj_ch_0000_whole.tpr --deffnm whole
echo 1 2 3 4 5 6 7 8 9 0 | gmx energy -f whole.edr -o energy_whole.xvg

# Rerun unwrap
gmx grompp -f rerun.mdp -c SC_traj_ch_0000.gro -p polymer_SC.top -o SC_traj_ch_0000.tpr
gmx mdrun -rerun SC_traj_ch_0000.xtc -s SC_traj_ch_0000.tpr --deffnm wrap
echo 1 2 3 4 5 6 7 8 9 0 | gmx energy -f wrap.edr -o energy_wrap.xvg

rm \#*
