#!/bin/bash
#SBATCH --partition=gpu
#SBATCH -n 16
#SBATCH -N 1
#SBATCH --gres=gpu:1
#SBATCH --mem-per-cpu=3G
#SBATCH --job-name=04-PVL04_40mon

source /opt/gromacs/gromacs-2024_x/bin/GMXRC.bash

MDP="NPT.mdp"
TOP="../02-RP_OPLS/PVLMon_residues_replicate.top"
GRO="../02-RP_OPLS/PVLMon_residues_replicate.gro"

gmx grompp -f $MDP -p $TOP -c $GRO -o new_topo.tpr --maxwarn 5 >& outgro.dat
gmx mdrun -s new_topo.tpr -v -noappend  >& outmd.dat
#gmx mdrun -s new_topo.tpr -v -noappend  >& outmd.dat
#echo 0|gmx trjconv -f traj_comp.part0001.xtc -o traj_unwrap.xtc -pbc whole -s new_topo.tpr

