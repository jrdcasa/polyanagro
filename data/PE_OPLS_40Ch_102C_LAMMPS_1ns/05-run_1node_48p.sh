#!/bin/bash
#SBATCH --partition=generic
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --cpus-per-task=1
#SBATCH --mem=4000M
#SBATCH --job-name=L01

module purge
#module use /dragofs/sw/easybuild/modules/all/
module load foss/2020a
module load LAMMPS

WD=`pwd`

#mpirun -np 48 lmp -in PE102_40Ch_residues_replicate_clean.inp -log PE102_40Ch_residues_replicate_clean.log >output.dat
lmp -restart2data restart.data.500000 restart.data.500000.data

