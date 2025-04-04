#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --mem=2000M
#SBATCH --job-name=PE03
#SBATCH --exclude=aoki[01,04-08]
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
WD=`pwd`
cd $SLURM_SUBMIT_DIR

echo "Start Job: `date`"
../../src/Single_chain_mc_aoki00.x input.txt
echo "End   Job: `date`"
