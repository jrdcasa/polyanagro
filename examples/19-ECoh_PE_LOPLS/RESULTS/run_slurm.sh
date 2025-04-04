#!/bin/bash
#SBATCH -p master
#SBATCH -n 8
#SBATCH --job-name="ECOH_01PE_450K"

source /optnfs/CODES_SRC/sandbox_common/bin/activate

# Aoki12 --> rxt_4090 (gromacs 2023)
#GMXFULLPATH="/optnfs/gromacs/gromacs-2023/gpu_rtx4090/bin/gmx"
# Aoki05 --> gtx780 (gromacs 2022) TPR seems to be compatible with the gromacs2023
#GMXFULLPATH="/optnfs/gromacs/gromacs-2022/gpu_gtx780/bin/gmx_gpu"
# AokiMaster --> master (gromacs 2023) 
GMXFULLPATH="/optnfs/gromacs/gromacs-2023/gromacs_cpu_aokimaster/bin/gmx"
source /optnfs/gromacs/gromacs-2023/gromacs_cpu_aokimaster/bin/GMXRC.bash

BASE="../../../01-PE_450K"
XTC="${BASE}/03-PROD_0000_0500ns_450K/traj_comp.part0001.xtc"
TPR="${BASE}/03-PROD_0000_0500ns_450K/new_topo.tpr"
TOPO="${BASE}/00-RP_TOPO_LOPLS/trajectory_topo_replicate.top"
FRACTION=0.25

gmx trjconv -f ${XTC} -b 250000 -skip 1000 -o trj_out.xtc

energy_cohesive -t trj_out.xtc --gmxfullpath ${GMXFULLPATH} --tpr ${TPR} --topo ${TOPO} --fraction_trj_avg ${FRACTION}
