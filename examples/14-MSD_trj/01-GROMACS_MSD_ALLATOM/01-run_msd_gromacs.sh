#!/bin/bash
source /opt/gromacs/gromacs-2024_x/bin/GMXRC.bash

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <xtc trajectory> <tpr file> <unwrap: 0:No 1:Yes>"
    exit 1
fi

if [[ "$3" != "0" && "$3" != "1" ]]; then
    echo "Error: Second Argument must be 0 or 1."
    exit 1
fi
TRJ=$1
TPR=$2

if [[ $3 == "1" ]]; then
    echo 0 | gmx trjconv -f ${TRJ} -s ${TPR} -pbc whole -o trajout_whole.xtc
    echo 0 | gmx trjconv -f trajout_whole.xtc -pbc nojump -o trajout_nojump.xtc
    TRJ=trajout_nojump.xtc
fi

#echo 0 | gmx msd   -f ${TRJ} -dt ${SAMPLE_DT_MSD} -s ${TPR} -o msd_allinternal_0000-0010ns
echo 0 | gmx msd   -f ${TRJ} -s ${TPR} -o msd_allinternal_0000-0010ns






