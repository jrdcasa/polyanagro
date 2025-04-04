#!/bin/bash
#SBATCH --partition=gpu_gtx780_16p
#SBATCH -n 8
#SBATCH -N 1
##SBATCH -w aoki04
#SBATCH --mem-per-cpu=3G
#SBATCH --job-name=msdint_PVL4_C1_500K

DIRXTC=../../../01-PVL4_C0_LOPLS/03-PREPARE_LIBRARY_POLYPLY_40Mon/04-MD_NPT_0000_0200ns_500K
PSFFILE=../../../01-PVL4_C0_LOPLS/04-PDB_LABELED/confout_unwrap_labeled.psf
TRJXTC=traj_comp.part0004.xtc 
TOPOLTPR=new_topo.tpr
SKIP=1000
# Use None for default start and end time
#USERBTIME=1500000 #ps
#USERETIME=1750000 #ps
USERBTIME=None #ps
USERETIME=None #ps
PDBOUT=01-PVL4_C0_LOPLS_trj.pdb
GMXSOURCE=/optnfs/gromacs/gromacs-2023/gromacs_cpu_aokimaster/bin/GMXRC.bash
#GMXSOURCE=/opt/gromacs/gromacs-2024_x/bin/GMXRC.bash

# ====================== MAIN ==========================
WK=`pwd`

echo "Job Start" >Info.log

if [[ -f ${GMXSOURCE} ]]; then
    source ${GMXSOURCE}
else
   echo "  Warning. GROMACS has not been sourced." 
   echo "  Warning. GROMACS has not been sourced." >>Info.log
fi

# Get first and last time in the trajectory
echo "   START Check trajectories `date`" >>Info.log
gmx check -f ${DIRXTC}/${TRJXTC} 2>tmp.dat
sed -i "s/\r/\n/g" tmp.dat  # Change ^M by return carriage
BTIME=`egrep "Reading frame" tmp.dat |head -1 |awk '{print $5}'`   #ps
ETIME=`egrep "Last frame" tmp.dat |head -1 |awk '{print $5}'`      #ps
NFRAMES=`egrep "Step" tmp.dat | awk '{print $2}'`
DT=`egrep "Step" tmp.dat | awk '{print $3}'`                               #ps
echo "   End Check trajectories `date`" >>Info.log

echo " =======================  TRAJECTORY INFO    =======================" >>Info.log
echo "    Trajectory: Start = ${BTIME} ps. Finnished at ${ETIME} ps. Delta_t = ${DT} ps. NFrames: ${NFRAMES}" >>Info.log
echo "        Start     = ${BTIME} ps" >>Info.log
echo "        Finnished = ${ETIME} ps" >>Info.log
echo "        Delta_t   = ${DT} ps" >>Info.log
echo "        NFrames   = ${NFRAMES}" >>Info.log
echo " ======================= END TRAJECTORY INFO =======================" >>Info.log
rm tmp.dat

# This converts the value of USERBTIME to lowercase.
if [[ ! "${USERBTIME,,}" == "none" ]]; then
    BTIME=${USERBTIME}    
fi
# This converts the value of USERETIME to lowercase.
if [[ ! "${USERETIME,,}" == "none" ]]; then
    ETIME=${USERETIME}    
fi

echo "" >>Info.log
echo " =======================    EXTRACTION INFO      =======================" >>Info.log
#NEWFRAMES=`echo ${NFRAMES} ${SKIP} | awk '{printf("%d",$1/$2)}'`
NEWDT=`echo ${DT} ${SKIP} | awk '{printf("%d",$1*$2)}'`
NEWFRAMES=`echo ${BTIME} ${ETIME} ${NEWDT} | awk '{printf("%d",($2-$1)/$3+1)}'`
echo "        Start           = ${BTIME} ps" >>Info.log
echo "        Finnished       = ${ETIME} ps" >>Info.log
echo "        Frames to write = ${NEWFRAMES}" >>Info.log
echo "        DT              = ${NEWDT} ps" >>Info.log
echo " =======================  END EXTRACTION INFO    =======================" >>Info.log

echo "" >>Info.log
echo "   START Extract to PDB trajectory `date`" >>Info.log
echo 0 0 | gmx trjconv -f ${DIRXTC}/${TRJXTC} -skip ${SKIP} -b ${BTIME} -e ${ETIME} -s ${DIRXTC}/${TOPOLTPR} -o ${PDBOUT}
echo "     Trajectory has been written to ${PDBOUT}" >>Info.log
echo "   END   Extract to PDB trajectory `date`" >>Info.log

cp ${PSFFILE} ./
echo "" >>Info.log
echo "   ${PSFFILE} copied" >>Info.log

rm \#*
echo "Job Done!!!!!" >>Info.log
