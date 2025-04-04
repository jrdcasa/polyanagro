#!/bin/bash
WK=`pwd`

LISTDIR=(/home/jramos/PycharmProjects/polyanagro/data/0100Ch-C500-000br02/T-450K/06-NPT-100ns-Parrinello/)
module load gromacs-2018.3
#
# This script requests the following tools installed and accesible by $PATH:
#         * gmxdump
#         * catdcd
#         * PreMCparallel_V08.x
#
#
#  The tree directory has to be as follows:
#    1) All calculations has to be in a directory with the following pattern:
#             [Index %2s]-[Temp]K (Example: 01-600K)
#    2) Inside this directory the subdirectories with the different restarts
#       has to follow the pattern:
#            R-[Index %2s]_[Initime %4s]-[Endtime %4s]ns 
#                   (Example: R-01_0000-0100ns)
#   3) The script produces an ANALYSIS_PREMC dir containing all results from
#      the analysis.
#
###############################################################################

LEN=${#LISTDIR[@]}
index=0
TOTAL=$LEN

for (( i=0; i<$LEN; i++ )); do

  idir=${LISTDIR[$i]}
  echo "=========================================="
  echo "   Processing...$idir. $index of $TOTAL   "
  echo "   `date`                                 "
  echo "=========================================="
  cd $idir
  XTCFILES=`ls -d *REST*/*.xtc`
  if [ ! -e traj.dcd ]; then
      gmx_mpi trjcat -f $XTCFILES -o trajout.xtc
      catdcd -o traj.dcd -xtc trajout.xtc
      rm trajout.xtc
  fi

  DT=`gmx_mpi dump -s 01-*/*.tpr 2>out.dat | egrep dt | awk {'print $3'}`
  XTCOUT=`gmx_mpi dump -s 01-*/*.tpr 2>out.dat | egrep nstxout-compressed | awk {'print $3'}`
  TIMEDUMP=`echo "$DT*$XTCOUT" |bc`
  echo "$TIMEDUMP ps" >timestepframes.dat

  cp ./00-GENERATE/namd.psf ./00-GENERATE/trappeUA.dat .
  /home/jramos/PycharmProjects/polyanagro/data/0100Ch-C500-000br02/T-450K/06-NPT-100ns-Parrinello/PreMCparallel_v09.x psf namd.psf traj.dcd trappeUA 0 1 0
  rm ./namd.psf ./trappeUA.dat
  mkdir -p ANALYSIS_PREMC
  mv *.dat ANALYSIS_PREMC

  index=`echo "$index+1"|bc`
  cd $WK
  echo "******************************************"
  echo "   Finnishing...$idir.                    "
  echo "   `date`                                 "
  echo "******************************************"

done
