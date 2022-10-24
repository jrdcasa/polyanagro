#!/bin/bash
# ========================== MAIN =====================================
module purge
WK=`pwd`
rm -rf XX-CHECK_POLYANAGRO
cd $WK
DATAFOLDER=`realpath ../data`

if [[ $# -eq 0 ]]; then
    runjob=0        # Run all examples
else
    runjob=$1
fi

echo "========== TEST EXAMPLES ==========" >${WK}/out_tests.log

source /home/jramos/Programacion/sandboxes/sandbox_common/bin/activate
mkdir -p XX-CHECK_POLYANAGRO

# Example 01 ************************************
if [[ $runjob -eq 0 || $runjob -eq 1 ]]; then
    echo "------------------------------" >>${WK}/out_tests.log
    echo "Example 01: Polyanagro: info_trj GROMACS" >>${WK}/out_tests.log
    echo "------------------------------" >>${WK}/out_tests.log
    TESTDIR=01-Info_trj_gromacs
    TESTDIR_TARGET=${WK}/XX-CHECK_POLYANAGRO/${TESTDIR}
    mkdir -p ${TESTDIR_TARGET}
    cd ${TESTDIR_TARGET}
    info_trj -t ${DATAFOLDER}/0003Ch-C020-002br04/RUN-001/traj_comp.xtc ${DATAFOLDER}/0003Ch-C020-002br04/RUN-002/traj_comp.part0002.xtc --topo ${DATAFOLDER}/0003Ch-C020-002br04/RUN-001/topol.tpr
    cd $WK
fi
# Example 01 ************************************

# Example 02 ************************************
if [[ $runjob -eq 0 || $runjob -eq 2 ]]; then
    echo "------------------------------" >>${WK}/out_tests.log
    echo "Example 02: Polyanagro: info_trj LAMMPS" >>${WK}/out_tests.log
    echo "------------------------------" >>${WK}/out_tests.log
    TESTDIR=02-Info_trj_lammps
    TESTDIR_TARGET=${WK}/XX-CHECK_POLYANAGRO/${TESTDIR}
    mkdir -p ${TESTDIR_TARGET}
    cd ${TESTDIR_TARGET}
    info_trj -t ${DATAFOLDER}/PE_OPLS_40Ch_102C_LAMMPS_1ns/dump.dcd --topo ${DATAFOLDER}/PE_OPLS_40Ch_102C_LAMMPS_1ns/restart.data.500000.data
    cd $WK
fi
# Example 02 ************************************

# Example 03 ************************************
if [[ $runjob -eq 0 || $runjob -eq 3 ]]; then
    echo "------------------------------" >>${WK}/out_tests.log
    echo "Example 03: Polyanagro: info_trj PDB" >>${WK}/out_tests.log
    echo "------------------------------" >>${WK}/out_tests.log
    TESTDIR=03-Info_trj_pdb
    TESTDIR_TARGET=${WK}/XX-CHECK_POLYANAGRO/${TESTDIR}
    mkdir -p ${TESTDIR_TARGET}
    cd ${TESTDIR_TARGET}
    info_trj -t ${DATAFOLDER}/PE_OPLS_40Ch_102C_LAMMPS_1ns/dump.dcd --topo ${DATAFOLDER}/PE_OPLS_40Ch_102C_LAMMPS_1ns/PE102_40Ch_residues_replicate.pdb
    cd $WK
fi
# Example 03 ************************************

# Example 04 ************************************
if [[ $runjob -eq 0 || $runjob -eq 4 ]]; then
    echo "------------------------------" >>${WK}/out_tests.log
    echo "Example 04: Polyanagro: replicate_polymer" >>${WK}/out_tests.log
    echo "------------------------------" >>${WK}/out_tests.log
    TESTDIR=04-Melt_100ChPE500_TrappeUA
    TESTDIR_TARGET=${WK}/XX-CHECK_POLYANAGRO/${TESTDIR}
    mkdir -p ${TESTDIR_TARGET}
    cp ./${TESTDIR}/listend2end.dat ${TESTDIR_TARGET}
    cp ./${TESTDIR}/backbone_idx.dat ${TESTDIR_TARGET}
    cd ${TESTDIR_TARGET}
    polymer_size -t ../../04-Melt_100ChPE500_TrappeUA/01-RESTART-0000-1000ns/traj_comp.xtc ../../04-Melt_100ChPE500_TrappeUA/02-RESTART-1000-2000ns/traj_comp.part0002.xtc ../../04-Melt_100ChPE500_TrappeUA/03-RESTART-2000-3000ns/traj_comp.part0003.xtc --topo ../../04-Melt_100ChPE500_TrappeUA/01-RESTART-0000-1000ns/topol.tpr --unwrap 1 --rg_massw --e2e listend2end.dat --e2acf --distributions --c2n backbone_idx.dat --bondorientation --isodf
    cd $WK
fi
# Example 04 ************************************

# Example 05 ************************************
if [[ $runjob -eq 0 || $runjob -eq 5 ]]; then
    echo "------------------------------" >>${WK}/out_tests.log
    echo "Example 05: Polyanagro: replicate_polymer Branch" >>${WK}/out_tests.log
    echo "------------------------------" >>${WK}/out_tests.log
    TESTDIR=05-Small_branchsystem
    TESTDIR_TARGET=${WK}/XX-CHECK_POLYANAGRO/${TESTDIR}
    mkdir -p ${TESTDIR_TARGET}
    cp ./${TESTDIR}/listend2end.dat ${TESTDIR_TARGET}
    cp ./${TESTDIR}/backbone_idx.dat ${TESTDIR_TARGET}
    cd ${TESTDIR_TARGET}
    polymer_size -t  ${DATAFOLDER}/0003Ch-C020-002br04/RUN-001/traj_comp.xtc  ${DATAFOLDER}/0003Ch-C020-002br04/RUN-002/traj_comp.part0002.xtc ${DATAFOLDER}/0003Ch-C020-002br04/RUN-003/traj_comp.part0003.xtc --topo ${DATAFOLDER}/0003Ch-C020-002br04/RUN-001/topol.tpr --unwrap 1 --rg_massw --unwrap 1 --rg_massw --e2e listend2end.dat --e2acf --distributions --c2n backbone_idx.dat --bondorientation --isodf 
    cd $WK
fi
# Example 05 ************************************

