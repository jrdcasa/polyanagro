#!/bin/bash
# FOR GROMACS indexes start at 1
# FOR VMD     indexes start at 0
# A pdb containing only the backbone original indices. This can be produce
# with egrep from a labelled pdb with Topology
# ls ../../01-PVL4_C0_LOPLS/04-PDB_LABELED/confout_unwrap_labeled.pdb
# setA --> C carbonyl
# setB --> O carbonyl
natch=603
nchains=20
deltaA=15
iniA=3

echo "[ setA ]" >index.ndx
for ((ich=0; ich<$nchains; ich++)); do
    # Calculate the global offset for the current chain
    offset=$((ich * natch))
    # Generate indices for the current chain
    # natch-1 to skip the last monomer which is a terminal
    for ((j=iniA; j<=natch-1; j+=deltaA)); do
        iat_global=$((offset + j))
        echo ${iat_global} 
        echo ${iat_global} >>index.ndx
    done
done

deltaB=15
iniB=4
echo "[ setB ]" >>index.ndx
for ((ich=0; ich<$nchains; ich++)); do
    # Calculate the global offset for the current chain
    offset=$((ich * natch))
    # Generate indices for the current chain
    # natch-1 to skip the last monomer which is a terminal
    for ((j=iniB; j<=natch-1; j+=deltaB)); do
        iat_global=$((offset + j))
        echo ${iat_global} 
        echo ${iat_global} >>index.ndx
    done
done

