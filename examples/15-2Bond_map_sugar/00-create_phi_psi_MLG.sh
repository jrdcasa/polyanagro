#!/bin/bash

monomers=24
#iatom=27 VMD
iatom=28 #GROMACS

echo "[ phi ]" > tmp_phi.dat
echo "[ psi ]" > tmp_psi.dat
echo "[ phi13 ]" >tmp_phi13.dat
echo "[ psi13 ]" >tmp_psi13.dat
echo "[ phi14 ]" >tmp_phi14.dat
echo "[ psi14 ]" >tmp_psi14.dat
echo "[ phipsi ]" >tmp_phipsi.dat
echo "[ allbbangles ]" >tmp_allbb.dat

for (( i=0; i<${monomers}-1; i++ )); do

    jatom=`echo ${iatom}-4 |bc -l`
    if [ $((i % 2)) -eq 0 ]; then
        katom=`echo ${jatom}-10 |bc -l`
    else
        katom=`echo ${jatom}-6 |bc -l`
    fi
    latom=`echo ${katom}-2 |bc -l`
    matom=`echo ${latom}-4 |bc -l`


    echo "===== $i ====="
    echo $iatom  $jatom  $katom  $latom >>tmp_phi.dat
    echo $jatom  $katom  $latom  $matom >>tmp_psi.dat
    if [ $((i % 2)) -eq 0 ]; then
        echo $iatom  $jatom  $katom  $latom >>tmp_phi13.dat
        echo $jatom  $katom  $latom  $matom >>tmp_psi13.dat
    else
        echo $iatom  $jatom  $katom  $latom >>tmp_phi14.dat
        echo $jatom  $katom  $latom  $matom >>tmp_psi14.dat
    fi
    echo $iatom  $jatom  $katom  $latom >>tmp_phipsi.dat
    echo $jatom  $katom  $latom  $matom >>tmp_phipsi.dat
    echo $iatom  $jatom  $katom  $latom >>tmp_allbb.dat
    echo $jatom  $katom  $latom  $matom >>tmp_allbb.dat

    iatom=`echo ${iatom}+21|bc -l`

done


echo "# Index start at 1 (GROMACS)" >dihedral_data_dist.ndx
cat tmp_phi.dat >>dihedral_data_dist.ndx
cat tmp_psi.dat >>dihedral_data_dist.ndx
cat tmp_phi13.dat >>dihedral_data_dist.ndx
cat tmp_psi13.dat >>dihedral_data_dist.ndx
cat tmp_phi14.dat >>dihedral_data_dist.ndx
cat tmp_psi14.dat >>dihedral_data_dist.ndx
cat tmp_phipsi.dat >>dihedral_data_dist.ndx
cat tmp_allbb.dat >>dihedral_data_dist.ndx

rm tmp_*dat

