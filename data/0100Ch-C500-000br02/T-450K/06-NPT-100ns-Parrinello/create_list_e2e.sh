#!/bin/bash

nchains=100
nbb=500
idx=0
echo "# ich head tail" > listend2end.dat
for (( ich=0; ich<nchains; ich++ )); do

    idxend=`echo $idx $nbb|awk '{print $1+($2-1)}'`
    echo $ich $idx $idxend >>listend2end.dat
    idx=`echo $idxend |awk '{print $1+1}'`
    
done
#backbone_idx.dat
rm -f backbone_idx.dat
idx=0
for (( ich=0; ich<nchains; ich++ )); do
    echo "[mol$ich]" >>backbone_idx.dat
    idxend=`echo $idx $nbb|awk '{print $1+($2-1)}'`
    for (( iat=$idx; iat<=$idxend; iat++ )); do 
        echo -n "$iat " >>backbone_idx.dat
    done
    echo "" >>backbone_idx.dat
    idx=`echo $idxend |awk '{print $1+1}'`
done
