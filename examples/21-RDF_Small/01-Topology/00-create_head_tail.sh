#!bin/bash
nchains=10
nchA=10
natchA=17
# First kind chain
idxs=15
idxe=16
HEADFILENAME=headtail.dat

echo "nmols: $nchains" >${HEADFILENAME}
echo "#ich idx-head idx-tail (indexes start at 0)" >>${HEADFILENAME}
for (( i=0; i<$nchA; i++ )) do
    echo "$i $idxs $idxe" >>${HEADFILENAME}
    idxs=`echo $idxs+$natchA|bc -l`
    idxe=`echo $idxe+$natchA|bc -l`
done


