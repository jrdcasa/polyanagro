#!bin/bash
nchains=25
natch=18
# First chain
idxs=1
idxe=10

echo "nmols: $nchains" >EVOH_headtail.dat
echo "#ich idx-head idx-tail (indexes start at 0)" >>EVOH_headtail.dat
for (( i=0; i<$nchains; i++ )) do
    echo "$i $idxs $idxe" >>EVOH_headtail.dat
    idxs=`echo $idxs+$natch|bc -l`
    idxe=`echo $idxe+$natch|bc -l`
done

