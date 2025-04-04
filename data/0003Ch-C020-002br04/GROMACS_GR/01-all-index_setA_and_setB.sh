#!bin/bash

natoms=72

echo "[setA]" >index_all.ndx
for ((i=1; i<=$natoms; i++)); do
    echo -n "$i " >>index_all.ndx
done
echo "" >>index_all.ndx

echo "[setB]" >>index_all.ndx
for ((i=1; i<=$natoms; i++)); do
    echo -n "$i " >>index_all.ndx
done
