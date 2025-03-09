#!/bin/bash
natoms=12060

echo "[ setA ]" >index_full.ndx
for (( i=1; $i<=$natoms; i++ )); do
   echo $i >>index_full.ndx
done
echo "[ setB ]" >>index_full.ndx
for (( i=1; $i<=$natoms; i++ )); do
   echo $i >>index_full.ndx
done
