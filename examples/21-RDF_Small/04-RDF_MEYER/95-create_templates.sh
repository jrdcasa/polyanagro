#!/bin/bash
NATOMS=170
echo ${NATOMS} >template_all.dat
for (( i=1; i<=$NATOMS; i++ )); do
    echo $i >>template_all.dat
done

