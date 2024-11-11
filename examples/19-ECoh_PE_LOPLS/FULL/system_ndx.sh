#!/bin/bash

natoms=18060
echo "[system]" >>index.ndx
for (( i=1; $i<=$natoms; i++ )); do
    echo $i >>index.ndx
done

natoms=18060
echo "[other]" >>index.ndx
for (( i=603; $i<=$natoms; i++ )); do
    echo $i >>index.ndx
done

