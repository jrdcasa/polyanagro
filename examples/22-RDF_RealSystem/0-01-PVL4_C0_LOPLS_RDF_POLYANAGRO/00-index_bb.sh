#!/bin/bash
natoms=12060

echo "[ setA ]" >index_bb.ndx
egrep "          0" confout_unwrap_labeled.psf |awk '{print $1}' >>index_bb.ndx
echo "[ setB ]" >>index_bb.ndx
egrep "          0" confout_unwrap_labeled.psf |awk '{print $1}' >>index_bb.ndx

echo "[ setA ]" >index_br.ndx
egrep "          1" confout_unwrap_labeled.psf |awk '{print $1}' >>index_br.ndx
echo "[ setB ]" >>index_br.ndx
egrep "          1" confout_unwrap_labeled.psf |awk '{print $1}' >>index_br.ndx
