#!bin/bash

natoms=72
natch=24
nchains=3

rm -rf index_bychain.ndx
for ((ich=1; ich<=$nchains; ich++)); do
    echo "[chain$ich]" >>index_bychain.ndx

    for ((i=1; i<=$natch; i++)); do
       j=`echo $i $natch $ich | awk '{print $1+$2*($3-1)}'`
       echo -n "$j " >>index_bychain.ndx
    done
    echo "" >>index_bychain.ndx
done
echo "" >>index_bychain.ndx

