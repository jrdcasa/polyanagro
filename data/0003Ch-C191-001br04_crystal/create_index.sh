#/bin/bash

nend=585
rm -f index.dat
for (( i=0; i<$nend; i++ )); do
   echo -n "$i " >>index.dat
   echo -n "$i " 
done
