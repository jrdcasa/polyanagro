#/bin/bash

NATBB=191
NBRANCHPERCHAIN=1
LENGTHBRANCH=1
NATCH=192
NCHAINS=500
BPLOCAL=(97)
BPGLOBAL=(97)

rm -f template.dat

INDEX=0

for ((ichain=1;ichain<=$NCHAINS;ichain++));do

   if [ $NBRANCHPERCHAIN == 0 ]; then  
      echo "ich: "$ichain $NATBB "L" 0 >>template.dat
   else
      echo "ich: "$ichain $NATBB "R" $NBRANCHPERCHAIN >>template.dat
   fi

   for ((ibr=0;ibr<$NBRANCHPERCHAIN;ibr++)); do
      BPG=`echo ${BPGLOBAL[$ibr]}+$INDEX | bc -l`
      echo ${BPG} $LENGTHBRANCH >>template.dat
   done

   INDEX=`echo "$INDEX+$NATCH" | bc -l`

done
