#!/bin/sh

#settings:
iterations=200
TREENAME=INSTRUCT

make

cp GenScript GenScript${TREENAME}

for scenario in 2 3 4 5;do

for rap_ in 1 2;do
for pT_ in 6 7 8;do

if [ $pT_ -eq 8 ] && [ $rap_ -eq 1 ]; then
continue
fi

./GenScript${TREENAME} ${rap_}rap ${pT_}pT ${iterations}iter ${scenario}scen --JOBNAME=${TREENAME}

done
done

done

rm GenScript${TREENAME}
