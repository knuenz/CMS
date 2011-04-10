#!/bin/sh

#settings:
iterations=2
TREENAME=INSTRUCT

make

cp GenScript GenScript${TREENAME}

for scenario in 4;do

for rap_ in 1;do
for pT_ in 6;do

if [ $pT_ -eq 8 ] && [ $rap_ -eq 1 ]; then
continue
fi

./GenScript${TREENAME} ${rap_}rap ${pT_}pT ${iterations}iter ${scenario}scen --JOBNAME=${TREENAME}

done
done

done

rm GenScript${TREENAME}
